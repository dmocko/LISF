!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: MODIS_lai_Mod
!
! !DESCRIPTION:
!   This module contains interfaces and subroutines to read and use
!        the MODIS MCD15A2H LAI data as a parameter input.
!
! !REVISION HISTORY:
!  16 Jun 2020    Wanshu Nie; initial specification
!   8 Nov 2021    David Mocko; Converted from DA obs to param input
!
module MODIS_lai_Mod
! !USES:
  use ESMF
  use map_utils

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: setup_MODIS_lai
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: MODIS_lai_struc
!EOP
  type, public:: MODIS_lai_dec

     character*100          :: odir
     logical                :: startMode
     integer                :: tsmooth
     integer                :: climofill
     integer                :: qcflag
     integer                :: nc, nr
     real                   :: gridDesci(50)
     real,    allocatable   :: rlat(:)
     real,    allocatable   :: rlon(:)
     integer, allocatable   :: n11(:)
     integer, allocatable   :: n12(:)
     integer, allocatable   :: n21(:)
     integer, allocatable   :: n22(:)
     real,    allocatable   :: w11(:)
     real,    allocatable   :: w12(:)
     real,    allocatable   :: w21(:)
     real,    allocatable   :: w22(:)
     real*8                 :: time1, time2
     type(proj_info)        :: mod_proj

  end type MODIS_lai_dec

  type(MODIS_lai_dec),allocatable :: MODIS_lai_struc(:)

contains

!BOP
!
! !ROUTINE: setup_MODIS_lai
! \label{setup_MODIS_lai}
!
! !INTERFACE:
  subroutine setup_MODIS_lai(n)
! !USES:
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_logmod

    implicit none

! !ARGUMENTS:
    integer                ::  n
!
! !DESCRIPTION:
!
!   This routine completes the runtime initializations and
!   creation of data structures required for handling MCD15A2H LAI data.
!
!   The arguments are:
!   \begin{description}
!    \item[k] number of observation state
!    \item[OBS\_State]   observation state
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  status
    real*8                 ::  cornerlat1, cornerlon1
    real*8                 ::  cornerlat2, cornerlon2

    if (.not.allocated(MODIS_lai_struc)) then
       allocate(MODIS_lai_struc(LIS_rc%nnest))
    endif

    call ESMF_ConfigFindLabel(LIS_config,"MCD15A2H LAI data directory:",&
         rc=status)
    call ESMF_ConfigGetAttribute(LIS_config,MODIS_lai_struc(n)%odir,&
         rc=status)
    call LIS_verify(status, 'MCD15A2H LAI data directory: is missing')

    call ESMF_ConfigFindLabel(LIS_config,&
       "MCD15A2H LAI apply temporal smoother between 8-day intervals:",&
       rc=status)
    call ESMF_ConfigGetAttribute(LIS_config,MODIS_lai_struc(n)%tsmooth,&
         rc=status)
    call LIS_verify(status,&
       'MCD15A2H LAI apply temporal smoother between 8-day intervals:'//&
       'is missing')

    call ESMF_ConfigFindLabel(LIS_config,&
       "MCD15A2H LAI apply climatological fill values:",rc=status)
    call ESMF_ConfigGetAttribute(LIS_config,MODIS_lai_struc(n)%climofill,&
         rc=status)
    call LIS_verify(status,&
       'MCD15A2H LAI apply climatological fill values: is missing')

    call ESMF_ConfigFindLabel(LIS_config,"MCD15A2H LAI apply QC flags:",&
         rc=status)
    call ESMF_ConfigGetAttribute(LIS_config,MODIS_lai_struc(n)%qcflag,&
         rc=status)
    call LIS_verify(status, 'MCD15A2H LAI apply QC flags: is missing')

    cornerlat1 = max(-59.9978927, nint((LIS_rc%gridDesc(n,4)+59.9978927)/0.00416667)*0.00416667-59.9978927-50*0.00416667)
    cornerlon1 = max(-179.9979167, nint((LIS_rc%gridDesc(n,5)+179.9979167)/0.00416667)*0.00416667-179.9979167-50*0.00416667)
    cornerlat2 = min(89.9979167, nint((LIS_rc%gridDesc(n,7)+59.9978927)/0.00416667)*0.00416667-59.9978927+50*0.00416667)
    cornerlon2 = min(179.9979167, nint((LIS_rc%gridDesc(n,8)+179.9979167)/0.00416667)*0.00416667-179.9979167+50*0.00416667)

    MODIS_lai_struc(n)%nc = nint((cornerlon2-cornerlon1)/0.00416667)+1
    MODIS_lai_struc(n)%nr = nint((cornerlat2-cornerlat1)/0.00416667)+1

    MODIS_lai_struc(n)%gridDesci     = 0
    MODIS_lai_struc(n)%gridDesci(1)  = 0
    MODIS_lai_struc(n)%gridDesci(2)  = MODIS_lai_struc(n)%nc
    MODIS_lai_struc(n)%gridDesci(3)  = MODIS_lai_struc(n)%nr
    MODIS_lai_struc(n)%gridDesci(4)  = cornerlat1
    MODIS_lai_struc(n)%gridDesci(5)  = cornerlon1
    MODIS_lai_struc(n)%gridDesci(6)  = 128
    MODIS_lai_struc(n)%gridDesci(7)  = cornerlat2
    MODIS_lai_struc(n)%gridDesci(8)  = cornerlon2
    MODIS_lai_struc(n)%gridDesci(9)  = 0.00416667
    MODIS_lai_struc(n)%gridDesci(10) = 0.00416667
    MODIS_lai_struc(n)%gridDesci(20) = 64

    if (LIS_isAtAfinerResolution(n,MODIS_lai_struc(n)%gridDesci(9))) then
       allocate(MODIS_lai_struc(n)%rlat(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(MODIS_lai_struc(n)%rlon(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       allocate(MODIS_lai_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(MODIS_lai_struc(n)%n12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(MODIS_lai_struc(n)%n21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(MODIS_lai_struc(n)%n22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(MODIS_lai_struc(n)%w11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(MODIS_lai_struc(n)%w12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(MODIS_lai_struc(n)%w21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(MODIS_lai_struc(n)%w22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       call bilinear_interp_input(MODIS_lai_struc(n)%gridDesci,LIS_rc%gridDesc,&
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            MODIS_lai_struc(n)%rlat, MODIS_lai_struc(n)%rlon,&
            MODIS_lai_struc(n)%n11, MODIS_lai_struc(n)%n12, &
            MODIS_lai_struc(n)%n21, MODIS_lai_struc(n)%n22, &
            MODIS_lai_struc(n)%w11, MODIS_lai_struc(n)%w12, &
            MODIS_lai_struc(n)%w21, MODIS_lai_struc(n)%w22)
    else
       allocate(MODIS_lai_struc(n)%n11(&
            MODIS_lai_struc(n)%nc*MODIS_lai_struc(n)%nr))
       call upscaleByAveraging_input(MODIS_lai_struc(n)%gridDesci, LIS_rc%gridDesc,&
            MODIS_lai_struc(n)%nc*MODIS_lai_struc(n)%nr, &
            LIS_rc%lnc(n)*LIS_rc%lnr(n), MODIS_lai_struc(n)%n11)
    endif

    call LIS_registerAlarm("MODIS LAI read alarm", 86400.0, 86400.0)

    MODIS_lai_struc(n)%startMode = .true.

    call map_set(PROJ_LATLON, MODIS_lai_struc(n)%gridDesci(4), MODIS_lai_struc(n)%gridDesci(5), &
          0.0, MODIS_lai_struc(n)%gridDesci(9), MODIS_lai_struc(n)%gridDesci(10), 0.0, &
          MODIS_lai_struc(n)%nc, MODIS_lai_struc(n)%nr, MODIS_lai_struc(n)%mod_proj)

  end subroutine setup_MODIS_lai

end module MODIS_lai_Mod

