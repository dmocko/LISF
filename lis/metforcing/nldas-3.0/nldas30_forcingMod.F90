!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module nldas30_forcingMod
!BOP
! !MODULE: nldas30_forcingMod
!
! !REVISION HISTORY:
! 16 May 2025: David Mocko, Initial Specification
!                           (derived from merra2_forcingMod.F90)
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the NLDAS-3 forcing data.
!  The data is global 1 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt nldas30\_struc}
!  includes the variables that specify the runtime options, and the
!  weights and neighbor information to be used for spatial interpolation.
!  They are described below:
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the ECMWF data
!  \item[nldas30time1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[nldas30time2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[nldas30dir]
!    Directory containing the input data
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for bilinear interpolation.
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[n113]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for n. neighbor interpolation.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for
!   temporal interpolation.
!  \end{description}
!
! !USES:
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_nldas30      !defines the native resolution of the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: nldas30_struc

!EOP
  type, public ::  nldas30_type_dec
     real         :: ts
     integer      :: ncold,nrold
     character(len=LIS_CONST_PATH_LEN) :: nldas30dir   ! NLDAS-3 Forcing Directory
     real*8       :: nldas30time1,nldas30time2
     logical      :: reset_flag

     integer                :: mi
     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)

     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)
     integer, allocatable   :: n113(:)
     integer                :: findtime1,findtime2
     logical                :: startFlag,dayFlag
     real, allocatable      :: nldasforc1(:,:,:,:),nldasforc2(:,:,:,:)

     integer           :: nvars
     real*8            :: ringtime
     integer           :: nIter,st_iterid,en_iterid

     real, allocatable :: metdata1(:,:,:)
     real, allocatable :: metdata2(:,:,:)
  end type nldas30_type_dec

  type(nldas30_type_dec), allocatable :: nldas30_struc(:)

contains

!BOP
!
! !ROUTINE: init_nldas30
! \label{init_nldas30}
!
! !REVISION HISTORY:
! 16 May 2025: David Mocko, Initial Specification
!                           (derived from init_merra2.F90)
!
! !INTERFACE:
subroutine init_nldas30(findex)

! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod

  implicit none
! !AGRUMENTS:
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for NLDAS-3
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_nldas30](\ref{readcrd_nldas30}) \newline
!     reads the runtime options specified for NLDAS-3 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
  real :: gridDesci(LIS_rc%nnest,50)
  integer :: n

  allocate(nldas30_struc(LIS_rc%nnest))

  do n = 1,LIS_rc%nnest
     nldas30_struc(n)%ncold = 11700
     nldas30_struc(n)%nrold = 6500
  enddo

  call readcrd_nldas30
  LIS_rc%met_nf(findex) = 8

  nldas30_struc%reset_flag = .false.

  do n = 1,LIS_rc%nnest
     nldas30_struc(n)%ts = 3600
     call LIS_update_timestep(LIS_rc,n,nldas30_struc(n)%ts)
  enddo

  gridDesci = 0

  do n = 1,LIS_rc%nnest
     gridDesci(n,1)  = 0
     gridDesci(n,2)  = nldas30_struc(n)%ncold
     gridDesci(n,3)  = nldas30_struc(n)%nrold
     gridDesci(n,4)  =  -52.005
     gridDesci(n,5)  = -168.995
     gridDesci(n,6)  = 128
     gridDesci(n,7)  =    7.005
     gridDesci(n,8)  =   71.995
     gridDesci(n,9)  =    0.01
     gridDesci(n,10) =    0.01
     gridDesci(n,20) = 0

     nldas30_struc(n)%mi = nldas30_struc(n)%ncold*nldas30_struc(n)%nrold

! Setting up weights for Interpolation
     if (trim(LIS_rc%met_interp(findex)).eq."bilinear") then
         allocate(nldas30_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         call bilinear_interp_input(n,gridDesci(n,:),                  &
                       nldas30_struc(n)%n111,nldas30_struc(n)%n121,    &
                       nldas30_struc(n)%n211,nldas30_struc(n)%n221,    &
                       nldas30_struc(n)%w111,nldas30_struc(n)%w121,    &
                       nldas30_struc(n)%w211,nldas30_struc(n)%w221)

     elseif (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
        allocate(nldas30_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(nldas30_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         call bilinear_interp_input(n,gridDesci(n,:),                  &
                       nldas30_struc(n)%n111,nldas30_struc(n)%n121,    &
                       nldas30_struc(n)%n211,nldas30_struc(n)%n221,    &
                       nldas30_struc(n)%w111,nldas30_struc(n)%w121,    &
                       nldas30_struc(n)%w211,nldas30_struc(n)%w221)

         allocate(nldas30_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(nldas30_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(nldas30_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(nldas30_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(nldas30_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(nldas30_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(nldas30_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(nldas30_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         call conserv_interp_input(n,gridDesci(n,:),                   &
                      nldas30_struc(n)%n112,nldas30_struc(n)%n122,     &
                      nldas30_struc(n)%n212,nldas30_struc(n)%n222,     &
                      nldas30_struc(n)%w112,nldas30_struc(n)%w122,     &
                      nldas30_struc(n)%w212,nldas30_struc(n)%w222)

     elseif (trim(LIS_rc%met_interp(findex)).eq."neighbor") then
        allocate(nldas30_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
        call neighbor_interp_input(n,gridDesci(n,:),nldas30_struc(n)%n113)

     else
        write(LIS_logunit,*) '[ERR] Interpolation option '//           &
                              trim(LIS_rc%met_interp(findex))//        &
                             ' for NLDAS-3 forcing is not supported'
        call LIS_endrun
     endif

     call LIS_registerAlarm("NLDAS-3 forcing alarm",86400.0,86400.0)
     nldas30_struc(n)%startFlag = .true.
     nldas30_struc(n)%dayFlag = .true.

     nldas30_struc(n)%nvars = 8

     allocate(nldas30_struc(n)%nldasforc1(1,nldas30_struc(n)%nvars,24, &
                                          LIS_rc%lnc(n)*LIS_rc%lnr(n)))
     allocate(nldas30_struc(n)%nldasforc2(1,nldas30_struc(n)%nvars,24, &
                                          LIS_rc%lnc(n)*LIS_rc%lnr(n)))

     nldas30_struc(n)%st_iterid = 1
     nldas30_struc(n)%en_iterId = 1
     nldas30_struc(n)%nIter = 1

     allocate(nldas30_struc(n)%metdata1(1,LIS_rc%met_nf(findex),       &
                                        LIS_rc%ngrid(n)))
     allocate(nldas30_struc(n)%metdata2(1,LIS_rc%met_nf(findex),       &
                                        LIS_rc%ngrid(n)))

     nldas30_struc(n)%metdata1 = 0
     nldas30_struc(n)%metdata2 = 0

     nldas30_struc(n)%nldasforc1 = LIS_rc%udef
     nldas30_struc(n)%nldasforc2 = LIS_rc%udef

  enddo   ! End nest loop

end subroutine init_nldas30

end module nldas30_forcingMod

