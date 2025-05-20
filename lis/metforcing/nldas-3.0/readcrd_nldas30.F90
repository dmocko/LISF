!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_nldas30
! \label{readcrd_nldas30}
!
! !REVISION HISTORY:
! 16 May 2025: David Mocko, Initial Specification
!                           (derived from readcrd_merra2.F90)
!
! !INTERFACE:
subroutine readcrd_nldas30
! !USES:
  use ESMF
  use LIS_logMod
  use LIS_coreMod, only        : LIS_rc, LIS_config
  use nldas30_forcingMod, only : nldas30_struc
!
! !DESCRIPTION:
!  This routine reads the options specific to NLDAS-3 forcing
!  from the LIS configuration file.
!
!EOP
  implicit none

  integer :: n,t,rc

  call ESMF_ConfigFindLabel(LIS_config,                                &
                            "NLDAS-3 forcing directory:",rc=rc)
  do n = 1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                          &
                                  nldas30_struc(n)%nldas30dir,rc=rc)
     call LIS_verify(rc,'NLDAS-3 forcing directory: not defined')
  enddo

  do n = 1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using NLDAS-3 forcing'
     write(LIS_logunit,*) '[INFO] NLDAS-3 forcing directory: ',&
                          trim(nldas30_struc(n)%nldas30dir)

     nldas30_struc(n)%nldas30time1 = 3000.0
     nldas30_struc(n)%nldas30time2 = 0.0
  enddo

end subroutine readcrd_nldas30

