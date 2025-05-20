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
! !MODULE: reset_nldas30
! \label{reset_nldas30}
!
! !REVISION HISTORY:
! 16 May 2025: David Mocko, Initial Specification
!                           (derived from reset_merra2.F90)
!
! !INTERFACE:
subroutine reset_nldas30
! !USES:
  use LIS_coreMod,  only   : LIS_rc
  use LIS_timeMgrMod, only : LIS_date2time
  use nldas30_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for NLDAS-3 forcing.
!
!EOP
  implicit none
  integer :: n

  do n = 1,LIS_rc%nnest
     nldas30_struc(n)%startFlag    = .true.
     nldas30_struc(n)%dayFlag      = .true.
     nldas30_struc(n)%nldas30time1 = 3000.0
     nldas30_struc(n)%nldas30time2 = 0.0
     nldas30_struc(n)%ringtime     = 0.0
     nldas30_struc(n)%reset_flag   = .true.
  enddo

end subroutine reset_nldas30

