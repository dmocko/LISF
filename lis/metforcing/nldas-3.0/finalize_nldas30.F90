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
! !ROUTINE: finalize_nldas30
! \label{finalize_nldas30}
!
! !REVISION HISTORY:
! 16 May 2025: David Mocko, Initial Specification
!                           (derived from finalize_merra2.F90)
!
! !INTERFACE:
subroutine finalize_nldas30(findex)

! !USES:
   use LIS_coreMod,     only    : LIS_rc
   use nldas30_forcingMod, only : nldas30_struc
!
! !DESCRIPTION:
!  Routine to cleanup NLDAS-3 forcing related memory allocations.
!
!EOP
   implicit none

   integer :: findex
   integer :: n

   do n = 1,LIS_rc%nnest
      select case(LIS_rc%met_interp(findex))

      case ("bilinear")
         deallocate(nldas30_struc(n)%n111)
         deallocate(nldas30_struc(n)%n121)
         deallocate(nldas30_struc(n)%n211)
         deallocate(nldas30_struc(n)%n221)
         deallocate(nldas30_struc(n)%w111)
         deallocate(nldas30_struc(n)%w121)
         deallocate(nldas30_struc(n)%w211)
         deallocate(nldas30_struc(n)%w221)

      case ("budget-bilinear")
         deallocate(nldas30_struc(n)%n111)
         deallocate(nldas30_struc(n)%n121)
         deallocate(nldas30_struc(n)%n211)
         deallocate(nldas30_struc(n)%n221)
         deallocate(nldas30_struc(n)%w111)
         deallocate(nldas30_struc(n)%w121)
         deallocate(nldas30_struc(n)%w211)
         deallocate(nldas30_struc(n)%w221)
         deallocate(nldas30_struc(n)%n112)
         deallocate(nldas30_struc(n)%n122)
         deallocate(nldas30_struc(n)%n212)
         deallocate(nldas30_struc(n)%n222)
         deallocate(nldas30_struc(n)%w112)
         deallocate(nldas30_struc(n)%w122)
         deallocate(nldas30_struc(n)%w212)
         deallocate(nldas30_struc(n)%w222)

      case ("neighbor")
         deallocate(nldas30_struc(n)%n113)

      case ("average")
         deallocate(nldas30_struc(n)%n111)
      end select
   enddo
   deallocate(nldas30_struc)

end subroutine finalize_nldas30
