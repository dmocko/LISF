!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: AC72_coldstart
! \label{AC72_coldstart}
!
! !REVISION HISTORY: 
!   04 NOV 2024, Louise Busschaert; initial implementation for AC72
!
! !INTERFACE:
subroutine AC72_coldstart(mtype)
  ! !USES:
  use AC72_lsmMod
  use ac_global, only:    GetCompartment_Layer, &
       GetCompartment_theta
  use LIS_coreMod, only: LIS_rc
  use LIS_logMod, only: LIS_logunit
  use LIS_timeMgrMod, only: LIS_date2time
  !
  ! !DESCRIPTION:
  !
  !  This routine initializes the AC72 state variables with predefined
  !  values over the entire domain.
  !  For now, only soil moisture can be predefined. For future runs with
  !  salinity, this variable might also be considered.
  !
  !EOP

  implicit none

  integer, intent(in) :: mtype

  integer :: t, l, n

  do n=1, LIS_rc%nnest
     if (trim(LIS_rc%startcode) .eq. "coldstart") then
        write(LIS_logunit,*) "MSG: AC72_coldstart -- cold-starting AC72"
        do t=1, LIS_rc%npatch(n,mtype)
           do l=1, AC72_struc(n)%ac72(t)%NrCompartments
              AC72_struc(n)%ac72(t)%smc(l) = AC72_struc(n)%init_smc(GetCompartment_Layer(l))
           enddo
        enddo
     endif

     LIS_rc%yr = LIS_rc%syr
     LIS_rc%mo = LIS_rc%smo
     LIS_rc%da = LIS_rc%sda
     LIS_rc%hr = LIS_rc%shr
     LIS_rc%mn = LIS_rc%smn
     LIS_rc%ss = LIS_rc%sss

     call LIS_date2time(LIS_rc%time, LIS_rc%doy, LIS_rc%gmt, LIS_rc%yr,      &
          LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
     write(LIS_logunit,*) "MSG: AC72_coldstart -- ",     &
          "Using the specified start time ", LIS_rc%time
  enddo
end subroutine AC72_coldstart
