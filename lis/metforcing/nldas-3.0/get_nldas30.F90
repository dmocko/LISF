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
! !ROUTINE: get_nldas30
! \label{get_nldas30}
!
! !REVISION HISTORY:
! 16 May 2025: David Mocko, Initial Specification
!                           (derived from get_merra2.F90)
!
! !INTERFACE:
subroutine get_nldas30(n,findex)
! !USES:
   use LIS_coreMod
   use LIS_timeMgrMod
   use LIS_logMod
   use LIS_metforcingMod
   use nldas30_forcingMod
   use LIS_constantsMod, only : LIS_CONST_PATH_LEN

   implicit none

! !ARGUMENTS:
   integer, intent(in) :: n
   integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates 1-hourly NLDAS-3 forcing.
!
!  The NLDAS-3 forcing data are organized into daily files, where
!  each daily file contains 24 one-hourly records of forcing fields.
!
!  In general, metforcing readers read the forcing data before the current
!  time, referred to as bookend1, and after the current time, referred to as
!  bookend2.  Then the readers temporally interpolate between bookend1 and
!  bookend2.  Here, each bookend contains 24 one-hourly records of forcing,
!  and, in general, the NLDAS-3 reader will be temporally interpolating
!  from one hour-interval to the next hour-interval, where both hour-intervals
!  are contained in bookend1.  Issues arise between the hours 23z of one day
!  and 1z of the next day, which are complicated by the size of LIS' running
!  time-step, say 15mn, 30mn, or 1hr.
!
!  Below are some examples to illustrate the timing logic of the NLDAS-3 reader.
!
!  \begin{verbatim}
!          ---*---|---*---|---*---|---*---|---*---|---*---|---*---|---*---
!  hour          21      22      23       0       1       2       3
!  hr_int         <---22--X--23---X--24---X---1---X---2---X---3--->
!
!
!  where:
!  hour is the hour UTC
!  hr_int is the hour-interval
!  * marks the valid point for the interval <--- hr_int --->
!
!  For example, interval 2, is from 1z to 2z, valid at 01:30z.
!  \end{verbatim}
!
!  First, consider the situation where the start time is
!  2005-11-01T21:00:00 and the time-step is 15mn.  Here bookend1
!  contains 01 Nov data and bookend2 contains 02 Nov data.
!
!  \begin{tabular}{lll}
!  time-step & time     & intervals                         \cr
!  1         & 21:15:00 & 21 of bookend1 and 22 of bookend1 \cr
!  2         & 21:30:00 & 22 of bookend1 and 23 of bookend1 \cr
!  3         & 21:45:00 & 22 of bookend1 and 23 of bookend1 \cr
!  4         & 22:00:00 & 22 of bookend1 and 23 of bookend1 \cr
!  5         & 22:15:00 & 22 of bookend1 and 23 of bookend1 \cr
!  6         & 22:30:00 & 23 of bookend1 and 24 of bookend1 \cr
!  7         & 22:45:00 & 23 of bookend1 and 24 of bookend1 \cr
!  8         & 23:00:00 & 23 of bookend1 and 24 of bookend1 \cr
!  9         & 23:15:00 & 23 of bookend1 and 24 of bookend1 \cr
!  10        & 23:30:00 & 24 of bookend1 and  1 of bookend2 \cr
!  11        & 23:45:00 & 24 of bookend1 and  1 of bookend2 \cr
!  12        & 00:00:00 & 24 of bookend1 and  1 of bookend2 \cr
!  13        & 00:15:00 & 24 of bookend1 and  1 of bookend2
!  \end{tabular}
!
!  At 00:30:00, the NLDAS-3 reader moves bookend2 to bookend1,
!  and reads 03 Nov as bookend2.
!
!  \begin{tabular}{lll}
!  14        & 00:30:00 & 1 of bookend1 and  2 of bookend1 \cr
!  15        & 00:45:00 & 1 of bookend1 and  2 of bookend1 \cr
!  16        & 01:00:00 & 1 of bookend1 and  2 of bookend1 \cr
!  17        & 01:15:00 & 1 of bookend1 and  2 of bookend1 \cr
!  18        & 01:30:00 & 2 of bookend1 and  3 of bookend1
!  \end{tabular}
!
!  Next, consider a similar situation where the start time is
!  2005-11-01T21:00:00 and the time-step is 30mn.  Here bookend1
!  contains 01 Nov data and bookend2 contains 02 Nov data.
!
!  \begin{tabular}{lll}
!  time-step & time     & intervals                         \cr
!  1         & 21:30:00 & 22 of bookend1 and 23 of bookend1 \cr
!  2         & 22:00:00 & 22 of bookend1 and 23 of bookend1 \cr
!  3         & 22:30:00 & 23 of bookend1 and 24 of bookend1 \cr
!  4         & 23:00:00 & 23 of bookend1 and 24 of bookend1 \cr
!  5         & 23:30:00 & 24 of bookend1 and  1 of bookend2 \cr
!  6         & 00:00:00 & 24 of bookend1 and  1 of bookend2
!  \end{tabular}
!
!  At 00:30:00, the NLDAS-3 reader moves bookend2 to bookend1,
!  and reads 03 Nov as bookend2.
!
!  \begin{tabular}{lll}
!  7         & 00:30:00 &  1 of bookend1 and  2 of bookend1 \cr
!  8         & 01:00:00 &  1 of bookend1 and  2 of bookend1 \cr
!  9         & 01 30:00 &  2 of bookend1 and  3 of bookend1
!  \end{tabular}
!
!  Finally, consider the situation where the start time is
!  2005-11-01T21:00:00 and the time-step is 1hr.  Here bookend1
!  contains 01 Nov data and bookend2 contains 02 Nov data.
!
!  \begin{tabular}{lll}
!  time-step & time     & intervals                         \cr
!  1         & 22:00:00 & 22 of bookend1 and 23 of bookend1 \cr
!  2         & 23:00:00 & 23 of bookend1 and 24 of bookend1 \cr
!  3         & 00:00:00 & 24 of bookend1 and  1 of bookend2
!  \end{tabular}
!
!  At 01:00:00, the NLDAS-3 reader moves bookend2 to bookend1,
!  and reads 03 Nov as bookend2.
!
!  \begin{tabular}{lll}
!  4         & 01:00:00 &  1 of bookend1 and  2 of bookend1 \cr
!  5         & 02:00:00 &  2 of bookend1 and  3 of bookend1
!  \end{tabular}
!
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    forcing dataset index
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    call to advance or retract time
!  \item[nldas30files](\ref{nldas30files}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_nldas30](\ref{read_nldas30}) \newline
!    call to read the NLDAS-3 data and perform spatial interpolation
!  \end{description}
!EOP
   integer           :: order
   integer           :: ferror
   character(len=LIS_CONST_PATH_LEN) :: filename
   integer           :: c,r,kk
   integer           :: yr1,mo1,da1,hr1,mn1,ss1,doy1
   integer           :: yr2,mo2,da2,hr2,mn2,ss2,doy2
   real*8            :: time1,time2,timenow
   real              :: gmt1,gmt2,ts1,ts2
   integer           :: gid

   integer           :: hr_int1,hr_int2
   integer           :: movetime  ! Flag to move bookend2 files to bookend1

! _________________________________________________________

   ! Please note that the timing logic has been tested only for
   ! these scenarios:
   !
   ! startime of 2005-11-01T00:30:00 with time-step of 15mn
   ! startime of 2005-11-01T00:00:00 with time-step of 30mn
   ! startime of 2005-11-01T00:00:00 with time-step of 1hr
   ! startime of 2005-11-02T00:00:00 with time-step of 15mn

   if (LIS_rc%nts(n).gt.3600) then   ! > 1-hr timestep
      write(LIS_logunit,*)                                              &
         '[ERR] When running LIS with NLDAS-3, the clock'
      write(LIS_logunit,*)                                              &
         '[ERR] should run with a timestep less than or'
      write(LIS_logunit,*) '[ERR] equal to one hour.'
      call LIS_endrun
   endif

   nldas30_struc(n)%findtime1 = 0
   nldas30_struc(n)%findtime2 = 0
   movetime = 0

   ! Initialize ts1 and ts2 timepoints at beginning of run:
   if ((LIS_get_nstep(LIS_rc,n).eq.1).or.(LIS_rc%rstflag(n).eq.1).or.   &
      nldas30_struc(n)%reset_flag) then
      nldas30_struc(n)%findtime1 = 1
      nldas30_struc(n)%findtime2 = 1
      LIS_rc%rstflag(n) = 0
      nldas30_struc(n)%reset_flag = .false.

      yr1 = LIS_rc%yr
      mo1 = LIS_rc%mo
      da1 = LIS_rc%da
      hr1 = 0
      mn1 = 0
      ss1 = 0
      ! initialize ringtime to tomorrow at 00:00z
      ts1 = 86400 ! 1 day
      call LIS_tick(nldas30_struc(n)%ringtime,doy1,gmt1,                &
         yr1,mo1,da1,hr1,mn1,ss1,ts1)
   endif

   !----------------------------------------------------------
   ! Determine current time
   !----------------------------------------------------------
   yr1 = LIS_rc%yr
   mo1 = LIS_rc%mo
   da1 = LIS_rc%da
   hr1 = LIS_rc%hr
   mn1 = LIS_rc%mn
   ss1 = 0
   call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,0.0)

   ! If current time >= time for when to move bookend values, 2 --> 1.
   if (timenow.ge.nldas30_struc(n)%ringtime) then
      nldas30_struc(n)%findtime2 = 1
      if (nldas30_struc(n)%findtime1.eq.0) then
         movetime = 1
      endif

      ! reset ringtime to tomorrow at 00:00z
      yr1 = LIS_rc%yr
      mo1 = LIS_rc%mo
      da1 = LIS_rc%da
      hr1 = 0
      mn1 = 0
      ss1 = 0
      ts1 = 86400 ! 1 day
      call LIS_tick(nldas30_struc(n)%ringtime,doy1,gmt1,                &
         yr1,mo1,da1,hr1,mn1,ss1,ts1)
   endif

   if (nldas30_struc(n)%findtime1.eq.1) then
      !----------------------------------------------------------
      ! Determine NLDAS-3 Forcing 1 Time
      !----------------------------------------------------------
      yr1 = LIS_rc%yr
      mo1 = LIS_rc%mo
      da1 = LIS_rc%da
      hr1 = LIS_rc%hr
      mn1 = LIS_rc%mn
      ss1 = 0
      ts1 = 0
      call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
   endif

   if (nldas30_struc(n)%findtime2.eq.1) then
      !----------------------------------------------------------
      ! Determine NLDAS-3 Forcing 2 Time
      !----------------------------------------------------------
      yr2 = LIS_rc%yr
      mo2 = LIS_rc%mo
      da2 = LIS_rc%da
      hr2 = LIS_rc%hr
      mn2 = 0
      ss2 = 0

      if (nldas30_struc(n)%findtime1.eq.1) then
         ! need tomorrow for bookend2
         ts2 = 86400
      endif

      call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
   endif

   ! Read NLDAS-3 - Bookend 1 files:
   if (nldas30_struc(n)%findtime1.eq.1) then
      order = 1
      do kk = nldas30_struc(n)%st_iterid,nldas30_struc(n)%en_iterid
         call nldas30files(nldas30_struc(n)%nldas30dir,            &
            yr1,mo1,da1,filename)

         call read_nldas30(n,order,mo1,findex,filename,                 &
            nldas30_struc(n)%nldasforc1(kk,:,:,:),ferror)
      enddo
   endif

   ! Read NLDAS-3 - Bookend 2 files (next day - store values):
   if (nldas30_struc(n)%findtime2.eq.1) then
      ! Move bookend 2 files to bookend 1 timeframe:
      if (movetime.eq.1) then
         nldas30_struc(n)%nldasforc1 = nldas30_struc(n)%nldasforc2
         nldas30_struc(n)%nldasforc2 = LIS_rc%udef
      endif

      order = 2
      do kk = nldas30_struc(n)%st_iterid,nldas30_struc(n)%en_iterid
         call nldas30files(nldas30_struc(n)%nldas30dir,            &
            yr2,mo2,da2,filename)

         call read_nldas30(n,order,mo2,findex,filename,                 &
            nldas30_struc(n)%nldasforc2(kk,:,:,:),ferror)
      enddo
   endif

   if (timenow.ge.nldas30_struc(n)%nldas30time2) then
      yr1 = LIS_rc%yr
      mo1 = LIS_rc%mo
      da1 = LIS_rc%da
      hr1 = LIS_rc%hr
      mn1 = 0
      ss1 = 0

      yr2 = LIS_rc%yr
      mo2 = LIS_rc%mo
      da2 = LIS_rc%da
      hr2 = LIS_rc%hr
      mn2 = 0
      ss2 = 0

      ts1 = 0 + LIS_rc%ts
      ts2 = 60*60 + LIS_rc%ts

      call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
      call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

      if (LIS_rc%hr.eq.23) then
         hr_int1 = 24
         hr_int2 = 1
      ! For all other hours (0-22Z):
      else
         hr_int1 = hr1+1
         hr_int2 = hr2+1
      endif

      ! Assign NLDAS-3 forcing fields to two LIS time-interp placeholders (metdata1,2):
      do r = 1,LIS_rc%lnr(n)
         do c = 1,LIS_rc%lnc(n)
            if (LIS_domain(n)%gindex(c,r).ne.-1) then
               nldas30_struc(n)%metdata1(:,:,LIS_domain(n)%gindex(c,r)) = &
                  nldas30_struc(n)%nldasforc1(:,:,hr_int1,&  ! Store hour: Current hour (same day)
                  (c+(r-1)*LIS_rc%lnc(n)))
               nldas30_struc(n)%metdata2(:,:,LIS_domain(n)%gindex(c,r)) = &
                  nldas30_struc(n)%nldasforc1(:,:,hr_int2,&  ! Store hour:  next hour (same day)
                  (c+(r-1)*LIS_rc%lnc(n)))
            endif
         enddo
      enddo

      ! Assign the hourly times:
      nldas30_struc(n)%nldas30time2 = time2
      nldas30_struc(n)%nldas30time1 = time1
   endif
end subroutine get_nldas30

!BOP
! !ROUTINE: nldas30files
! \label{nldas30files}
!
! !INTERFACE:
subroutine nldas30files(nldas30dir,yr,mo,da,filename)

! !USES:
   !none

   implicit none
! !ARGUMENTS:
   character(len=*), intent(in)  :: nldas30dir
   integer, intent(in)           :: yr,mo,da
   character(len=*), intent(out) :: filename

! !DESCRIPTION:
!   This subroutine puts together the daily netCDF NLDAS-3 file name.
!
!  The arguments are:
!  \begin{description}
!  \item[nldas30dir]
!    Name of the NLDAS-3 directory
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[filename]
!   name of the timestamped daily file
!  \end{description}
!
!EOP

   character*4 :: cyear
   character*2 :: cmonth
   character*8 :: cdate
   character*7 :: dir

   write(unit=cyear, fmt='(i4.4)') yr
   write(unit=cmonth,fmt='(i2.2)') mo
   write(unit=cdate, fmt='(i4.4,i2.2,i2.2)') yr,mo,da

   dir = '/'//cyear//cmonth

   filename = trim(nldas30dir)//dir//'/NLDAS_FOR0010_H.A'//cdate//      &
      '.030.beta.nc'

end subroutine nldas30files
