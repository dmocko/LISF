!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_MODIS_lai
! \label{read_MODIS_lai}
!
! !REVISION HISTORY:
!  16 Jun 2020    Wanshu Nie; initial specification
!
! !INTERFACE:
subroutine read_MODIS_lai(n, wt1, wt2, array1, array2)
! !USES:
  use ESMF
!  use LIS_mpiMod
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod,    only : LIS_logunit, LIS_verify
  use MODIS_lai_Mod, only : MODIS_lai_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  real                :: wt1
  real                :: wt2
  real, intent(inout) :: array1(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real, intent(inout) :: array2(LIS_rc%lnc(n)*LIS_rc%lnr(n))

!
! !DESCRIPTION:
!  reads the MCD15A2H LAI observations from NETCDF files.

!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest
!  \item[wt1] interpolation weight at time1
!  \item[wt2] interpolation weight at time2
!  \item[array1] LAI at time1
!  \item[array2] LAI at time2
!  \end{description}
!
!EOP
  character*100          :: fname1,fname2,climofile1,climofile2
  integer                :: cyr,cmo,cda,chr,cmn,css,cdoy
  real                   :: ts
  integer                :: t,c,r,count,found
  real                   :: cgmt
  real*8                 :: time
  logical                :: alarmCheck,file_exists
  real                   :: laiobs1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                   :: laiobs2(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                   :: timenow

  alarmCheck = LIS_isAlarmRinging(LIS_rc,"MODIS LAI read alarm")

  if (MODIS_lai_struc(n)%tsmooth.eq.1) then
     if (alarmCheck.or.MODIS_lai_struc(n)%startMode) then
        MODIS_lai_struc(n)%startMode = .false.

        cyr  = LIS_rc%yr
        cmo  = LIS_rc%mo
        cda  = LIS_rc%da
        cdoy = LIS_rc%doy
        chr  = 0
        cmn  = 0
        css  = 0
        ts   = -86400.0

        file_exists = .false.
        count = 0
        do while((.not.file_exists).and.(count.lt.8))
           call create_MODIS_lai_filename(MODIS_lai_struc(n)%odir, &
                cyr,cdoy,fname1,climofile1)

           inquire(file=fname1,exist=file_exists)
           if (file_exists) then
              call LIS_tick(MODIS_lai_struc(n)%time1,cdoy,cgmt,cyr,cmo,cda, &
                   chr,cmn,css,0.0)
              exit;
           else
              !go back a day till 8 days
              call LIS_tick(MODIS_lai_struc(n)%time1,cdoy,cgmt,cyr,cmo,cda, &
                   chr,cmn,css,ts)
              count = count + 1
           endif
        enddo

        if (count.ne.8) then
! Find the next data location
           if (cdoy.lt.361) then
              call LIS_tick(MODIS_lai_struc(n)%time2,cdoy,cgmt,cyr,cmo,cda, &
                   chr,cmn,css,(-8.0)*ts)
           else
              if (mod(cyr,4).ne.0) then
                 call LIS_tick(MODIS_lai_struc(n)%time2,cdoy,cgmt,cyr,cmo,cda, &
                      chr,cmn,css,(-5.0)*ts)
              else
                 call LIS_tick(MODIS_lai_struc(n)%time2,cdoy,cgmt,cyr,cmo,cda, &
                      chr,cmn,css,(-6.0)*ts)
              endif
           endif
           call create_MODIS_lai_filename(MODIS_lai_struc(n)%odir, &
                cyr,cdoy,fname2,climofile2)

           write(LIS_logunit,*) '[INFO] Reading ',trim(fname1)
           call read_MODIS_LAI_data(n,fname1,climofile1,laiobs1)

           write(LIS_logunit,*) '[INFO] Reading ',trim(fname2)
           call read_MODIS_LAI_data(n,fname2,climofile2,laiobs2)
           found = 1
        else
           found = 0
           laiobs1 = LIS_rc%udef
           laiobs2 = LIS_rc%udef
        endif
     endif

! Not smooth - just read in new data every 8 days
  else
     wt1 = 1.0
     wt2 = 0.0
     if (alarmCheck.or.MODIS_lai_struc(n)%startMode) then
        MODIS_lai_struc(n)%startMode = .false.

        call create_MODIS_lai_filename(MODIS_lai_struc(n)%odir, &
             LIS_rc%yr,LIS_rc%doy,fname1,climofile1)

        inquire(file=fname1,exist=file_exists)
        if (file_exists) then
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname1)
           call read_MODIS_LAI_data(n,fname1,climofile1,laiobs1)
        else
           write(LIS_logunit,*) '[WARN] Missing LAI file: ',trim(fname1)
        endif

        do t=1,LIS_rc%ntiles(n)
           c = LIS_domain(n)%tile(t)%col
           r = LIS_domain(n)%tile(t)%row
           array1(t) = laiobs1(c,r)
        enddo
        array2 = LIS_rc%udef
     else
        array1 = LIS_rc%udef
        array2 = LIS_rc%udef
     endif
  endif

  if (MODIS_lai_struc(n)%tsmooth.eq.1) then
! Interpolate between two data points every day
     timenow = float(LIS_rc%hr)*3600 + 60.0*LIS_rc%mn + LIS_rc%ss
     alarmCheck = (mod(timenow,691200.0).eq.0)

     call LIS_tick(time,cdoy,cgmt,LIS_rc%yr,LIS_rc%mo,LIS_rc%da, &
                                  LIS_rc%hr,LIS_rc%mn,LIS_rc%ss,0.0)
     wt2 = (time - MODIS_lai_struc(n)%time1) / &
           (MODIS_lai_struc(n)%time2-MODIS_lai_struc(n)%time1)
     wt1 = 1.0 - wt2

     if (alarmCheck) then
        if (found.eq.1) then

           do t=1,LIS_rc%ntiles(n)
              c = LIS_domain(n)%tile(t)%col
              r = LIS_domain(n)%tile(t)%row
              array1(t) = laiobs1(c,r)
           enddo

           do t=1,LIS_rc%ntiles(n)
              c = LIS_domain(n)%tile(t)%col
              r = LIS_domain(n)%tile(t)%row
              array2(t) = laiobs2(c,r)
           enddo
        else
           array1 = LIS_rc%udef
           array2 = LIS_rc%udef
        endif
     endif
  endif
!  deallocate(locallai)

end subroutine read_MODIS_lai

!BOP
!
! !ROUTINE: read_MODIS_LAI_data
! \label{read_MODIS_LAI_data}
!
! !INTERFACE:
subroutine read_MODIS_LAI_data(n, fname, climofile, laiobs_ip)
!
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use MODIS_lai_Mod, only : MODIS_lai_struc

  implicit none
!
! !INPUT PARAMETERS:
!
  integer                       :: n
  character(len=*)              :: fname
  character(len=*)              :: climofile
  real                          :: laiobs_ip(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real*8                        :: cornerlat(2),cornerlon(2)

! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine reads the MCD15A2H LAI file and applies the data
!  quality flags to filter the data.
!
!  The arguments are:
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the MCD15A2H LAI file
!  \item[climofile]    Generated MCD152AH LAI climatology file
!  \item[laiobs\_ip]   MCD15A2H LAI data processed to the LIS domain
!  \end{description}
!
!EOP

!--------------Wanshu -----------------------
  integer                 :: lat_off,lon_off
  integer                 :: lai(MODIS_lai_struc(n)%nc,MODIS_lai_struc(n)%nr)
  integer                 :: flag(MODIS_lai_struc(n)%nc,MODIS_lai_struc(n)%nr)
  real                    :: lai_flagged(MODIS_lai_struc(n)%nc,MODIS_lai_struc(n)%nr)
  real                    :: lai_in(MODIS_lai_struc(n)%nc*MODIS_lai_struc(n)%nr)
  logical*1               :: lai_data_b(MODIS_lai_struc(n)%nc*MODIS_lai_struc(n)%nr)
  logical*1               :: laiobs_b_ip(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: laiobs_climo_ip(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                 :: c,r,t
  integer                 :: nid,ios
  integer                 :: laiid,flagid

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LIS_verify(ios,'Error opening file '//trim(fname))

  ios = nf90_inq_varid(nid,'Lai_500m',laiid)
  call LIS_verify(ios,'Error nf90_inq_varid: Lai_500m')

  ios = nf90_inq_varid(nid,'FparLai_QC',flagid)
  call LIS_verify(ios,'Error nf90_inq_varid: flag')

  cornerlat(1) = MODIS_lai_struc(n)%gridDesci(4)
  cornerlon(1) = MODIS_lai_struc(n)%gridDesci(5)
  cornerlat(2) = MODIS_lai_struc(n)%gridDesci(7)
  cornerlon(2) = MODIS_lai_struc(n)%gridDesci(8)

  lai_data_b = .false.

  lat_off = nint((cornerlat(1) +  89.9979167)/0.00416667)+1
  lon_off = nint((cornerlon(1) + 179.9979167)/0.00416667)+1
      
  ios = nf90_get_var(nid, laiid, lai, &
       start=(/lon_off,lat_off/), &
       count=(/MODIS_lai_struc(n)%nc,MODIS_lai_struc(n)%nr/))
  call LIS_verify(ios, 'Error nf90_get_var: Lai_500m')

!  if ((cornerlon(1).lt.-120.0).and.(cornerlat(1).gt.35.0)) then
!      print *,' '
!      print *,'corners and nr/nc and lat_off/lon_off'
!      print *,cornerlat(1),cornerlat(2),cornerlon(1),cornerlon(2)
!      print *,MODIS_lai_struc(n)%nr,MODIS_lai_struc(n)%nc
!      print *,lat_off,lon_off
!  endif

  ios = nf90_get_var(nid, flagid, flag, &
       start=(/lon_off,lat_off/), &
       count=(/MODIS_lai_struc(n)%nc,MODIS_lai_struc(n)%nr/))
  call LIS_verify(ios, 'Error nf90_get_var: flag')

  ios = nf90_close(ncid=nid)
  call LIS_verify(ios,'Error closing file '//trim(fname))

  do r = 1,MODIS_lai_struc(n)%nr
     do c = 1,MODIS_lai_struc(n)%nc
        if (MODIS_lai_struc(n)%qcflag.eq.1) then !apply QC flag
           if (lai(c,r).gt.0.and.lai(c,r).le.100) then
              if ((MOD(flag(c,r),2).eq.0).and.(flag(c,r).le.62)) then
                 lai_flagged(c,r) = lai(c,r)*0.1
              else
                 lai_flagged(c,r) = LIS_rc%udef
               endif
            else
              lai_flagged(c,r) = LIS_rc%udef
           endif
        else  ! no QC flag applied
           if ((lai(c,r).gt.0).and.(lai(c,r).le.100)) then
              lai_flagged(c,r) = lai(c,r)*0.1
           else
              lai_flagged(c,r) = LIS_rc%udef
           endif
        endif
     enddo
  enddo

  do r = 1,MODIS_lai_struc(n)%nr
     do c = 1,MODIS_lai_struc(n)%nc
        lai_in(c+(r-1)*MODIS_lai_struc(n)%nc) = lai_flagged(c,r)
        if (lai_flagged(c,r).ne.LIS_rc%udef) then
           lai_data_b(c+(r-1)*MODIS_lai_struc(n)%nc) = .true.
        else
           lai_data_b(c+(r-1)*MODIS_lai_struc(n)%nc) = .false.
        endif
     enddo
  enddo

  if (LIS_isAtAfinerResolution(n,MODIS_lai_struc(n)%gridDesci(9))) then
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!--------------------------------------------------------------------------
     call bilinear_interp(LIS_rc%gridDesc(n,:),&
          lai_data_b, lai_in, laiobs_b_ip, laiobs_ip, &
          MODIS_lai_struc(n)%nc*MODIS_lai_struc(n)%nr, &
          LIS_rc%lnc(n)*LIS_rc%lnr(n), &
          MODIS_lai_struc(n)%rlat,MODIS_lai_struc(n)%rlon,&
          MODIS_lai_struc(n)%w11,MODIS_lai_struc(n)%w12,&
          MODIS_lai_struc(n)%w21,MODIS_lai_struc(n)%w22,&
          MODIS_lai_struc(n)%n11,MODIS_lai_struc(n)%n12,&
          MODIS_lai_struc(n)%n21,MODIS_lai_struc(n)%n22,LIS_rc%udef,ios)
  else
     call upscaleByAveraging(MODIS_lai_struc(n)%nc*MODIS_lai_struc(n)%nr,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n), &
          LIS_rc%udef, MODIS_lai_struc(n)%n11,&
          lai_data_b,lai_in,laiobs_b_ip,laiobs_ip)
  endif

  if (MODIS_lai_struc(n)%climofill.eq.1) then

     write(LIS_logunit,*) '[INFO] Opening climo file ',trim(climofile)
     ios = nf90_open(path=trim(climofile),mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error opening file '//trim(climofile))

     ios = nf90_inq_varid(nid,'Lai_500m',laiid)
     call LIS_verify(ios,'Error nf90_inq_varid: Lai_500m')

     cornerlat(1) = MODIS_lai_struc(n)%gridDesci(4)
     cornerlon(1) = MODIS_lai_struc(n)%gridDesci(5)
     cornerlat(2) = MODIS_lai_struc(n)%gridDesci(7)
     cornerlon(2) = MODIS_lai_struc(n)%gridDesci(8)

     lai_data_b = .false.

     lat_off = nint((cornerlat(1) +  89.9979167)/0.00416667)+1
     lon_off = nint((cornerlon(1) + 179.9979167)/0.00416667)+1

     ios = nf90_get_var(nid, laiid, lai, &
          start=(/lon_off,lat_off/), &
          count=(/MODIS_lai_struc(n)%nc,MODIS_lai_struc(n)%nr/))

     call LIS_verify(ios,'Error nf90_get_var: Lai_500m')
     ios = nf90_close(ncid=nid)
     call LIS_verify(ios,'Error closing file '//trim(fname))

     do r = 1,MODIS_lai_struc(n)%nr
        do c = 1,MODIS_lai_struc(n)%nc
           if ((lai(c,r).gt.0).and.(lai(c,r).le.100)) then
              lai_flagged(c,r) = lai(c,r)*0.1
           else
              lai_flagged(c,r) = LIS_rc%udef
           endif
        enddo
     enddo

     do r = 1,MODIS_lai_struc(n)%nr
        do c = 1,MODIS_lai_struc(n)%nc
           lai_in(c+(r-1)*MODIS_lai_struc(n)%nc) = lai_flagged(c,r)
           if (lai_flagged(c,r).ne.LIS_rc%udef) then
              lai_data_b(c+(r-1)*MODIS_lai_struc(n)%nc) = .true.
           else
              lai_data_b(c+(r-1)*MODIS_lai_struc(n)%nc) = .false.
           endif
        enddo
     enddo

     if (LIS_isAtAfinerResolution(n,MODIS_lai_struc(n)%gridDesci(9))) then
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!--------------------------------------------------------------------------
        call bilinear_interp(LIS_rc%gridDesc(n,:),&
             lai_data_b, lai_in, laiobs_b_ip, laiobs_climo_ip, &
             MODIS_lai_struc(n)%nc*MODIS_lai_struc(n)%nr, &
             LIS_rc%lnc(n)*LIS_rc%lnr(n), &
             MODIS_lai_struc(n)%rlat,MODIS_lai_struc(n)%rlon,&
             MODIS_lai_struc(n)%w11,MODIS_lai_struc(n)%w12,&
             MODIS_lai_struc(n)%w21,MODIS_lai_struc(n)%w22,&
             MODIS_lai_struc(n)%n11,MODIS_lai_struc(n)%n12,&
             MODIS_lai_struc(n)%n21,MODIS_lai_struc(n)%n22,LIS_rc%udef,ios)
     else
        call upscaleByAveraging(MODIS_lai_struc(n)%nc*MODIS_lai_struc(n)%nr,&
             LIS_rc%lnc(n)*LIS_rc%lnr(n), &
             LIS_rc%udef, MODIS_lai_struc(n)%n11,&
             lai_data_b,lai_in, laiobs_b_ip, laiobs_climo_ip)
     endif

     do r = 1,LIS_rc%lnr(n)
        do c = 1,LIS_rc%lnc(n)
           if ((laiobs_ip(c,r).eq.-9999.0).and. &
               (laiobs_climo_ip(c,r).ne.-9999.0)) then
              laiobs_ip(c,r) = laiobs_climo_ip(c,r)
           endif
        enddo
     enddo
  endif
#endif

end subroutine read_MODIS_LAI_data


!BOP
! !ROUTINE: create_MODIS_lai_filename
! \label{create_MODIS_lai_filename}
!
! !INTERFACE:
subroutine create_MODIS_lai_filename(ndir, yr, doy, filename, climofile)
! !USES:

  implicit none
! !ARGUMENTS:
  character(len=*)  :: filename
  character(len=*)  :: climofile
  integer           :: yr, doy
  character(len=*)  :: ndir
!
! !DESCRIPTION:
!  This subroutine creates the MCD15A2H LAI filename
!  based on the time and date
!
!  The arguments are:
!  \begin{description}
!  \item[ndir] name of the MCD15A2H LAI data directory
!  \item[yr]  current year
!  \item[doy]  current day of the year
!  \item[filename] Generated MCD15A2H LAI filename
!  \item[climofile] Generated MCD152AH LAI climatology file
!  \end{description}
!EOP

  character(len=4) :: fyr
  character(len=3) :: fdoy

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  filename = trim(ndir)//'/'//trim(fyr)//'/MCD15A2H.006_LAI_'//&
             trim(fyr)//trim(fdoy)//'.nc4'

  climofile = trim(ndir)//'/MCD15A2H.006_LAI_YYYY'//&
       trim(fdoy)//'.nc4'

end subroutine create_MODIS_lai_filename

