!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_nldas30
! \label{read_nldas30}
!
! !REVISION HISTORY:
! 16 May 2025: David Mocko, Initial Specification
!                           (derived from read_merra2.F90)
!
! !INTERFACE:
subroutine read_nldas30(n,order,month,findex,filename,nldasforc,ferror)
! !USES:
   use LIS_logMod
   use LIS_FORC_AttributesMod
!<debug -- jim testing>
   !use LIS_coreMod,       only  : LIS_rc,LIS_domain,LIS_masterproc
   use LIS_coreMod,       only  : LIS_rc,LIS_domain,LIS_masterproc,LIS_ews_ind, LIS_nss_ind,LIS_localPet
!</debug -- jim testing>
   use LIS_metforcingMod, only  : LIS_forc
   use nldas30_forcingMod, only : nldas30_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif

   implicit none
! !ARGUMENTS:
   integer, intent(in)          :: n
   integer, intent(in)          :: order
   integer, intent(in)          :: month
   integer, intent(in)          :: findex
   character(len=*), intent(in) :: filename
   real, intent(inout)          :: nldasforc(nldas30_struc(n)%nvars,24, &
      LIS_rc%lnc(n)*LIS_rc%lnr(n))
   integer, intent(out)         :: ferror

! !DESCRIPTION:
!  For the given time, reads parameters from NLDAS-3 data,
!  transforms into 8 LIS forcing parameters and interpolates
!  to the LIS domain. \newline
!
! NLDAS-3 FORCING VARIABLES (unless noted, fields are 1-hr): \newline
!  1. T 2m    Temperature interpolated to 2 metres [$K$] \newline
!  2. q 2m    Instantaneous specific humidity interpolated to 2 metres[$kg/kg$] \newline
!  3. swdn    Downward shortwave flux at the ground [$W/m^2$] \newline
!  4. lwdn    Downward longwave radiation at the ground [$W/m^2$] \newline
!  5. u 10m   Instantaneous zonal wind interpolated to 10 metres [$m/s$] \newline
!  6. v 10m   Instantaneous meridional wind interpolated to 10 metres[$m/s$] \newline
!  7. ps      Instantaneous Surface Pressure [$Pa$] \newline
!  8. preacc  Total precipitation [$mm/s$] \newline
!
!  The arguments are:
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous
!    1 hourly instance, order=2, read the next 1 hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the daily NLDAS-3 forcing file
!  \item[tscount]
!    time step count
!  \item[ferror]
!    return error code (0 indicates success)
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!EOP

   integer   :: ftn
   integer   :: tmpId,qId,uId,vId,psId,rainfId,swdnId,lwdnId
   integer   :: nr_index,nc_index
   logical   :: file_exists
   integer   :: c,r,t,iret
   integer   :: mo
   logical   :: read_lnd

   real      :: tair(LIS_rc%lnc(n),LIS_rc%lnr(n),24)
   real      :: qair(LIS_rc%lnc(n),LIS_rc%lnr(n),24)
   real      :: uwind(LIS_rc%lnc(n),LIS_rc%lnr(n),24)
   real      :: vwind(LIS_rc%lnc(n),LIS_rc%lnr(n),24)
   real      :: ps(LIS_rc%lnc(n),LIS_rc%lnr(n),24)
   real      :: rainf(LIS_rc%lnc(n),LIS_rc%lnr(n),24)
   real      :: swdn(LIS_rc%lnc(n),LIS_rc%lnr(n),24)
   real      :: lwdn(LIS_rc%lnc(n),LIS_rc%lnr(n),24)

!<debug -- jim testing>
   integer :: sC, sR, cC, cR
!</debug -- jim testing>
!_______________________________________________________________________

#if (defined USE_NETCDF3)
   write(LIS_logunit,*) "[ERR] NLDAS-3 reader requires NetCDF4"
   call LIS_endrun
#endif

#if (defined USE_NETCDF4)
   ferror = 0
   mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

   ! Read NLDAS-3 fields
   inquire(file=filename,exist=file_exists)
   if (file_exists) then
      write(LIS_logunit,*) '[INFO] Reading NLDAS-3 file (bookend,',     &
         order,' ... ',trim(filename)
      call LIS_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,  &
         ncid=ftn),'nf90_open failed in read_nldas30')

!<debug -- jim testing>
      ! I need to compute the size of each subdomain.
      ! For now, that will be lnc,lnr.
      ! I need indices into global grid for the start array.
      sC = LIS_ews_ind(n,LIS_localPet+1) + 1 - 1
      sR = LIS_nss_ind(n,LIS_localPet+1) + 1 - 1
      cC = LIS_rc%lnc(n)
      cR = LIS_rc%lnr(n)
!</debug -- jim testing>

      call LIS_verify(nf90_inq_varid(ftn,'Tair',tmpId), &
         'nf90_inq_varid failed for Tair in read_nldas30')
      call LIS_verify(nf90_inq_varid(ftn,'Qair',qId), &
         'nf90_inq_varid failed for Qair in read_nldas30')
      call LIS_verify(nf90_inq_varid(ftn,'PSurf',psId), &
         'nf90_inq_varid failed for PSurf in read_nldas30')
      call LIS_verify(nf90_inq_varid(ftn,'LWdown',lwdnId), &
         'nf90_inq_varid failed for LWdown in read_nldas30')
      call LIS_verify(nf90_inq_varid(ftn,'SWdown',swdnId), &
         'nf90_inq_varid failed for SWdown in read_nldas30')
      call LIS_verify(nf90_inq_varid(ftn,'Wind_E',uId), &
         'nf90_inq_varid failed for Wind_E in read_nldas30')
      call LIS_verify(nf90_inq_varid(ftn,'Wind_N',vId), &
         'nf90_inq_varid failed for Wind_N in read_nldas30')
      call LIS_verify(nf90_inq_varid(ftn,'Rainf',rainfId), &
         'nf90_inq_varid failed for Rainf in read_nldas30')

      call LIS_verify(nf90_get_var(ftn, tmpId, tair, &
         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
         'nf90_get_var failed for tair in read_nldas30')
      call LIS_verify(nf90_get_var(ftn, qId, qair, &
         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
         'nf90_get_var failed for qair in read_nldas30')
      call LIS_verify(nf90_get_var(ftn, psId, ps, &
         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
         'nf90_get_var failed for ps in read_nldas30')
      call LIS_verify(nf90_get_var(ftn, lwdnId, lwdn, &
         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
         'nf90_get_var failed for lwdn in read_nldas30')
      call LIS_verify(nf90_get_var(ftn, swdnId, swdn, &
         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
         'nf90_get_var failed for swdn in read_nldas30')
      call LIS_verify(nf90_get_var(ftn, uId, uwind, &
         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
         'nf90_get_var failed for uwind in read_nldas30')
      call LIS_verify(nf90_get_var(ftn, vId, vwind, &
         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
         'nf90_get_var failed for vwind in read_nldas30')
      call LIS_verify(nf90_get_var(ftn, rainfId, rainf, &
         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
         'nf90_get_var failed for rainf in read_nldas30')

      call LIS_verify(nf90_close(ftn),'failed to close in read_nldas30')

      call interp_nldas30_var(n,findex,month,tair,  1,.false.,nldasforc)
      call interp_nldas30_var(n,findex,month,qair,  2,.false.,nldasforc)
      call interp_nldas30_var(n,findex,month,ps,    7,.false.,nldasforc)
      call interp_nldas30_var(n,findex,month,lwdn,  4,.false.,nldasforc)
      call interp_nldas30_var(n,findex,month,swdn,  3,.false.,nldasforc)
      call interp_nldas30_var(n,findex,month,uwind, 5,.false.,nldasforc)
      call interp_nldas30_var(n,findex,month,vwind, 6,.false.,nldasforc)
      call interp_nldas30_var(n,findex,month,rainf, 8, .true.,nldasforc)
   else
      write(LIS_logunit,*) '[ERR] ',trim(filename)//' does not exist'
      call LIS_endrun
   endif
#endif
end subroutine read_nldas30

!BOP
!
! !ROUTINE: interp_nldas30_var
! \label{interp_nldas30_var}
!
! !INTERFACE:
#if 1
subroutine interp_nldas30_var(n,findex,month,input_var,var_index, &
   pcp_flag,nldasforc)

! !USES:
   use LIS_coreMod
   use LIS_logMod
   use LIS_spatialDownscalingMod
   use nldas30_forcingMod, only : nldas30_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif
   implicit none

! !ARGUMENTS:
   integer, intent(in)    :: n
   integer, intent(in)    :: findex
   integer, intent(in)    :: month
   real,    intent(in)    :: input_var(LIS_rc%lnc(n), LIS_rc%lnr(n),24)
   integer, intent(in)    :: var_index
   logical, intent(in)    :: pcp_flag
   real,    intent(inout) :: nldasforc(nldas30_struc(n)%nvars,24, &
      LIS_rc%lnc(n)*LIS_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine spatially interpolates a NLDAS-3 field
!  to the LIS running domain
!
!EOP
   integer   :: t,c,r,k,iret
   integer   :: doy
   integer   :: ftn
   integer   :: pcp1Id,pcp2Id,pcp3Id,pcp4Id,pcp5Id,pcp6Id
   real      :: f (LIS_rc%lnc(n)*LIS_rc%lnr(n))
   logical*1 :: lb(LIS_rc%lnc(n)*LIS_rc%lnr(n))
   logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
   integer   :: input_size
   logical   :: scal_read_flag
! _____________________________________________________________

   input_size = LIS_rc%lnc(n)*LIS_rc%lnr(n)

   do t = 1,24
      lb = .true.
      do r = 1, LIS_rc%lnr(n) !nldas30_struc(n)%nrold
         do c = 1, LIS_rc%lnc(n) !nldas30_struc(n)%ncold
            k = c+(r-1)*LIS_rc%lnc(n) !nldas30_struc(n)%ncold
            f(k) = input_var(c,r,t)
            if (f(k).eq.-9999.0) then
               f(k)  = LIS_rc%udef
               lb(k) = .false.
            endif
         enddo
      enddo

      if (pcp_flag.and. &
         trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
         call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,            &
            nldasforc(var_index,t,:),                                 &
            input_size,LIS_rc%lnc(n)*LIS_rc%lnr(n),                   &
            LIS_domain(n)%lat,LIS_domain(n)%lon,                      &
            nldas30_struc(n)%w112,nldas30_struc(n)%w122,              &
            nldas30_struc(n)%w212,nldas30_struc(n)%w222,              &
            nldas30_struc(n)%n112,nldas30_struc(n)%n122,              &
            nldas30_struc(n)%n212,nldas30_struc(n)%n222,              &
            LIS_rc%udef,iret)
      elseif ((trim(LIS_rc%met_interp(findex)).eq."bilinear").or. &
              (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear")) then
            call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,           &
               nldasforc(var_index,t,:),                                 &
               input_size,LIS_rc%lnc(n)*LIS_rc%lnr(n),                   &
               LIS_domain(n)%lat,LIS_domain(n)%lon,                      &
               nldas30_struc(n)%w111,nldas30_struc(n)%w121,              &
               nldas30_struc(n)%w211,nldas30_struc(n)%w221,              &
               nldas30_struc(n)%n111,nldas30_struc(n)%n121,              &
               nldas30_struc(n)%n211,nldas30_struc(n)%n221,              &
               LIS_rc%udef,iret)
      elseif (trim(LIS_rc%met_interp(findex)).eq."neighbor") then
         call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,           &
            nldasforc(var_index,t,:),input_size,                      &
            LIS_rc%lnc(n)*LIS_rc%lnr(n),                              &
            LIS_domain(n)%lat,LIS_domain(n)%lon,                      &
            nldas30_struc(n)%n113,LIS_rc%udef,iret)
      else
         write(LIS_logunit,*) '[ERR] Spatial interpolation option '// &
            trim(LIS_rc%met_interp(findex))// &
            ' not supported for NLDAS-3'
         call LIS_endrun
      endif
   enddo

end subroutine interp_nldas30_var
#endif
#if 0
subroutine interp_nldas30_var(n,findex,month,input_var,var_index, &
   pcp_flag,nldasforc)

! !USES:
   use LIS_coreMod
   use LIS_logMod
   use LIS_spatialDownscalingMod
   use nldas30_forcingMod, only : nldas30_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif
   implicit none

! !ARGUMENTS:
   integer, intent(in)    :: n
   integer, intent(in)    :: findex
   integer, intent(in)    :: month
   real,    intent(in)    :: input_var(LIS_rc%lnc(n), LIS_rc%lnr(n),24)
   integer, intent(in)    :: var_index
   logical, intent(in)    :: pcp_flag
   real,    intent(inout) :: nldasforc(nldas30_struc(n)%nvars,24, &
      LIS_rc%lnc(n)*LIS_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine spatially interpolates a NLDAS-3 field
!  to the LIS running domain
!
!EOP
   integer   :: t,c,r,k,iret
   integer   :: doy
   integer   :: ftn
   integer   :: pcp1Id,pcp2Id,pcp3Id,pcp4Id,pcp5Id,pcp6Id
   real      :: f (LIS_rc%lnc(n)*LIS_rc%lnr(n))
   logical*1 :: lb(LIS_rc%lnc(n)*LIS_rc%lnr(n))
   logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
   integer   :: input_size
   logical   :: scal_read_flag
! _____________________________________________________________

!<debug -- jim testing>
logical :: first = .true.
!</debug -- jim testing>
   input_size = LIS_rc%lnc(n)*LIS_rc%lnr(n)

   do t = 1,24
      lb = .true.
      do r = 1, LIS_rc%lnr(n) !nldas30_struc(n)%nrold
         do c = 1, LIS_rc%lnc(n) !nldas30_struc(n)%ncold
            k = c+(r-1)*LIS_rc%lnc(n) !nldas30_struc(n)%ncold
            f(k) = input_var(c,r,t)
            if (f(k).eq.1.e+15) then
!<debug -- jim testing>
write(unit=LIS_logunit,fmt=*) 'GREP: Should not be here.'
!</debug -- jim testing>
               f(k)  = LIS_rc%udef
               lb(k) = .false.
            endif
            if (f(k).eq.-9999.0) then
!<debug -- jim testing>
if ( first ) then
write(unit=LIS_logunit,fmt=*) 'GREP: This is expected.'
first = .false.
endif
!</debug -- jim testing>
               f(k)  = LIS_rc%udef
               lb(k) = .false.
            endif
         enddo
      enddo
      nldasforc(var_index,t,:) = f(:)
   enddo

end subroutine interp_nldas30_var
#endif
