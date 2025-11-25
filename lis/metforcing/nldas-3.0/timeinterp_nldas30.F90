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
! !ROUTINE: timeinterp_nldas30
! \label{timeinterp_nldas30}
!
! !REVISION HISTORY:
! 16 May 2025: David Mocko, Initial Specification
!                           (derived from timeinterp_merra2.F90)
!
! !INTERFACE:
subroutine timeinterp_nldas30(n,findex)

! !USES:
   use ESMF
   use LIS_FORC_AttributesMod
   use LIS_coreMod,          only : LIS_rc,LIS_domain,LIS_localPet
   use LIS_metforcingMod,    only : LIS_FORC_Base_State
   use LIS_constantsMod,     only : LIS_CONST_SOLAR
   use LIS_logMod,           only : LIS_logunit,LIS_verify,LIS_endrun
   use LIS_forecastMod,      only : LIS_get_iteration_index
   use nldas30_forcingMod,   only : nldas30_struc

   implicit none

! !ARGUMENTS:
   integer, intent(in):: n
   integer, intent(in):: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model
!  timestep. Downward shortwave radiation is interpolated using a
!  zenith-angled based approach. Precipitation and longwave radiation
!  are not temporally interpolated, and the previous 1 hourly value
!  is used. All other variables are linearly interpolated between
!  the 1 hourly blocks.
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!  \end{description}
!EOP
   integer :: t,k,kk
   integer :: index1
   real    :: wt1,wt2
   integer            :: status
   integer            :: mfactor,m
   type(ESMF_Field)   :: tmpField,q2Field,uField,vField
   type(ESMF_Field)   :: psurfField,pcpField,swdField,lwdField
   real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
   real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:)

   !-----------------------------------------------------------------------
   !  Interpolate Data in Time
   !-----------------------------------------------------------------------

   wt1 = (nldas30_struc(n)%nldas30time2-LIS_rc%time)/                   &
      (nldas30_struc(n)%nldas30time2-nldas30_struc(n)%nldas30time1)
   wt2 = 1.0 - wt1
   if (wt1.gt.1.01) then
      wt1 = 1.0
      wt2 = 0.0
   endif

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),tmpField,&
      rc=status)
   call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(1),q2Field,&
      rc=status)
   call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1),swdField,&
      rc=status)
   call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LWdown%varname(1),lwdField,&
      rc=status)
   call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uField,&
      rc=status)
   call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vField,&
      rc=status)
   call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(1),psurfField,&
      rc=status)
   call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),pcpField,&
      rc=status)
   call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')

   call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
   call LIS_verify(status)

   call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
   call LIS_verify(status)

   call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
   call LIS_verify(status)

   call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
   call LIS_verify(status)

   call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
   call LIS_verify(status)

   call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
   call LIS_verify(status)

   call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
   call LIS_verify(status)

   call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
   call LIS_verify(status)

   mfactor = LIS_rc%nensem(n) / nldas30_struc(n)%nIter

   do k = 1,LIS_rc%ntiles(n)/mfactor
      do m = 1,mfactor
         t = m + (k-1)*mfactor
         index1 = LIS_domain(n)%tile(t)%index
         kk = LIS_get_iteration_index(n,k,index1,mfactor)
         if ((nldas30_struc(n)%metdata1(kk,3,index1).ne.LIS_rc%udef).and.&
             (nldas30_struc(n)%metdata2(kk,3,index1).ne.LIS_rc%udef)) then
             if (nldas30_struc(n)%metdata1(kk,3,index1).lt.0.0) then
!                write(unit=LIS_logunit,fmt=*) '[WARN] SWdown negative'
!                write(unit=LIS_logunit,fmt=*) '[WARN] metdata1 = ',     &
!                      nldas30_struc(n)%metdata1(kk,3,index1)
!                write(unit=LIS_logunit,fmt=*) '[WARN] Resetting to zero'
                nldas30_struc(n)%metdata1(kk,3,index1) = 0.0
             endif
             if (nldas30_struc(n)%metdata2(kk,3,index1).lt.0.0) then
!                write(unit=LIS_logunit,fmt=*) '[WARN] SWdown negative'
!                write(unit=LIS_logunit,fmt=*) '[WARN] metdata2 = ',     &
!                      nldas30_struc(n)%metdata2(kk,3,index1)
!                write(unit=LIS_logunit,fmt=*) '[WARN] Resetting to zero'
                nldas30_struc(n)%metdata2(kk,3,index1) = 0.0
             endif
             swd(t) = (nldas30_struc(n)%metdata1(kk,3,index1) * wt1) +  &
                      (nldas30_struc(n)%metdata2(kk,3,index1) * wt2)
            if (swd(t).gt.LIS_CONST_SOLAR) then
               write(unit=LIS_logunit,fmt=*)                            &
                  '[WARN] sw radiation too high!!'
               write(unit=LIS_logunit,fmt=*) '[WARN] it is',swd(t)
               write(unit=LIS_logunit,fmt=*) '[WARN] nldas30data1 = ',  &
                  nldas30_struc(n)%metdata1(kk,3,index1)
               write(unit=LIS_logunit,fmt=*) '[WARN] nldas30data2 = ',  &
                  nldas30_struc(n)%metdata2(kk,3,index1)
               write(unit=LIS_logunit,fmt=*) '[WARN] wt1=',wt1,'wt2=',wt2
               swd(t) = LIS_CONST_SOLAR
               write(unit=LIS_logunit,fmt=*) '[WARN] forcing set to ',swd(t)
            endif
         endif

         if ((swd(t).ne.LIS_rc%udef).and.(swd(t).lt.0.0)) then
            if (swd(t).gt.-0.00001) then
               swd(t) = 0.0
            else
               write(LIS_logunit,*)                                     &
                  '[ERR] timeinterp_nldas30 -- Stopping because ',    &
                  'forcing not udef but lt0,'
               write(LIS_logunit,*)'[ERR] timeinterp_nldas30 -- ',      &
                  t,swd(t),' (',LIS_localPet,')'
               write(LIS_logunit,*) &
                  nldas30_struc(n)%metdata2(kk,3,index1),    &
                  nldas30_struc(n)%metdata1(kk,3,index1),    &
                  wt1,wt2
               call LIS_endrun
            endif
         endif

         if (swd(t).gt.LIS_CONST_SOLAR) then
            swd(t) = nldas30_struc(n)%metdata2(kk,3,index1)
         endif
      enddo
   enddo

   !-----------------------------------------------------------------------
   ! precip variable - constant rate over the NLDAS-3 hour
   !-----------------------------------------------------------------------
   do k = 1,LIS_rc%ntiles(n)/mfactor
      do m = 1,mfactor
         t = m + (k-1)*mfactor
         index1 = LIS_domain(n)%tile(t)%index
         kk = LIS_get_iteration_index(n,k,index1,mfactor)
         ! Total precip is a constant amount.
         ! Based on the first/current time and not on the latter time.
         pcp(t) = nldas30_struc(n)%metdata1(kk,8,index1)
         pcp(t) = pcp(t) / (60.0*60.0)

         if ( pcp(t) .lt. 0 ) then
            pcp(t) = 0
         endif
      enddo
   enddo

   !-----------------------------------------------------------------------
   ! LW down - linearly interpolate
   !-----------------------------------------------------------------------
   do k = 1,LIS_rc%ntiles(n)/mfactor
      do m = 1,mfactor
         t = m + (k-1)*mfactor
         index1 = LIS_domain(n)%tile(t)%index
         kk = LIS_get_iteration_index(n,k,index1,mfactor)
         lwd(t) = nldas30_struc(n)%metdata1(kk,4,index1)
         lwd(t) = (nldas30_struc(n)%metdata1(kk,4,index1) * wt1) + &
                  (nldas30_struc(n)%metdata2(kk,4,index1) * wt2)
      enddo
   enddo

   !-----------------------------------------------------------------------
   ! Bookend 1 for everything else
   !-----------------------------------------------------------------------
   do k = 1,LIS_rc%ntiles(n)/mfactor
      do m = 1,mfactor
         t = m + (k-1)*mfactor
         index1 = LIS_domain(n)%tile(t)%index
         kk = LIS_get_iteration_index(n,k,index1,mfactor)
         tmp(t) = nldas30_struc(n)%metdata1(kk,1,index1)
         q2(t)  = nldas30_struc(n)%metdata1(kk,2,index1)
         uwind(t) = nldas30_struc(n)%metdata1(kk,5,index1)
         vwind(t) = nldas30_struc(n)%metdata1(kk,6,index1)
         psurf(t) = nldas30_struc(n)%metdata1(kk,7,index1)
      enddo
   enddo

end subroutine timeinterp_nldas30
