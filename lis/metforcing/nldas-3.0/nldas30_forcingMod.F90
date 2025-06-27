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
   use bounding_box_mod

   implicit none

   PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
   public :: init_nldas30 !defines the native resolution of the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
   public :: nldas30_struc

!EOP
   type, public ::  nldas30_type_dec
      real         :: ts
      integer      :: ncold,nrold
      character(len=LIS_CONST_PATH_LEN) :: nldas30dir ! NLDAS-3 Forcing Directory
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

      type(bounding_box_type) :: bb
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
! 04 Jun 2025: James Geiger, add support to read subsets of NLDAS-3 domain
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
!   \item[neighbor\_interp\_input](\ref{neighbor_interp_input}) \newline
!    computes the neighbor, weights for nearest-neighbor interpolation
!   \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!    computes the neighbors for upscaling by averaging
!  \end{description}
!EOP
   real :: gridDesci(50)
   integer :: n
   real, allocatable, dimension(:) :: nldas30_lon, nldas30_lat
   type(bounding_box_type) :: bb

   LIS_rc%met_nf(findex) = 8

   allocate(nldas30_struc(LIS_rc%nnest))

   call readcrd_nldas30

   nldas30_struc(:)%ncold = 11700
   nldas30_struc(:)%nrold = 6500
   ! ncold and nrold are independent of nest
   nldas30_struc(:)%mi = nldas30_struc(1)%ncold*nldas30_struc(1)%nrold

   nldas30_struc(:)%reset_flag = .false.
   nldas30_struc(:)%startFlag = .true.
   nldas30_struc(:)%dayFlag = .true.
   nldas30_struc(:)%nvars = 8
   nldas30_struc(:)%st_iterid = 1
   nldas30_struc(:)%en_iterId = 1
   nldas30_struc(:)%nIter = 1

   nldas30_struc(:)%ts = 3600
   do n = 1,LIS_rc%nnest
      call LIS_update_timestep(LIS_rc,n,nldas30_struc(n)%ts)
   enddo

   ! NLDAS-3 domain grid description
   ! ncold and nrold are independent of nest
   gridDesci = 0
   gridDesci(1)  = 0
   gridDesci(2)  = nldas30_struc(1)%ncold
   gridDesci(3)  = nldas30_struc(1)%nrold
   gridDesci(4)  =    7.005
   gridDesci(5)  = -168.995
   gridDesci(6)  = 128
   gridDesci(7)  =   71.995
   gridDesci(8)  =  -52.005
   gridDesci(9)  =    0.01
   gridDesci(10) =    0.01
   gridDesci(20) = 0

   ! ncold and nrold are independent of nest
   allocate(nldas30_lon(nldas30_struc(1)%ncold))
   allocate(nldas30_lat(nldas30_struc(1)%nrold))

   call nldas30_earth_coords(gridDesci, nldas30_lon, nldas30_lat)

   do n = 1,LIS_rc%nnest
      ! ncold and nrold are independent of nest
      bb = find_bounding_box(nldas30_struc(1)%ncold, nldas30_struc(1)%nrold, nldas30_lat, nldas30_lon, minval(LIS_domain(n)%lat), maxval(LIS_domain(n)%lat), minval(LIS_domain(n)%lon), maxval(LIS_domain(n)%lon))

      gridDesci(2)  = bb%NLON
      gridDesci(3)  = bb%NLAT
      gridDesci(4)  = nldas30_lat(bb%i_llat)
      gridDesci(5)  = nldas30_lon(bb%i_llon)
      gridDesci(7)  = nldas30_lat(bb%i_ulat)
      gridDesci(8)  = nldas30_lon(bb%i_ulon)

      nldas30_struc(n)%bb = bb

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
         call bilinear_interp_input(n,gridDesci,            &
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
         call bilinear_interp_input(n,gridDesci,            &
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
         call conserv_interp_input(n,gridDesci,              &
            nldas30_struc(n)%n112,nldas30_struc(n)%n122,     &
            nldas30_struc(n)%n212,nldas30_struc(n)%n222,     &
            nldas30_struc(n)%w112,nldas30_struc(n)%w122,     &
            nldas30_struc(n)%w212,nldas30_struc(n)%w222)

      elseif (trim(LIS_rc%met_interp(findex)).eq."neighbor") then
         allocate(nldas30_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         call neighbor_interp_input(n,gridDesci,nldas30_struc(n)%n113)

      elseif (trim(LIS_rc%met_interp(findex)).eq."average") then
         allocate(nldas30_struc(n)%n111(bb%NLON*bb%NLAT))
         call upscaleByAveraging_input(gridDesci,&
            LIS_rc%gridDesc(n,:), bb%NLON*bb%NLAT,&
            LIS_rc%lnc(n)*LIS_rc%lnr(n), nldas30_struc(n)%n111)

      else
         write(LIS_logunit,*) '[ERR] Interpolation option '//           &
            trim(LIS_rc%met_interp(findex))//        &
            ' for NLDAS-3 forcing is not supported'
         call LIS_endrun
      endif

      allocate(nldas30_struc(n)%nldasforc1(1,nldas30_struc(n)%nvars,24, &
         LIS_rc%lnc(n)*LIS_rc%lnr(n)))
      allocate(nldas30_struc(n)%nldasforc2(1,nldas30_struc(n)%nvars,24, &
         LIS_rc%lnc(n)*LIS_rc%lnr(n)))

      allocate(nldas30_struc(n)%metdata1(1,LIS_rc%met_nf(findex),       &
         LIS_rc%ngrid(n)))
      allocate(nldas30_struc(n)%metdata2(1,LIS_rc%met_nf(findex),       &
         LIS_rc%ngrid(n)))

      nldas30_struc(n)%metdata1 = 0
      nldas30_struc(n)%metdata2 = 0

      nldas30_struc(n)%nldasforc1 = LIS_rc%udef
      nldas30_struc(n)%nldasforc2 = LIS_rc%udef
   enddo

   deallocate(nldas30_lon)
   deallocate(nldas30_lat)
end subroutine init_nldas30

!BOP
!
! !ROUTINE: nldas30_earth_coords
! \label{nldas30_earth_coords}
!
! !REVISION HISTORY:
! 04 Jun 2025: James Geiger, initial specification
!
! !INTERFACE:
subroutine nldas30_earth_coords(gridDesci, nldas3_lon, nldas3_lat)
! !USES:
!  none
!
   implicit none
!
! !AGRUMENTS:
   real :: gridDesci(50)
   real, dimension(gridDesci(2)) :: nldas3_lon
   real, dimension(gridDesci(3)) :: nldas3_lat
!
! !DESCRIPTION:
!  Compute an array of latitudes and an array of longitudes
!  for the full NLDAS-3 domain.
!
!  \begin{description}
!  \item[gridDesci]
!    grid description array for the NLDAS-3 domain
!  \item[nldas3_lon]
!    Array of longitude values for the NLDAS-3 domain
!  \item[nldas3_lat]
!    Array of latitude values for the NLDAS-3 domain
!  \end{description}
!
!EOP

   integer :: i, NC, NR
   real :: res, start

   NC = gridDesci(2)
   NR = gridDesci(3)

   res = gridDesci(9)
   start = gridDesci(5)
   do i = 1, NC
      nldas3_lon(i) = start + (i-1)*res
   enddo

   res = gridDesci(10)
   start = gridDesci(4)
   do i = 1, NR
      nldas3_lat(i) = start + (i-1)*res
   enddo
end subroutine nldas30_earth_coords

end module nldas30_forcingMod
