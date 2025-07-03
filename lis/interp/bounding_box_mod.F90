!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module bounding_box_mod
!BOP
!
! !MODULE: bounding_box_mod
!
! !DESCRIPTION:
! This module contains the data structures and routines for computing
! a bounding box around given latitude and longitude values.
!
! Example:
! Say that you have a large metforcing dataset.
! The grid for this metforcing dataset is the input domain.
! The running grid for a given LIS process is the output domain.
! This module will find the smallest sub-domain of the input domain
! that contains the output domain.  I.e., it will find the smallest
! sub-domain of the metforcing dataset that contains the LIS running
! domain for that process.  This allows the metforcing reader to read
! a subset of the input metforcing data to be spatially interpolated
! to the LIS running domain for that process.
!
! Caveats:
! This has been tested only on equidistant cylindrical projections
! for the input and output grids.
!
! !REVISION HISTORY:
! 05 Jun 2024: James Geiger, Initial Specification
!                                 
! !USES:
! none
!
! EOP
type bounding_box_type
   integer :: i_llat ! lower latitude index
   integer :: i_ulat ! upper latitude index
   integer :: i_llon ! lower longitude index
   integer :: i_ulon ! upper longitude index
   integer :: NLAT ! number of latitude points
   integer :: NLON ! number of longitude points
end type

contains

!BOP
!
! !ROUTINE: find_bounding_box
! \label{find_bounding_box
!
! !REVISION HISTORY:
! 05 Jun 2024: James Geiger, Initial Specification
!
! !INTERFACE:
function find_bounding_box(NC, NR, ilat, ilon, min_out_lat, max_out_lat, min_out_lon, max_out_lon) result(bb)
! !USES:
   use LIS_logMod, only : LIS_logunit, LIS_endrun

   implicit none

! !ARGUMENTS:
   integer, intent(in) :: NC, NR
   real, dimension(NR), intent(in) :: ilat
   real, dimension(NC), intent(in) :: ilon
   real, intent(in) :: min_out_lat, max_out_lat, min_out_lon, max_out_lon

! !DESCRIPTION:
! This function finds the bounding box within the input domain
! that contains the given output domain.  This function returns
! a bounding box datatype.
!
! The arguments are:
! \begin{description}
! \item[NC]
!    number of columns in the input domain
! \item[NR]
!    number of rows in the input domain
! \item[ilat]
!    array of latitude values for the input domain
! \item[ilon]
!    array of longitude values for the input domain
! \item[min\_out\_lat]
!    minimum latitude value for the output domain
! \item[max\_out\_lat]
!    maximum latitude value for the output domain
! \item[min\_out\_lon]
!    minimum longitude value for the output domain
! \item[max\_out\_lon]
!    maximum longitude value for the output domain
! \end{description}
!
!EOP

   type(bounding_box_type) :: bb
   real :: min_olat, max_olat, min_olon, max_olon
   real, parameter :: epsilon = 0.0001

   min_olat = min_out_lat + epsilon
   max_olat = max_out_lat - epsilon
   min_olon = min_out_lon + epsilon
   max_olon = max_out_lon - epsilon

   if ( min_olat < ilat(1) .or. &
        max_olat > ilat(NR) .or. &
        min_olon < ilon(1) .or. &
        max_olon > ilon(NC) ) then
      write(LIS_logunit,*) '[ERR] The output domain is outside the input domain.  (with epsilon ', epsilon, ')'
      write(LIS_logunit,'(a,2f14.8)') '[ERR] min_out_lat, ilat(1) ', min_out_lat, ilat(1)
      write(LIS_logunit,'(a,2f14.8)') '[ERR] max_out_lat, ilat(NR) ' , max_out_lat, ilat(NR)
      write(LIS_logunit,'(a,2f14.8)') '[ERR] min_out_lon, ilon(1) ', min_out_lon, ilon(1)
      write(LIS_logunit,'(a,2f14.8)') '[ERR] max_out_lon, ilon(NC) ', max_out_lon, ilon(NC)
      call LIS_endrun
   endif

   bb%i_llat = find_lower(NR, ilat, min_olat)
   bb%i_ulat = find_upper(NR, ilat, max_olat)
   bb%i_llon = find_lower(NC, ilon, min_olon)
   bb%i_ulon = find_upper(NC, ilon, max_olon)
   bb%NLAT = bb%i_ulat - bb%i_llat + 1
   bb%NLON = bb%i_ulon - bb%i_llon + 1

   contains

   function find_lower(Ni, input, output) result(i)
      implicit none

      integer :: Ni
      real, dimension(Ni) :: input
      real :: output

      integer :: i

      i = 1
      do while ( input(i) <= output )
         i = i + 1
      enddo 
      i = i - 1
   end function find_lower

   function find_upper(Ni, input, output) result(i)
      implicit none

      integer :: Ni
      real, dimension(Ni) :: input
      real :: output

      integer :: i

      i = Ni
      do while ( input(i) >= output )
         i = i - 1
      enddo 
      i = i + 1
   end function find_upper
end function find_bounding_box
end module bounding_box_mod
