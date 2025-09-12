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
! for the input grid.
!
! !REVISION HISTORY:
! 05 Jun 2025: James Geiger, Initial Specification
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
! !ROUTINE: remove_zeros
! \label{remove_zeros}
!
! !REVISION HISTORY:
! 12 Sep 2025: James Geiger, Initial Specification
!
! !INTERFACE:
function remove_zeros(N, x) result(r)
! !USES:
!  none

   implicit none

! !ARGUMENTS:
   integer, intent(in) :: N
   integer, dimension(N), intent(in) :: x
   integer, pointer, dimension(:) :: r

! !DESCRIPTION:
! This function takes a 1-D array of integers and creates an output
! array by removing the elements that are 0.
!
! The arguments are:
! \begin{description}
! \item[N]
!    size of the input array
! \item[x]
!    input array
! \item[r]
!    output array
! \end{description}
!
!EOP

   integer :: i, c

   allocate(r(count(x /= 0)))
   c = 1
   do i = 1, N
      if ( x(i) /= 0 ) then
         r(c) = x(i)
         c = c + 1
      endif
   enddo
end function remove_zeros

!BOP
!
! !ROUTINE: find_bounding_box_bilinear
! \label{find_bounding_box_bilinear}
!
! !REVISION HISTORY:
! 17 Jul 2025: James Geiger, Initial Specification
!
! !INTERFACE:
function find_bounding_box_bilinear(INC, INR, ilat, ilon, ONCxONR, n111, n121, n211, n221) result(bb)
! !USES:
! none

   implicit none

! !ARGUMENTS:
   integer, intent(in) :: INC, INR, ONCxONR
   real, dimension(INR), intent(in) :: ilat
   real, dimension(INC), intent(in) :: ilon
   integer, dimension(ONCxONR), intent(in) :: n111, n121, n211, n221

! !DESCRIPTION:
! This function finds the bounding box within the input domain
! that contains the given output domain.  This function returns
! a bounding box datatype.
!
! The arguments are:
! \begin{description}
! \item[INC]
!    number of columns in the input domain
! \item[INR]
!    number of rows in the input domain
! \item[ilat]
!    array of latitude values for the input domain
! \item[ilon]
!    array of longitude values for the input domain
! \item[ONCxONR]
!    maximum number of grid-cells in the output domain
!    (number of columns times number of rows)
! \item[n111]
!    n111 array from bilinear_interp_input
! \item[n121]
!    n121 array from bilinear_interp_input
! \item[n211]
!    n211 array from bilinear_interp_input
! \item[n221]
!    n221 array from bilinear_interp_input
! \end{description}
!
!EOP

   integer, parameter :: BUFSIZE = 0
   type(bounding_box_type) :: bb

   integer :: min_c, max_c, min_r, max_r
   integer :: min_r111, max_r111, min_c111, max_c111
   integer :: min_r121, max_r121, min_c121, max_c121
   integer :: min_r211, max_r211, min_c211, max_c211
   integer :: min_r221, max_r221, min_c221, max_c221

   integer, allocatable, target, dimension(:) :: tmp_array

   tmp_array = remove_zeros(ONCxONR, n111)
   tmp_array = tmp_array - 1
   min_c111 = minval(mod(tmp_array, INC) + 1)
   max_c111 = maxval(mod(tmp_array, INC) + 1)
   min_r111 = minval(tmp_array / INC + 1)
   max_r111 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   tmp_array = remove_zeros(ONCxONR, n121)
   tmp_array = tmp_array - 1
   min_c121 = minval(mod(tmp_array, INC)+1)
   max_c121 = maxval(mod(tmp_array, INC)+1)
   min_r121 = minval(tmp_array / INC + 1)
   max_r121 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   tmp_array = remove_zeros(ONCxONR, n211)
   tmp_array = tmp_array - 1
   min_c211 = minval(mod(tmp_array, INC)+1)
   max_c211 = maxval(mod(tmp_array, INC)+1)
   min_r211 = minval(tmp_array / INC + 1)
   max_r211 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   tmp_array = remove_zeros(ONCxONR, n221)
   tmp_array = tmp_array - 1
   min_c221 = minval(mod(tmp_array, INC)+1)
   max_c221 = maxval(mod(tmp_array, INC)+1)
   min_r221 = minval(tmp_array / INC + 1)
   max_r221 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   min_r = min(min_r111, min_r121, min_r211, min_r221)
   max_r = max(max_r111, max_r121, max_r211, max_r221)
   min_c = min(min_c111, min_c121, min_c211, min_c221)
   max_c = max(max_c111, max_c121, max_c211, max_c221)

   bb%i_llat = min_r
   bb%i_ulat = max_r
   bb%i_llon = min_c
   bb%i_ulon = max_c
   bb%i_llat = max(min_r-BUFSIZE, 1)
   bb%i_ulat = min(max_r+BUFSIZE, INR)
   bb%i_llon = max(min_c-BUFSIZE, 1)
   bb%i_ulon = min(max_c+BUFSIZE, INC)
   bb%NLAT = bb%i_ulat - bb%i_llat + 1
   bb%NLON = bb%i_ulon - bb%i_llon + 1
end function find_bounding_box_bilinear

!BOP
!
! !ROUTINE: find_bounding_box_budget_bilinear
! \label{find_bounding_box_budget_bilinear}
!
! !REVISION HISTORY:
! 14 Aug 2025: James Geiger, Initial Specification
!
! !INTERFACE:
function find_bounding_box_budget_bilinear(INC, INR, ilat, ilon, ONCxONR, n111, n121, n211, n221, n112, n122, n212, n222) result(bb)
! !USES:
! none

   implicit none

! !ARGUMENTS:
   integer, intent(in) :: INC, INR, ONCxONR
   real, dimension(INR), intent(in) :: ilat
   real, dimension(INC), intent(in) :: ilon
   integer, dimension(ONCxONR), intent(in) :: n111, n121, n211, n221
   integer, dimension(ONCxONR,25), intent(in) :: n112, n122, n212, n222

! !DESCRIPTION:
! This function finds the bounding box within the input domain
! that contains the given output domain.  This function returns
! a bounding box datatype.
!
! The arguments are:
! \begin{description}
! \item[INC]
!    number of columns in the input domain
! \item[INR]
!    number of rows in the input domain
! \item[ilat]
!    array of latitude values for the input domain
! \item[ilon]
!    array of longitude values for the input domain
! \item[ONCxONR]
!    maximum number of grid-cells in the output domain
!    (number of columns times number of rows)
! \item[n111]
!    n111 array from bilinear_interp_input
! \item[n121]
!    n121 array from bilinear_interp_input
! \item[n211]
!    n211 array from bilinear_interp_input
! \item[n221]
!    n221 array from bilinear_interp_input
! \item[n112]
!    n112 array from conserv_interp_input
! \item[n122]
!    n122 array from conserv_interp_input
! \item[n212]
!    n212 array from conserv_interp_input
! \item[n222]
!    n222 array from conserv_interp_input
! \end{description}
!
!EOP

   integer, parameter :: BUFSIZE = 0
   type(bounding_box_type) :: bb

   integer :: min_c, max_c, min_r, max_r
   integer :: min_r111, max_r111, min_c111, max_c111
   integer :: min_r121, max_r121, min_c121, max_c121
   integer :: min_r211, max_r211, min_c211, max_c211
   integer :: min_r221, max_r221, min_c221, max_c221
   integer :: min_r112, max_r112, min_c112, max_c112
   integer :: min_r122, max_r122, min_c122, max_c122
   integer :: min_r212, max_r212, min_c212, max_c212
   integer :: min_r222, max_r222, min_c222, max_c222

   integer, allocatable, target, dimension(:) :: tmp_array

   tmp_array = remove_zeros(ONCxONR, n111)
   tmp_array = tmp_array - 1
   min_c111 = minval(mod(tmp_array, INC) + 1)
   max_c111 = maxval(mod(tmp_array, INC) + 1)
   min_r111 = minval(tmp_array / INC + 1)
   max_r111 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   tmp_array = remove_zeros(ONCxONR, n121)
   tmp_array = tmp_array - 1
   min_c121 = minval(mod(tmp_array, INC) + 1)
   max_c121 = maxval(mod(tmp_array, INC) + 1)
   min_r121 = minval(tmp_array / INC + 1)
   max_r121 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   tmp_array = remove_zeros(ONCxONR, n211)
   tmp_array = tmp_array - 1
   min_c211 = minval(mod(tmp_array, INC) + 1)
   max_c211 = maxval(mod(tmp_array, INC) + 1)
   min_r211 = minval(tmp_array / INC + 1)
   max_r211 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   tmp_array = remove_zeros(ONCxONR, n221)
   tmp_array = tmp_array - 1
   min_c221 = minval(mod(tmp_array, INC) + 1)
   max_c221 = maxval(mod(tmp_array, INC) + 1)
   min_r221 = minval(tmp_array / INC + 1)
   max_r221 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   tmp_array = remove_zeros(ONCxONR*25, reshape(n112,(/ONCxONR*25/)))
   tmp_array = tmp_array - 1
   min_c112 = minval(mod(tmp_array, INC) + 1)
   max_c112 = maxval(mod(tmp_array, INC) + 1)
   min_r112 = minval(tmp_array / INC + 1)
   max_r112 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   tmp_array = remove_zeros(ONCxONR*25, reshape(n122,(/ONCxONR*25/)))
   tmp_array = tmp_array - 1
   min_c122 = minval(mod(tmp_array, INC) + 1)
   max_c122 = maxval(mod(tmp_array, INC) + 1)
   min_r122 = minval(tmp_array / INC + 1)
   max_r122 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   tmp_array = remove_zeros(ONCxONR*25, reshape(n212,(/ONCxONR*25/)))
   tmp_array = tmp_array - 1
   min_c212 = minval(mod(tmp_array, INC) + 1)
   max_c212 = maxval(mod(tmp_array, INC) + 1)
   min_r212 = minval(tmp_array / INC + 1)
   max_r212 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   tmp_array = remove_zeros(ONCxONR*25, reshape(n222,(/ONCxONR*25/)))
   tmp_array = tmp_array - 1
   min_c222 = minval(mod(tmp_array, INC) + 1)
   max_c222 = maxval(mod(tmp_array, INC) + 1)
   min_r222 = minval(tmp_array / INC + 1)
   max_r222 = maxval(tmp_array / INC + 1)

   min_r = min(min_r111, min_r121, min_r211, min_r221, min_r112, min_r122, min_r212, min_r222)
   max_r = max(max_r111, max_r121, max_r211, max_r221, max_r112, max_r122, max_r212, max_r222)
   min_c = min(min_c111, min_c121, min_c211, min_c221, min_c112, min_c122, min_c212, min_c222)
   max_c = max(max_c111, max_c121, max_c211, max_c221, max_c112, max_c122, max_c212, max_c222)

   bb%i_llat = min_r
   bb%i_ulat = max_r
   bb%i_llon = min_c
   bb%i_ulon = max_c
   bb%i_llat = max(min_r-BUFSIZE, 1)
   bb%i_ulat = min(max_r+BUFSIZE, INR)
   bb%i_llon = max(min_c-BUFSIZE, 1)
   bb%i_ulon = min(max_c+BUFSIZE, INC)
   bb%NLAT = bb%i_ulat - bb%i_llat + 1
   bb%NLON = bb%i_ulon - bb%i_llon + 1
end function find_bounding_box_budget_bilinear

!BOP
!
! !ROUTINE: find_bounding_box_average
! \label{find_bounding_box_average}
!
! !REVISION HISTORY:
! 17 Jul 2025: James Geiger, Initial Specification
!
! !INTERFACE:
function find_bounding_box_average(INC, INR, ilat, ilon, LNC, LNR, olat, olon, n111) result(bb)
! !USES:
! none

   implicit none

! !ARGUMENTS:
   integer, intent(in) :: INC, INR, LNC, LNR
   real, dimension(INR), intent(in) :: ilat
   real, dimension(INC), intent(in) :: ilon
   real, dimension(LNC*LNR), intent(in) :: olat
   real, dimension(LNC*LNR), intent(in) :: olon
   integer, dimension(INC*INR), intent(in) :: n111

! !DESCRIPTION:
! This function finds the bounding box within the input domain
! that contains the given output domain.  This function returns
! a bounding box datatype.
!
! The arguments are:
! \begin{description}
! \item[INC]
!    number of columns in the input domain
! \item[INR]
!    number of rows in the input domain
! \item[ilat]
!    array of latitude values for the input domain
! \item[ilon]
!    array of longitude values for the input domain
! \item[LNC]
!    number of columns in the output domain
! \item[LNR]
!    number of rows in the output domain
! \item[olat]
!    array of latitude values for the output domain
! \item[olon]
!    array of longitude values for the output domain
! \item[n111]
!    n111 array from upscaleByAveraging_input
! \end{description}
!
!EOP

   integer, parameter :: BUFSIZE = 0
   type(bounding_box_type) :: bb

   integer :: min_c, max_c, min_r, max_r
   integer :: min_r111, max_r111, min_c111, max_c111
   real :: min_lat, max_lat, min_lon, max_lon

   integer :: tmp_size, i, c
   real, allocatable, dimension(:) :: tmp_array_r
   integer, dimension(1) :: ml, fl

   tmp_size = count(n111 > 0)
   allocate(tmp_array_r(tmp_size))
   c = 1
   do i = 1, INC*INR
      if ( n111(i) > 0 ) then
         tmp_array_r(c) = ilon(mod((i - 1), INC) + 1)
         c = c + 1
      endif
   enddo
   min_lon = minval(tmp_array_r)
   max_lon = maxval(tmp_array_r)

   c = 1
   do i = 1, INC*INR
      if ( n111(i) > 0 ) then
         tmp_array_r(c) = ilat((i - 1) / INC + 1)
         c = c + 1
      endif
   enddo
   min_lat = minval(tmp_array_r)
   max_lat = maxval(tmp_array_r)
   deallocate(tmp_array_r)

   min_r = find_lower(INR, ilat, min_lat)
   max_r = find_upper(INR, ilat, max_lat) 
   min_c = find_lower(INC, ilon, min_lon)
   max_c = find_upper(INC, ilon, max_lon)

   bb%i_llat = min_r
   bb%i_ulat = max_r
   bb%i_llon = min_c
   bb%i_ulon = max_c
   bb%i_llat = max(min_r-BUFSIZE, 1)
   bb%i_ulat = min(max_r+BUFSIZE, INR)
   bb%i_llon = max(min_c-BUFSIZE, 1)
   bb%i_ulon = min(max_c+BUFSIZE, INC)
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
end function find_bounding_box_average

!BOP
!
! !ROUTINE: find_bounding_box_neighbor
! \label{find_bounding_box_bilinear}
!
! !REVISION HISTORY:
! 5 Sep 2025: James Geiger, Initial Specification
!
! !INTERFACE:
function find_bounding_box_neighbor(INC, INR, ilat, ilon, ONCxONR, n113) result(bb)
! !USES:
! none

   implicit none

! !ARGUMENTS:
   integer, intent(in) :: INC, INR, ONCxONR
   real, dimension(INR), intent(in) :: ilat
   real, dimension(INC), intent(in) :: ilon
   integer, dimension(ONCxONR), intent(in) :: n113

! !DESCRIPTION:
! This function finds the bounding box within the input domain
! that contains the given output domain.  This function returns
! a bounding box datatype.
!
! The arguments are:
! \begin{description}
! \item[INC]
!    number of columns in the input domain
! \item[INR]
!    number of rows in the input domain
! \item[ilat]
!    array of latitude values for the input domain
! \item[ilon]
!    array of longitude values for the input domain
! \item[ONCxONR]
!    maximum number of grid-cells in the output domain
!    (number of columns times number of rows)
! \item[n113]
!    n113 array from neighbor_interp_input
! \end{description}
!
!EOP

   integer, parameter :: BUFSIZE = 0
   type(bounding_box_type) :: bb

   integer :: min_c, max_c, min_r, max_r
   integer :: min_r113, max_r113, min_c113, max_c113

   integer, allocatable, target, dimension(:) :: tmp_array

   tmp_array = remove_zeros(ONCxONR, n113)
   tmp_array = tmp_array - 1
   min_c113 = minval(mod(tmp_array, INC) + 1)
   max_c113 = maxval(mod(tmp_array, INC) + 1)
   min_r113 = minval(tmp_array / INC + 1)
   max_r113 = maxval(tmp_array / INC + 1)
   deallocate(tmp_array)

   min_r = min_r113
   max_r = max_r113
   min_c = min_c113
   max_c = max_c113

   bb%i_llat = min_r
   bb%i_ulat = max_r
   bb%i_llon = min_c
   bb%i_ulon = max_c
   bb%i_llat = max(min_r-BUFSIZE, 1)
   bb%i_ulat = min(max_r+BUFSIZE, INR)
   bb%i_llon = max(min_c-BUFSIZE, 1)
   bb%i_ulon = min(max_c+BUFSIZE, INC)
   bb%NLAT = bb%i_ulat - bb%i_llat + 1
   bb%NLON = bb%i_ulon - bb%i_llon + 1
end function find_bounding_box_neighbor
end module bounding_box_mod
