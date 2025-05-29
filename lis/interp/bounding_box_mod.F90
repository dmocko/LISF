module bounding_box_mod
type bounding_box_type
   integer :: i_llat ! lower latitude index
   integer :: i_ulat ! upper latitude index
   integer :: i_llon ! lower longitude index
   integer :: i_ulon ! upper longitude index
   integer :: NLAT ! number of latitude points
   integer :: NLON ! number of longitude points
end type

contains

function find_bounding_box(NC, NR, ilat, ilon, min_olat, max_olat, min_olon, max_olon) result(bb)
   implicit none

   integer :: NC, NR
   real, dimension(NC) :: ilon
   real, dimension(NR) :: ilat
   real :: min_olat, max_olat, min_olon, max_olon
   type(bounding_box_type) :: bb

   if ( min_olat < ilat(1) .or. &
        max_olat > ilat(NR) .or. &
        min_olon < ilon(1) .or. &
        max_olon > ilon(NC) ) then
        print*, '[ERR] output domain is outside input domain'
        stop 666
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
