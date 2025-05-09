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
! !ROUTINE: AGRMET_readpcpcntm
!  \label{AGRMET_readpcpcntm}
! 
! !REVISION HISTORY:
!
!    1  nov 05  Sujay Kumar, Initial specification
!   19  Dec 07  Marv Freimund, Simplify filename creation
!   11  Mar 10  Chris Franks, Changed program names in messages to LIS.
!   29  May 15  Ryan Ruhge, Updated file reads
!
! !INTERFACE:
subroutine AGRMET_readpcpcntm(n)
! !USES: 
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_coreMod,       only : LIS_rc, LIS_ews_halo_ind, LIS_ewe_halo_ind,&
                                LIS_nss_halo_ind, LIS_nse_halo_ind, LIS_localPet
  use LIS_logMod,        only : LIS_logunit, LIS_abort, LIS_endrun, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use AGRMET_forcingMod, only : agrmet_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This routine reads the spreading radii used in the barnes analyses
!  for AGRMET precipitation processing. The data is in polar stereographic
!  projection at 48km resolution. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!   
!  The routines invoked are: 
!  \begin{description}
!   \item[lis\_abort](\ref{LIS_abort}) \newline
!    program aborts in case of error in file open or read
!   \end{description}
!
!EOP

  logical       :: exists
  character(len=LIS_CONST_PATH_LEN) :: name
  character(len=LIS_CONST_PATH_LEN) :: message(20)
  integer       :: ftn
  real          :: data_in(LIS_rc%gnc(n), LIS_rc%gnr(n))
  integer       :: istat

  name = trim(agrmet_struc(n)%climodir)//'/pcp_spread_radii.1gd4r'
  inquire( file = trim(name) , exist = exists)
  if ( exists ) then
     !write(LIS_logunit,*)' '
     write(LIS_logunit,*)'[INFO] READING ', trim(name)
     ftn= LIS_getNextUnitNumber()
     open(ftn, file=trim(name), access='direct',&
          status='old', form="unformatted", recl=LIS_rc%gnr(n)*LIS_rc%gnc(n)*4)
     
     read(ftn, rec=1, iostat=istat) data_in
     agrmet_struc(n)%radius(:,:) =  data_in(&
        LIS_ews_halo_ind(n,LIS_localPet+1):&
        LIS_ewe_halo_ind(n,LIS_localPet+1), &
        LIS_nss_halo_ind(n,LIS_localPet+1): &
        LIS_nse_halo_ind(n,LIS_localPet+1))

     call LIS_releaseUnitNumber(ftn)
     
  else
     write(LIS_logunit,*) '[ERR] LIS:  error opening file ',trim(name)
     write(LIS_logunit,*) '[ERR]  file does not exist'
     write(LIS_logunit,*) '[ERR]  this could indicate serious data discrepancies'
     write(LIS_logunit,*) '[ERR]  LIS will be aborted'
     message(1) = 'program:  LIS'
     message(2) = '  routine: AGRMET_readpcpcntm'
     message(3) = '  error opening file '//trim(name)
     message(4) = '  file does not exist'
     message(5) = '  this could indicate serious data discrepancies'
     call LIS_abort( message )
     call LIS_endrun
  endif

end subroutine AGRMET_readpcpcntm

