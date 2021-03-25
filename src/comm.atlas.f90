MODULE ATLCOMM
  use types

  implicit none
!

!
!.... 1996 Apr - CHANGED maxmol, maxmeq, maxloc TO CONFORM TO
!
!                BOB'S VALUES IN XNFPELSYN DATED 3Jan95
!
!.... maxd   = maximum number of depths in the atmosphere
!.... maxmeq = maximum number of molecular equations
!.... maxloc = maximum number of molecular components
!.... maxmol = maximum number of molecules
!.... maxmu  = maximum number of angles for the radiation field
!.... maxnu  = the maximum number of frequencies

  integer, parameter :: maxd = 360
  integer, parameter :: maxloc = 600
  integer, parameter :: maxmeq = 30
  integer, parameter :: maxmol = 200
  integer, parameter :: maxmu = 35

! for control
      character(30) fileODFinput, filemodelinput, filefluxinput
      character(30) filePTgrid, filestartmodel, filefluxmodel
      character(30) fodfodf, fmodelodf, ffluxodf
!----- names for the current input to consider: 
      character(30) INPUTfile, modelfile, ODFfile


!  for frequency grid (bin grid)
!      common /freset/  freset, rcoset, inifreset, nulo, nuhi, numnu
  integer   isubbin, ibin 
  integer    numt, numpres 

  integer  nulo, nuhi, numnu
  integer, allocatable :: ifpnch(:)
  integer(kind=2), allocatable :: kapw(:,:,:,:), kapwvar(:,:,:,:,:)
  real(kind=dp), allocatable  ::  freset(:), rcoset(:), inifreset(:)

!-- for communication / parallel calc model

  real(kind=dp), allocatable :: sdarrays(:,:), sendscalar(:)
  real(kind=dp), allocatable :: recarrays(:,:), recscalar(:)
 
!  for ODF's T and P gird
!      common / tabptv / tabkap, tabp, tabt, tabv, iftabk, nt, np, nv

  real(kind=dp), allocatable::   tabp(:), tabt(:)
  real(kind=dp), allocatable:: sbvalues(:,:), sbwt(:,:) ! boarders of subbins, with of subbin
  integer(kind=2), allocatable:: tabrosskap(:,:)


!  for parallel flux calculations 

  real(kind=dp), allocatable :: glwave(:), allwave(:), glsurfin(:,:), allsurfin(:,:) 
  logical :: freqsplit   

!-------------------------------------------------------------------------
!  for pressure calculations
!------------------------------------------------------------------------
!
!      common /pzerob/  knu, pcon, pradk, pradk0, pturb0, pzero, raden
  real(kind=dp)  knu(maxd), pcon, pradk(maxd), pradk0, pturb0
  real(kind=dp)  pzero, raden(maxd)
  real(kind=dp)  ptotal(maxd)
  real(kind=dp)  accrad(maxd), prad(maxd)

!--- for model calculations, radiap and tcorr routines:

  real(kind=dp) :: hflux(maxd), rdabh(maxd), rjmins(maxd), rdiagj(maxd)


!-------------------------------------------------------------------



!-----------------------------------------------------------------------


CONTAINS


subroutine set_freq(isub, ifre, nt, np, ivar) 

  implicit none

  integer, intent(in) :: ivar,  isub, ifre, nt, np

  allocate(freset(ifre))
  allocate(rcoset(ifre))
  allocate(inifreset(ifre+1))
  
  if(ivar .eq.0)  allocate(kapw(isub, ifre,np, nt ) )
  if(ivar .eq. 1) allocate(kapwvar(isub, ifre,np, nt, 5))

!------------------------
  allocate(tabp(np))
  allocate(tabt(nt))

  allocate(sbvalues(isub,ifre))
  allocate(sbwt(isub,ifre))

!------------------------

end subroutine set_freq

subroutine close_freq
  implicit none

  deallocate(freset)
  deallocate(rcoset)
  deallocate(inifreset)

  if(allocated(kapw))  deallocate(kapw)
  if(allocated(kapwvar)) deallocate(kapwvar)


!--------------------
  deallocate(tabp)
  deallocate(tabt)

  deallocate(sbvalues)
  deallocate(sbwt)
  
end subroutine close_freq


subroutine alloc_comm(numa) 
 implicit none
 integer, intent(in) :: numa
 integer :: ierr(3)
 integer :: numb
 numb = 9  
 ierr = -1 
 allocate( sdarrays(numa,numb), stat = ierr(1)  ) 
 allocate( recarrays(numa,numb), stat  = ierr(2) ) 
 allocate( sendscalar(numb), stat  = ierr(3)  ) 

 if (maxval(ierr) .gt. 0) then 
   print*, 'array allocation in alloc_comm failed'
   stop
 end if 


end subroutine alloc_comm

subroutine dealloc_comm
  implicit none 

  if(allocated(sdarrays)) deallocate(sdarrays) 
  if(allocated(recarrays)) deallocate(recarrays)
  if(allocated(sendscalar)) deallocate(sendscalar)


end subroutine dealloc_comm 


subroutine set_splitf(nf, numa) 
  implicit none
  integer, intent(in) :: numa, nf
  
  allocate(glsurfin(nf, numa))
  allocate(allsurfin(nf, numa))
  allocate(glwave(nf))
  allocate(allwave(nf))

  glsurfin = 0.0d0
  glwave = 0.0d0 
end subroutine set_splitf


subroutine close_splitf
  implicit none
  
  if(allocated(glsurfin)) deallocate(glsurfin)
  if(allocated(allsurfin)) deallocate(allsurfin)
  if(allocated(glwave)) deallocate(glwave)
  if(allocated(allwave)) deallocate(allwave)

end subroutine close_splitf




END MODULE ATLCOMM

