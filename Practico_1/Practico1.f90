program cosmo100

implicit none
real:: dl, da,  dchi , chi
real:: H_0=100., c=299792., z
integer :: i 
external dchi


open(unit=25,file='cosmos.csv',status='unknown')

write(25,*) 'z',',', 'coord_comovil',',', 'dist_luminosa',',', 'dist_angular '
do i = 1 , 100 
z = z +0.001*i

call qromb(dchi,0.,z,chi)
write(25,*) z ,',',  chi*c/H_0 ,',', chi*(1+z)*c/H_0 ,',', (chi*c)/((1+z)*H_0)
end do

close(25)


end program cosmo100

!____________FUNCION INTEGRANDDO____________________________________________________
function dchi(z) 
implicit none
real:: z, dchi
real, parameter::  o_mat=0.3 , o_lambda=0.7
dchi=1./sqrt(o_mat*(1.+z)**3+o_lambda)
end function dchi

!____________SUBRUTINAS INCLUIDAS___________________________________________________
include 'qromb.f' 
include 'polint.f' 
include 'trapzd.f'
!include 'locate.f'
