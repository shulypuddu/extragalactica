!::MODULO:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
module lista
implicit none
integer :: i , j, k! e= eliptica, s= espiral, l= lenticular
integer, parameter :: long=5725 , bin=8 ! long cantidad de galaxias , bin = cant de bines 
real, parameter :: rad=4.*atan(1.)/180, pi=4.*atan(1.) , H0=70.,c=299792.
!-0.651179075
!2.8375802
end module lista

!::PROGRAMA:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
program hist
use lista
use numerical
implicit none
real::  logdens, rho, ty, dr, rho_medio
real, parameter:: rho_min=-0.651179075, rho_max=2.8375802
real, dimension(bin)::  e, s, l
integer, dimension(bin):: total
dr = (rho_max-rho_min)/real(bin) 


open(unit=15, file='guardar.dat',status='old')
open(unit=25, file='bines.dat',status='unknown')
read(15,*)
e=0.
s=0.
l=0.
do i =1,long
read(15,*)logdens,rho,ty
j= int((logdens-rho_min)/dr)+1

if (ty==4) cycle

total(j)=total(j)+1

if(ty==1) s(j)=s(j)+1
if (ty==2) l(j)=l(j)+1
if(ty==3) e(j)=e(j) +1
end do

do k=1,bin
    e(k)=e(k)/real(total(k))
    l(k)=l(k)/real(total(k))
    s(k)=s(k)/real(total(k))
    rho_medio = rho_min+dr*real(2*k-1)*0.5   

write(25,*) rho_medio, total(k),e(k),l(k),s(k)
end do 

close(25)
close(15)
end program hist
