!::MODULOS:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
module lista
implicit none
integer :: i , j, k, x, m, n, num,ty,ig
integer, parameter :: l=56  ! l cantidad de cumulos 
real, parameter :: pi2= 1.57079632, rad=0.0174533, pi=3.14159265359, H0=70.
real, allocatable, dimension(:) :: ra, dec, r_mpc
real:: cls, ra_0, dec_0,ra_1, dec_1, y,redshift_cls, dist


end module lista

!::PROGRAMA:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
program practico
use lista
implicit none
real:: ra_d, dec_d, r_rad, rho, theta, tita
external rho, theta

open(unit=15, file='centros.dat',status='old')
open(unit=25, file='galaxias.dat',status='old')
open(unit=35,file='guardar.dat',status='unknown')
read(15,*)
read(25,*)

do i = 1, l !voy cúmulo a cúmulo
read(15,*) cls, num, redshift_cls,y,y
dist= H0*redshift_cls
allocate(ra(num),dec(num),r_mpc(num))

do j=1,num !necesitamos tener leidos TODAS las galaxias antes de calcular las cosas
read(25,*) x,x,ra_d,dec_d,ty 
ra(j)=ra_d*rad
dec(j)=dec_d*rad
end do

do k=1,num 
ra_0= ra(k) !ra_0 y dec_0 para ubicar
dec_0=dec(k)
    do j=1,num
    if(k==j) cycle
    ra_1= ra(j)
    dec_1=dec(j)
    tita = theta(ra_0,dec_0,ra_1,dec_1)
    r_mpc(i) = tan(tita)*dist ! valeria nos dijo q no hace falta usar cosmos, verificar q si esta bien esta forma de pasar las distancias
    end do  !do q calcula la distancias de una galaxia con el resto
end do  !do q selecciona cada galaxia

deallocate(ra,dec,r_mpc)
end do !do que recorre cada cumulo

close(35)
close(25)
close(15)
end program practico

!::FUNCIONES:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
function rho(r10) !calculo la densidad de galaxias por unidad de area
use lista
real:: rho , r10
rho = 10 / (pi*r10**2)
end function rho

function theta(ra0,dec0,ra1,dec1) !calculo la distancia ANGULAR entre UN PAR DE GALAXIAS
use lista 
implicit none
real:: ra0, dec0, ra1, dec1, cos_theta, theta
cos_theta=  cos(pi2-dec(i))*cos(pi2-dec0) - sin(pi2-dec(i))*sin(pi2-dec0)*cos(ra(i)-ra0)
theta= acos(cos_theta)
return
end function theta

!::SUBRUTINAS:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
