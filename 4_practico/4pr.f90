! gfortran -o dressler nrtype.f90 nrutil.f90 numerical.f90 4pr.f90
! al correrlo no solo corro mi archivo sino todos los que contiene los modulos externos
!::MODULOS:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
module lista
implicit none
integer :: i , j, k, x, m, n, num,ig 
integer, parameter :: l=56  ! l cantidad de cumulos 
real, parameter :: rad=4.*atan(1.)/180, pi=4.*atan(1.) , H0=70.,c=299792.
real, allocatable, dimension(:) :: ra, dec, r_mpc ! tipo morfologico
integer, allocatable, dimension(:) :: indice,ty!para la subrutina index_sp
real:: cls, ra_0, dec_0,ra_1, dec_1, y,redshift_cls, dist

end module lista

!::PROGRAMA:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
program dressler
use lista
use numerical
implicit none
real:: ra_d, dec_d, logdens, denso, rho, theta, tita, d10
external rho, theta

open(unit=15, file='centros.dat',status='old')
open(unit=25, file='galaxias.dat',status='old')
open(unit=35,file='guardar.dat',status='unknown')
read(15,*)
read(25,*)
write(35,*)'logrho,rho,tipo'

do i = 1, l !voy cúmulo a cúmulo
    read(15,*) cls, num, redshift_cls,y,y
    dist= redshift_cls*c/H0
    allocate(ra(num),dec(num),r_mpc(num),ty(num),indice(num))

    do j=1,num !necesitamos leer TODAS las galaxias antes de calcular las cosas
        read(25,*) x,x,ra_d,dec_d,ty(j) 
        ra(j)=ra_d*rad
        dec(j)=dec_d*rad
    end do
    
    r_mpc=99999 ! para que si lo saltea, al ordenar de menor a mayor queda afuera.
    do k=1,num 
        ra_0= ra(k) !ra_0 y dec_0 para ubicar
        dec_0=dec(k)
        do j=1,num
            if(k==j) cycle
            ra_1= ra(j) !ra_1 y dec_1 para recorrer
            dec_1=dec(j)
            r_mpc(j) = theta(ra_0,dec_0,ra_1,dec_1,dist) !distancia ya en mpc 
        end do  !do q calcula la distancias de una galaxia con el resto
      !do q selecciona cada galaxia
    call indexx_sp(r_mpc,indice)
    d10=r_mpc(indice(10))
    denso= rho(d10)
    logdens = log10(denso)
    write(35,*) logdens,',',denso,',',ty(k)
    end do
    deallocate(ra,dec,r_mpc,ty,indice)
end do !do que recorre cada cumulo

close(35)
close(25)
close(15)
end program dressler

!::FUNCIONES:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
function rho(r10) !calculo la densidad de galaxias por unidad de area
use lista
real:: rho , r10
rho = 10 / (pi*r10**2)
end function rho

function theta(ra0,dec0,ra1,dec1, r) !calculo la distancia EN MPC entre UN PAR DE GALAXIAS
use lista 
implicit none
real:: ra0, dec0, ra1, dec1, cos_theta, theta, r
cos_theta=  sin(dec1)*sin(dec0) + cos(dec1)*cos(dec0)*cos(ra1-ra0)
theta=tan(acos(cos_theta))*r 
return
end function theta
