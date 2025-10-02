module lista
implicit none
integer, parameter :: f=6000, g=685035! f para cosmos.dat, g para galaxies.dat
real, parameter :: pi=3.1415
!si cambiamos la long de la tabla tenemos q cambiar f a mano, no hay forma de automatizarlo
real, dimension(f) :: redshift, coord_comovil, dist_luminosa, dist_angular
real::  z, chi, dL,  dLmax,zmax
end module lista


program limpieza 
use lista
implicit none
integer ::  i, j, k
real ::  petro_r, red_r , r50, corr_ab=0.010
real ::  petro_s_ext, petro_abs, vmax
real :: distL, red, vol, zmin=0.
real, dimension(5):: rk_p, rks_p

external distL, red, vol

open(unit=25,file='cosmos.dat',status='old')
read(25,*)  ! Salta la primera línea (header)
do i = 2 , f
read(25,*) redshift(i), coord_comovil(i), dist_luminosa(i), dist_angular(i) 
end do 
close(25)

open(unit=30,file='datos.dat',status='old')
do i=1,4
    read(30,*)z,petro_r,red_r ,r50,(rk_p(k),k=1,5),(rks_p(k),k=1,5)
    !red == reddening o extinción
    !
    if (z>0.15) cycle
    if (r50<1.5) cycle
    petro_s_ext = petro_r -red_r
    if (petro_s_ext<14.5 .OR. petro_s_ext>17.77) cycle

    dL=distL(z)

    petro_abs = petro_s_ext-5*log10(dL)-25+corr_ab + rks_p(3)
    
    dLmax= 10.**((17.77-petro_abs)/5-5)   
    zmax=red(dLmax)
    call qromb(vol,zmin,zmax,vmax)

    write(35,*)z,',',petro_s_ext,',',petro_abs,',',rks_p(3),',',zmax,',', vmax
    !write(35,*)vmax




end do
close(30)

end program limpieza


!_______________________________________FUNCION INTERPOLACION___________________________________________________
! dL(z)
function distL(x)
use lista
implicit none
real::x, distL
integer :: i, j
call locate(redshift,f,x,j)
distL= ((dist_luminosa(j+1)-dist_luminosa(j))/(redshift(j+1)-redshift(j)))*(x-redshift(j)) + dist_luminosa(j)

end function distL
! z(dL)
function red(x)
use lista
implicit none
real::x, red
integer :: i, j
call locate(dist_luminosa,f,x,j)
red= ((redshift(j+1)-redshift(j))/(dist_luminosa(j+1)-dist_luminosa(j)))*(x-dist_luminosa(j)) + redshift(j)
end function red
!_______________________________________FUNCION INTEGRANDO___________________________________________________

function vol(x)
!x==z

implicit none
real, parameter :: Om = 0.3, ODark =0.7 
real :: vol,x, distL
external distL

vol=(distL(x)**2)/(((1+x)**2)*sqrt(Om*((1+x)**3)+ODark))

end function vol

include 'locate.f'
include 'qromb.f'
include 'trapzd.f'
include 'polint.f'