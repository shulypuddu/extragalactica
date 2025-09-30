module lista
implicit none
integer, parameter :: f=6000, g=685035, pi=3.1415 ! f para cosmos.dat, g para galaxies.dat
!si cambiamos la long de la tabla tenemos q cambiar f a mano, no hay forma de automatizarlo
real, dimension(f) :: redshift, coord_comovil, dist_luminosa, dist_angular
real::  z, chi, dL
end module lista


program limpieza 
use lista
implicit none
integer ::  i, j, k
real ::  petro_r, red_r , r50, corr_ab=0.010
real ::  petro_s_ext, petro_abs, coord, Vmax
real, dimension(5):: rk_p, rks_p


external coord

open(unit=25,file='cosmos.dat',status='unknown')
read(25,*)  ! Salta la primera línea (header)
do i = 2 , f
read(25,*) redshift(i), coord_comovil(i), dist_luminosa(i), dist_angular(i) 
end do 
close(25)

open(unit=30,file='datos.dat',status='old')

do i=1,g
    read(30,*)z,petro_r,red_r ,r50,(rk_p(k),k=1,5),(rks_p(k),k=1,5)
    !red == reddening o extinción
    !
    if (z>0.15) cycle
    if (r50<1.5) cycle
    petro_s_ext = petro_r -red_r
    if (petro_s_ext<14.5 .OR. petro_s_ext>17.77) cycle

    chi=coord(z)
    dL= (1.+z)*chi
    petro_abs = petro_s_ext-5*log10(dL)-25+corr_ab
    

    write(35,*)z,',',petro_s_ext,',',petro_abs,',',rks_p(3),',',
end do
close(30)

end program limpieza


!_______________________________________FUNCION INTERPOLACION___________________________________________________
function coord(x)
use lista
implicit none
real::x, coord
integer :: i, j
call locate(redshift,f,x,j)
coord= ((coord_comovil(j+1)-coord_comovil(j))/(redshift(j+1)-redshift(j)))*(x-redshift(j)) + coord_comovil(j)

end function coord

include 'locate.f'
