program galaxies
implicit none
integer :: i, j
integer,parameter :: f=100, g=20000 ! f para cosmos.dat, g para galaxies.dat
!si cambiamos la long de la tabla tenemos q cambiar f a mano, no hay forma de automatizarlo
character(len=18),dimension(g) :: ID 
!guardamos los ID de las galaxias como texto xq tienen una longitud de caracteres mayor a la de un integer y es una paja usar un real
!real :: ra,dec,z=0.125
real, dimension(f) :: redshift,coord_comovil,dist_luminosa,dist_angular


open (unit=25,file='cosmos.dat',status='old')


do i = 1 , f
read(25,*) redshift(i), coord_comovil(i), dist_luminosa(i), dist_angular(i) 
end do 
close(25)

call locate(redshift,f,z,j)
print *, 'El indice j es: ', j
print *, redshift(j), coord_comovil(j), dist_luminosa(j), dist_angular(j)

open (unit=30,file='galaxies.dat',status='old')
read(30,*) ID, ra, dec, 







end program galaxies

include 'locate.f'

