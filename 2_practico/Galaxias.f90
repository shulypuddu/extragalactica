! filepath: /mnt/sda2/extragalactica/2_practico/pr2.f90
module lista
implicit none
integer, parameter :: f=100, g=20000 ! f para cosmos.dat, g para galaxies.dat
!si cambiamos la long de la tabla tenemos q cambiar f a mano, no hay forma de automatizarlo
real, dimension(f) :: redshift, coord_comovil, dist_luminosa, dist_angular
real, dimension(g) ::  z, chi, dL
end module lista

program galaxies
use lista
implicit none
integer :: i, j
character(len=18),dimension(g) :: ID 
!guardamos los ID de las galaxias como texto xq tienen una longitud de caracteres mayor a la de un integer y es una paja usar un real
real :: coord0, dL0, coord, z0

real, dimension(g) :: ra, dec,fracDeV_r, velDisp
real, dimension(g) :: petroMag_u, petroMag_g, petroMag_r, petroMag_i, petroMag_z
real, dimension(g) :: modelMag_u, modelMag_g, modelMag_r, modelMag_i, modelMag_z
real, dimension(g) :: extinction_u, extinction_g, extinction_r, extinction_i, extinction_z
real, dimension(g) :: petroRad_u, petroRad_g, petroRad_r, petroRad_i, petroRad_z
real, dimension(g) :: petroR50_u, petroR50_g, petroR50_r, petroR50_i, petroR50_z
real, dimension(g) :: abs_petroMag_u, abs_petroMag_g, abs_petroMag_r, abs_petroMag_i, abs_petroMag_z
real, dimension(g) :: abs_modelMag_u, abs_modelMag_g, abs_modelMag_r, abs_modelMag_i, abs_modelMag_z

external coord

open (unit=25,file='cosmos.dat',status='old')
read(25,*)  ! Salta la primera línea (header)
do i = 2 , f
read(25,*) redshift(i), coord_comovil(i), dist_luminosa(i), dist_angular(i) 
end do 
close(25)

open (unit=30,file='tabla.dat',status='old')
read(30,*)  ! Salta la primera línea (header)
do i = 1 , g
  read(30,*) ID(i),ra(i),dec(i),z(i),fracDeV_r(i),velDisp(i), &
    petroMag_u(i),petroMag_g(i),petroMag_r(i),petroMag_i(i),petroMag_z(i), &
    modelMag_u(i),modelMag_g(i),modelMag_r(i),modelMag_i(i),modelMag_z(i), &
    extinction_u(i),extinction_g(i),extinction_r(i),extinction_i(i),extinction_z(i), &
    petroRad_u(i),petroRad_g(i),petroRad_r(i),petroRad_i(i),petroRad_z(i), &
    petroR50_u(i),petroR50_g(i),petroR50_r(i),petroR50_i(i),petroR50_z(i)
end do
close(30)

do i = 1 , g
chi(i)=coord(z(i))
dL(i)= (1+z(i))*chi(i)
abs_petroMag_u(i)=petroMag_u(i)-5*log10(dL(i))-25-extinction_u(i)-0.036
abs_petroMag_g(i)=petroMag_g(i)-5*log10(dL(i))-25-extinction_g(i)+0.012
abs_petroMag_r(i)=petroMag_r(i)-5*log10(dL(i))-25-extinction_r(i)+0.010
abs_petroMag_i(i)=petroMag_i(i)-5*log10(dL(i))-25-extinction_i(i)+0.028
abs_petroMag_z(i)=petroMag_z(i)-5*log10(dL(i))-25-extinction_z(i)+0.040
end do




print *, ID(2),ra(2), dec(2),z(2), dL(2), abs_petroMag_u(2)


!open (unit=35,file='galaxias.dat',status='old')
!do i = 1, g
!write(35,*) ' ',chi(i)
!end do 
!close(35)


end program galaxies

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

