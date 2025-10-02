! filepath: /mnt/sda2/extragalactica/2_practico/pr2.f90
module lista
implicit none
integer, parameter :: f=6000, g=20000, pi=3.1415 ! f para cosmos.dat, g para galaxies.dat
!si cambiamos la long de la tabla tenemos q cambiar f a mano, no hay forma de automatizarlo
real, dimension(f) :: redshift, coord_comovil, dist_luminosa, dist_angular
real::  z, chi, dL
end module lista

program galaxies
use lista
implicit none
integer :: i, j,  k
character(len=21):: ID 
!guardamos los ID de las galaxias como texto xq tienen una longitud de caracteres mayor a la de un integer
real :: coord

real :: fracDeV_r, velDisp, ra, dec, mu_sup
real, dimension(5) :: petroMag, extinction, modelMag,petroMag_s_ext,modelMag_s_ext , corr_ab=(/-0.036,0.012,0.010,0.028,0.040/)
real :: petroRad_u, petroRad_g, petroRad_r, petroRad_i, petroRad_z
real :: petroR50_r, petroR90_r, c_par, color_u_r, color_g_r
real, dimension(5) :: abs_petroMag, abs_modelMag

external coord

open (unit=25,file='cosmos.dat',status='unknown')
read(25,*)  ! Salta la primera línea (header)
do i = 2 , f
read(25,*) redshift(i), coord_comovil(i), dist_luminosa(i), dist_angular(i) 
end do 
close(25)

open (unit=30,file='tabla.dat',status='old')
read(30,*)  ! Salta la primera línea (header)
write(35,*) 'ra',',','dec',',','z',',','c_par',',','fracDeV_r',',','velDisp',',', &
  'abs_model_mag_u',',','abs_model_mag_g',',','abs_model_mag_r',',','abs_model_mag_i',',','abs_model_mag_z',',', &
  'abs_petro_mag_u',',','abs_petro_mag_g',',','abs_petro_mag_r',',','abs_petro_mag_i',',','abs_petro_mag_z',',', &
  'color_u_r',',','color_g_r',',','mu_sup',',','petroR50_r',',','petroR90_r'

do i = 1 , g
  read(30,*) ID,ra,dec,z,fracDeV_r,velDisp, (petroMag(k),k=1,5), (modelMag(k),k=1,5), &
    (extinction(k),k=1,5),petroRad_u,petroRad_g,petroRad_r,petroRad_i,petroRad_z, &
    petroR50_r,petroR90_r
    
  do k =1,5
    modelMag_s_ext(k)=modelMag(k)-extinction(k)  !mag corregidas por extincion
    petroMag_s_ext(k)=petroMag(k)-extinction(k)
  end do 
    
  if (petroRad_r<1.5) cycle

  if (petroMag_s_ext(3)<14.5 .OR. petroMag_s_ext(3)>17.77) cycle
  

  chi=coord(z)
  dL= (1.+z)*chi

  do k=1,5
    abs_modelMag(k)=modelMag_s_ext(k)-5*log10(dL)-25+corr_ab(k)
    abs_petroMag(k)=petroMag_s_ext(k)-5*log10(dL)-25+corr_ab(k)

  end do

  c_par=petroR90_r/petroR50_r
  color_g_r=abs_modelMag(2)-abs_modelMag(3)
  color_u_r=abs_modelMag(1)-abs_modelMag(3)
  mu_sup= petroMag(3) + 2.5*log10(2*pi*(petroR50_r)**2) + corr_ab(3)
  petroR50_r= petroR50_r* ( (pi/180.)/3600.) * dA*1000  !de arcsec a radianes a kpc 
  petroR90_r= petroR90_r* ( (pi/180.)/3600.) * dA*1000  

  write(35,*) ra,',',dec,',',z,',',c_par,',',fracDeV_r,',', velDisp ,',',&
  (abs_modelMag(k),',',k=1,5),(abs_petroMag(k),',',k=1,5),color_u_r,',',color_g_r,',',&
  mu_sup,',',petroR50_r,',',petroR90_r 
  !me lo escribe en otro archivo llamado for.35 o for.(nro unidad), para leerlo cambiar a extension csv 
  !antes de usarlo quitar todos los espacios.

end do

close(30)


end program galaxies

!_________________________FUNCION INTERPOLACION___________________________________________
function coord(x)
use lista
implicit none
real::x, coord
integer :: i, j
call locate(redshift,f,x,j)
coord= ((coord_comovil(j+1)-coord_comovil(j))/(redshift(j+1)-redshift(j)))*(x-redshift(j)) + coord_comovil(j)

end function coord

include 'locate.f'
