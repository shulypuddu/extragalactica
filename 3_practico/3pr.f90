!::MODULOS:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
module lista 
implicit none 
integer, parameter :: f=6000 , g = 685035 , b = 20 !f para cosmos, g para datos, b cant de bines
real, parameter :: pi=3.1415
real, dimension(f):: redshift, d_luminosa
real :: z, chi, dL
end module lista 

!::PROGRAMA:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

program luminosidad
use lista
implicit none
integer i, j, k, n, m
real :: petro_r, red_r, r50, corr_ab=0.010
real :: petro_s_ext, petro_abs, dm , l_phi
real :: distL, red, vol, dLmax
real :: vmax, vmin, w, zmin, zmax, m_medio
real,dimension(b) ::  phi
real, dimension(5) :: rk_p, rks_p
real, parameter :: rmin=-23. , rmax=-16.

external distL, red, vol
call read_cosmo
zmin=0.
n=0
dm = (rmax-rmin)/real(b) 
phi=0.

open(unit=30, file='datos.dat',status='old')
write(35,*)'z,dL,petro_abs,rks_p,dLmax,zmax,vmax'
do i=1,g
    read(30,*)z,petro_r,red_r,r50,(rk_p(k),k=1,5),(rks_p(k),k=1,5)
    if (z>0.15 .or. z<=0.0) cycle
    if (r50<1.5) cycle
    petro_s_ext = petro_r -red_r
    if (petro_s_ext<14.5 .OR. petro_s_ext>17.77) cycle
    dL=distL(z)
    
    petro_abs=petro_s_ext - 25. - 5.*log10(dL)+corr_ab + rks_p(3)
    if(petro_abs >= rmax.or.petro_abs<=rmin) cycle
    dLmax=10**(-0.2*(petro_abs-rmax+25.))
    zmax=red(dLmax)
    n=n+1
    call qromb(vol,zmin,zmax,vmax) 
    w=1./vmax
    m= int((petro_abs - rmin)/dm)+1
    phi(m)= phi(m)+1.*w
write(35,*)z,',',dL,',',petro_abs,',',rks_p(3),',',dLmax,',',zmax,',',vmax
end do 
close(30)

open(unit=50, file='fun_lum.csv',status='unknown')
write(50,*)'m_medio,phi_hist'

do j = 1,b
    phi(j)=phi(j)*dm
    m_medio = rmin+dm*real(2*j-1)*0.5   
    l_phi= log10(phi(j))
    !print*, rmin, m_medio, phi(j)
    write(50,*) m_medio ,',', l_phi
end do
close(50)


end program luminosidad

!::FUNCIONES:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!--------------------- FUNCIÓN INTERPOLACIÓN dL(z) -----------------------------
function distL(x)
    use lista
    implicit none
    real :: x, distL
    integer :: j 
    call locate(redshift,f,x,j)

    distL = ((d_luminosa(j+1)-d_luminosa(j))/(redshift(j+1)-redshift(j)))*&
            (x-redshift(j)) + d_luminosa(j)
    return 
end function distL

!--------------------- FUNCIÓN INTERPOLACIÓN z(dL) -----------------------------
function red(x)
    use lista
    implicit none
    real:: x, red
    integer :: j
    call locate(d_luminosa,f,x,j)

    red = ((redshift(j+1)-redshift(j))/(d_luminosa(j+1)-d_luminosa(j)))* &
            (x-d_luminosa(j)) + redshift(j)
    return 
end function red

!--------------------- FUNCIÓN VOLUMEN V(z) -----------------------------
function vol(x)
    real :: x, vol, om, ol
    real :: arg1, arg2, distL
    external distL
    ol=0.7
    om=0.3
    arg1=distL(x)**2*(1.+x)**(-2)
    arg2=sqrt(om*(1.+x)**3+ol)
    vol=arg1/arg2
    return
end function vol

!::SUBRUTINAS:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include 'locate.f'
include 'qromb.f'
include 'trapzd.f'
include 'polint.f'
!-------------------------------- LEER COSMOLOGIA --------------------------------
subroutine read_cosmo
        use lista
        implicit none
        integer :: i
        real :: xx
        open(10,file='cosmos.dat',status='old')
        read(10,*)
        do i=1,f
          read(10,*)redshift(i),xx,d_luminosa(i),xx
        end do
        close(10)
end subroutine read_cosmo 