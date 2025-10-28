
!::::::::: PROGRAMA ::::::::::::::::::::::::::::::::::::::::::::::::::::
program leer
character,parameter :: nom1='sample1.dat'
character,parameter :: nom2='sample2.dat'
character,parameter :: nom3='sample3.dat'
integer,parameter:: n1=5712 ,n2= 508330,n3= 27005

real::cls,ra,dec,z,stellarmass
real::sfr,ssfr,Dn4000,class,Zc,sigma,Mvir,Rvir,R200,w,dist,Pe,Ps    
real::Mr01,ur,mur,kr50,C,OH      
real::ig,type,clase,smass,u,g,r,Mcoldgas,Mhotgas,Tau,bt      
integer:: i,j,k 
character:: a

open(unit=15,file='sample1.dat',status='old')
open(unit=115,file='muestra1.dat',status='unknown')
read(15,*)a
write(115,*)a
do i=1,n1
read(15,*) cls,ra,dec,z,Mr01,ur,mur,kr50,C,stellarmass,sfr,ssfr,Dn4000,class,Zc,&
    & sigma,Mvir,Rvir,R200,w,dist,Pe,Ps    

write(115,*) cls,ra,dec,z,Mr01,ur,mur,kr50,C,stellarmass,sfr,ssfr,Dn4000,class,Zc, &
    & sigma,Mvir,Rvir,R200,w,dist,Pe,Ps
end do 

close(115)
close(15)


open(unit=25,file='sample2.dat',status='old')
open(unit=125,file='muestra2.dat',status='unknown')
read(25,*)a
write(125,*)a
do i=1,n3
read(25,*) ra,dec,z,Mr01,ur,mur,kr50,C, stellarmass,sfr,ssfr,Dn4000,OH,w,Pe,Ps
write(125,*)ra,dec,z,Mr01ur,mur,kr50,C, stellarmass,sfr,ssfr,Dn4000,OH,w,Pe,Ps
end do 


close(125)
close(25)



open(unit=35,file='sample3.dat',status='old')
open(unit=135,file='muestra3.dat',status='unknown')
read(35,*)a
write(135,*)a
do i=1,n3
read(35,*) cls,ig,type,clase,smass,u,g,r,sfr,Mcoldgas,Mhotgas,Tau,OH,bt,ssfr    
write(135,*) cls,ig,type,clase,smass,u,g,r,sfr,Mcoldgas,Mhotgas,Tau,OH,bt,ssfr    

end do 

close(135)
close(35)


end program leer
!::::::::: FUNCIONES ::::::::::::::::::::::::::::::::::::::::::::::::::::

