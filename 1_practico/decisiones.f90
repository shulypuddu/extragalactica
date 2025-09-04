program decisiones
implicit none
real:: x,y
read(*,*) x,y
if (y<0.) print*, 'y es estrictamente negativo'
if (x==0) then
	print*, 'x es nulo'
else if (x>0.) then
	print*, 'x es positivo '
else if (x<0.) then 
	print*, 'x es negativo'
end if 

if ((x>0. .AND. y>0.) .OR. (x<0. .AND. y<0.)) then
	print*, 'x*y es positivo'
else 
	print*, 'x*y es negativo o nulo'
end if

end program decisiones
