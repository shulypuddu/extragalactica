program suma
implicit none !anula las definiciones previas de las variables
real :: a , b , s !los dos puntos VAN JUNTOS

!print*, 'Suma  ' ! escribe en pantalla 
!write(*,*) 'Hola fortran' !tambien escribe en pantalla ES EL Q VAMOS A USAR MAS SEGUIDO 

print*, 'Ingrese dos numeros: '
read(*,*) a , b
s=a+b
write(*,*) 'Suma= ', s

end program suma
