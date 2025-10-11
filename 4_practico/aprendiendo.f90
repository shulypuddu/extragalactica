program practico
real, allocatable, dimension(:) :: x
double precision :: a !doble precision ocupa 64 bits en vez de 32
real(8) :: b !otra forma de escribir el de doble precision
real :: k1 
real :: k2 = selected_real_kind(14)


a=1.20805978
b=1.20805978
k1=1.20805978
k2=1.20805978

print*, a , b, k1, k2
end program practico