PROGRAM test
USE determinant

real, dimension(2,2) :: input 
real :: output

input = reshape((/ 12.254, 43.872, 98.399, 0.010 /), shape(input))



call find_determinant(input, output)

print *, output


END PROGRAM test