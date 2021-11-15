PROGRAM determinant
IMPLICIT NONE


real :: start, finish
real, dimension(2,2) :: a
real :: b





a = reshape((/ 12.254, 43.872, 98.399, 0.010 /), shape(a))

b = a(1,1) * a(2,2) - a(1,2) * a(2,1)


    
call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish
print *, b






END PROGRAM determinant