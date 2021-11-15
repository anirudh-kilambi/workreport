real*8 function DET(aa)
real*8 aa(:,:)
real*8 tmp,c(size(aa,dim=1),size(aa,dim=2))
real*8 max
	integer i,j,k,l,m,num(size(aa,dim=1))
	n=size(aa,dim=1)
	det=1.	
	do k=1,n
		max=aa(k,k);num(k)=k;
		do i=k+1,n 
			if(abs(max)<abs(aa(i,k))) then
				max=aa(i,k)
				num(k)=i
			endif
		enddo
		if (num(k)/=k) then
			do l=k,n 
				tmp=aa(k,l)
				aa(k,l)=aa(num(k),l)
				aa(num(k),l)=tmp
			enddo
			det=-1.*det
		endif
		do m=k+1,n
			c(m,k)=aa(m,k)/aa(k,k)
			do l=k,n 
				aa(m,l)=aa(m,l)-c(m,k)*aa(k,l)
			enddo
		enddo !There we made matrix triangular!	
	enddo

	do i=1,n
	det=det*aa(i,i)
	enddo
	return
end function

integer, dimension(10,10) :: mat_a


mat_a = reshape((/10,	11,	12,	13,	14,	15,	16,	17,	18,	19,	20,	21,	22,	23,	24,	25,	26,	27,	28,	29,	30,31,32,33,34,/), shape(mat_a)

output = det(mat_a)