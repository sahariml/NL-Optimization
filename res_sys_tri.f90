subroutine res_sys_tri(x,a,b)
	double precision som
	double precision, dimension(n) ::b,x
	double precision, dimension(n,n) ::a
	integer i,j,k
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	x(1)=b(1)/a(1,1)
	do i=2,n
		som=0
		do j=1,i-1
			som=som+(a(i,j)*x(j))
		enddo
		x(i)=(b(i)-som)/a(i,i)
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	b=x
	a = TRANSPOSE(a)
	x(n)=b(n)/a(n,n)
	do i=n-1,1,-1
		som=0
		do j=i+1,n
			som=som+(a(i,j)*x(j))
		enddo
		x(i)=(b(i)-som)/a(i,i)
	enddo
end subroutine res_sys_tri
