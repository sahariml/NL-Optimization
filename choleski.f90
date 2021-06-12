subroutine choleski(a,l)
double precision, dimension(n,n) ::a,l
integer i,j
!!!!!!!!!!!!!!!!!!!!!!!!!!!
l(:,:)=0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,n
	l(i,i)=sqrt(a(i,i)-somme(l(i,:)*l(i,:),i-1))
	do j=i+1,n
		l(j,i)=(a(i,j)-somme(l(i,:)*l(j,:),i-1))/(l(i,i))
	enddo
enddo
end subroutine choleski
