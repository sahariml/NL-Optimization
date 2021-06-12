implicit none
double precision :: tdebut,tfin
real ::d,f
integer :: i,j,n,m
real, dimension(:,:),allocatable :: A,B,L
open(10,file='matrice.txt')
open(20,file='matrice_out.txt')
n=3
allocate(A(1:n,1:n),L(1:n,1:n))
do i=1,n
read (10,*)A(i,:)
enddo
print *,"la matrice A:"
write (*,fmt='(3F10.5)'),A
!write(*,*)A
!write(20,*)A
call cholesky(A,L,n)
print *,"la matrice L:"
write (*,fmt='(3F10.5)'),L
write (20,fmt='(3F10.5)'),L
!!!!!!!!!!!!!!!!!!!!!!!!
!do i=1,n
!	L(i,i)=sqrt(A(i,i)-somme(L(i,:)*L(i,:),i-1))
!	do j=i+1,n
!		L(j,i)=(A(i,j)-somme(L(i,:)*L(j,:),i-1))/(L(i,i))
!	enddo
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write (*,fmt='(4F8.5)'),L
!write (20,*),"LLLLLLLLLLLLLLLLLLLLL"
!write (20,*),L
 close(20)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call cpu_time(tdebut)
contains
subroutine cholesky(B,L,m)
real, dimension(n,n):: B,L
integer :: i,j,m
!!!!!!!!!!!!!!!!!!!!!!!!
!write (*,fmt='(4F8.5)'),B
do i=1,m
	L(i,i)=sqrt(B(i,i)-somme(L(i,:)*L(i,:),i-1))
	do j=i+1,m
		L(j,i)=(B(i,j)-somme(L(i,:)*L(j,:),i-1))/(L(i,i))
	enddo
enddo
end subroutine cholesky
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function somme(V,m)
real :: somme
real, dimension(n) :: V
integer i,m
somme=0.0
if (m<1) then
	somme=0.0
else
	do i=1,m
!print *,"vi",V(i),somme
		somme=somme+V(i)
	enddo
endif
end function
end
