program steepest
implicit none
integer :: k
double precision :: tk
double precision, dimension(:), allocatable :: xk,x1,dk
double precision, parameter ::alpha=dble(100),eps=dble(1e-5)
integer, parameter :: n=2,kmax=6000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=10,file='x.txt',form='unformatted')
open(unit=20,file='y.txt')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(xk(n),x1(n),dk(n))
x1(1)=-1.2
x1(2)=1.0
print *,'norm=',g(x1),norm(x1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print *,'f(',x1(1),',',x1(2),')=',f(x1)
print *,'g(',x1(1),',',x1(2),')=',g(x1)
xk=x1
ALG: do k=1,kmax
	dk=-g(xk)
	if (norm(g(xk))>eps) then
		call dichotomie(xk,dk,tk)
		xk=xnew(xk,dk,tk)
		print *,'xk =',xk
		print *,'t =',tk
		print *,'gk,norm =',g(xk),norm(g(xk))
		write(10),xk(1)
		write(20,*),xk(2)
	else
		exit ALG
	endif
enddo ALG
print *,'alpha =',alpha
print *,'La solution appchee est :',xk
print *,'Le nombre d''iter: ',k
 close(10)
 close(20)
deallocate(xk,x1,dk)
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!
function f(x)
double precision :: f
double precision, dimension(n) :: x
	f=alpha*(x(1)**2-x(2))**2+(x(1)-1.0)**2
end function
function g(x)
!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision, dimension(n) :: x,g
	g(1)=alpha*4.0*(x(1)**2-x(2))*x(1)+2.0*(x(1)-1.0)
	g(2)=2.0*alpha*(x(2)-x(1)**2)
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!
function norm(x)
double precision :: norm
double precision, dimension(n) :: x
	norm=sqrt(x(1)**2+x(2)**2)
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION theta(x,d,t)
	double precision theta,t
	double precision,dimension(n) :: x,d
	theta=f(x+t*d)
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION thetap(x,d,t)
	double precision thetap,t
	double precision,dimension(n) :: x,d
	thetap=dot_product(g(x+t*d),d)
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION xnew(x,d,t)
	double precision t
	double precision,dimension(n) :: xnew,x,d
	xnew=x+(t*d)
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION mid(x,y)
	double precision x,y,mid
	mid=(x+y)/dble(2.0)
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dichotomie(x,d,t)
implicit none
double precision, dimension(n) :: x,d
double precision :: a,b,t
double precision, parameter :: tol=1e-5,delta=dble(10.0)
integer i
a=dble(0.0)
b=1e-3
do i=1,100
	if (thetap(x,d,b)>=dble(0.0)) then
		exit 	
	else	
		b=b*delta
	endif
end do
RL: do i=1,150
	t=mid(a,b)
	if (dabs(thetap(x,d,t))<tol) then
		exit RL
	elseif (thetap(x,d,t)<dble(0.0)) then
		a=t
	else
		b=t
	endif
enddo RL
end subroutine dichotomie
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program steepest
