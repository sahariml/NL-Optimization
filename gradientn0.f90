program gradientn0
implicit none
double precision tk,zero,tola
integer, parameter::n=2, nrl=50,nal=20000
double precision, dimension(n) ::xk,dk
integer i,casrl
xk(1)=dble(-1.2)
xk(2)=dble(1.0)
tola=dble(1e-5)
write(*,'(a69)',advance='no'),'Recherche lineaire:1->Armijo, 2->Goldestein, 3->Wolfe, 4->Dichotomie :'
read *,casrl
!!!!!!!!!!!!!!!!!!!!!Creation d'un fichier donnees!!!!!!!!!!!!!!!!!!!
open(1,FILE='steepest.dat')
open(2,FILE='x.dat')
open(3,FILE='y.dat')
print *,'f=,g=',f(xk),g(xk)
write(*,'(1F15.5 )'),xk
read *,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Debut programme
al:do i=1,nal
	if(norm(g(xk))<tola) exit al
	print *,'Iteration num =>',i
	print *,'------------------------------------------------'
	print *,'tk=',tk
	print *,'xk=',xk
	print *,'|gk|=',norm(g(xk))
	print *,'fk=',f(xk)
	write (1,*) i
	write (1,*),tk
	write (1,'(A1)',advance='no') '('
	write (1,'(1F15.5 )',advance='no'),xk(1)
	write (1,'(A1)',advance='no') ','
	write (1,'(1F15.5)',advance='no'),xk(n)
	write (1,'(A1)') ')'
	write (1,'(1F15.5)'),f(xk)
	write (1,'(1E15.5)'),norm(g(xk))
	write (2,'(1F15.5)'),xk(1)
	write (3,'(1F15.5)'),xk(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	dk=-g(xk)
	select case(casrl)
	case(1)
		call rl_armijo(tk,xk,dk)
	case(2)
		call rl_goldestein(tk,xk,dk)
	case(3)
		call rl_wolfe(tk,xk,dk)
	case(4)
		call rl_dichotomie(tk,xk,dk)
	endselect
	xk=xnew(tk,xk,dk)
enddo al
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!Les sous programmes
!!!!!!!!!!!!!!!!!!!!!!!!!!!la fonction test
FUNCTION f(x)
	double precision f
	double precision,dimension(n) :: x
	f=(10*(x(1)**2-x(2))**2)+((x(1)-1)**2)
end function
!!!!!!!!!!!!!!!!!!!!!!!!Le gradient
FUNCTION g(x)
	double precision,dimension(n) :: x,g
	g(1)=(40*x(1)*((x(1)**2)-x(2)))+(2*x(1)-2)
	g(2)=(-20*(x(1)**2))+(20*x(2))
end function
!!!!!!!!!!!!!!!!!theta
FUNCTION theta(t,x,d)
	double precision theta,t
	double precision,dimension(n) :: x,d
	theta=f(x+t*d)
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!theta'
FUNCTION thetap(t,x,d)
	double precision thetap,t
	double precision,dimension(n) :: x,d
	thetap=dot_product(g(x+t*d),d)
end function
!!!!!!!!!!!!!!!!!!!xnew
FUNCTION xnew(t,x,d)
	double precision t
	double precision,dimension(n) :: xnew,x,d
	xnew=x+(t*d)
end function
!!!!!!!!!!!!!!!!!!!!!Norme
FUNCTION norm(x)
	double precision norm
	double precision,dimension(n) :: x
	norm=sqrt(dot_product(x,x))
end function
!!!!!!!!!!!!!!!!!!!!!RL Amijo
subroutine rl_armijo(t,x,d)
	double precision, parameter ::m=.3d0
	double precision t,t1
	double precision,dimension(n)::x,d
	integer i
	t1=0.01
	rl: do i=1,nrl
		t=(t1**i)
		if (theta(t,x,d) <= (m*t*thetap(zero,x,d))+theta(zero,x,d)) exit rl
	enddo rl
10 end subroutine rl_armijo
!!!!!!!!!!!!!!!!!!!!!!!RL Goldestein!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rl_goldestein(t,x,d)
	integer i
	double precision, parameter ::m1=.3d0,m2=.7d0,alf=10d0,lambda=.5d0
	double precision t,td,tg
	double precision,dimension(n)::x,d
	t=0.001
	td=0
	tg=0
	rl:do i=1,nrl
		if ((theta(t,x,d) <= (m1*t*thetap(zero,x,d))+theta(zero,x,d)).and.&
		&(theta(t,x,d) >= (m2*t*thetap(zero,x,d))+theta(zero,x,d))) exit rl
		if(theta(t,x,d) > (m1*t*thetap(zero,x,d))+theta(zero,x,d)) then
			td=t
		end if
		if (theta(t,x,d) < (m2*t*thetap(zero,x,d))+theta(zero,x,d)) then
			tg=t
		end if
		if (td==0) then
			t=alf*tg
		else
			t=td*lambda+((1-lambda)*tg) !t=(tg+td)/2
		end if
	enddo rl
end subroutine rl_goldestein
!!!!!!!!!!!!!!!!!!!RL Wolfe!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rl_wolfe(t,x,d)
	integer i
	double precision, parameter ::m1=.3d0,m2=.7d0,alf=10d0,lambda=.5d0
	double precision t,t1,td,tg
	double precision,dimension(n)::x,d
	t=0.02 !0.2
	td=0
	tg=0
	rl:do i=1,nrl
		if ((theta(t,x,d) <= (m1*t*thetap(zero,x,d))+theta(zero,x,d)).and.(thetap(t,x,d) >= (m2*thetap(zero,x,d)))) exit rl
		if(theta(t,x,d) > (m1*t*thetap(zero,x,d))+theta(zero,x,d)) then
			td=t
		end if
		if ((theta(t,x,d) <= (m1*t*thetap(zero,x,d))+theta(zero,x,d)).and.(thetap(t,x,d)<(m2*thetap(zero,x,d)))) then
			tg=t
		end if
		if (td==0) then
			t=alf*tg
		else
			t=td*lambda+((1-lambda)*tg)
		end if
	enddo rl
end subroutine rl_wolfe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rl_dichotomie(t,x,d)
double precision a,b,lmbda,mu,t,epsilon,l,mid
double precision,dimension(n)::x,d
integer nrl,i
l=1e-8
epsilon=0.001
a=0d0
b=0.02
nrl=100
do i=1,nrl
if ((b-a)<l) exit
  mid=(a+b)/2
  t=mid
  lmbda=mid-epsilon
  mu=lmbda+(2*epsilon)
  if (theta(lmbda,x,d)<theta(mu,x,d)) then
   b=mu
  else
   a=lmbda
  end if
end do
endsubroutine rl_dichotomie
end
