! programme d'optmisation non linéaire sans contraintes du type Newton-Choleski
!avec RL de Armijo, Goldestein, Wolf
! (c) Sahari Med Lamine 08/2003
!----------------------------------------------------------------------------------------------
!---------------------------------Partie declaration-------------------------------------------
implicit none
double precision tk,zero,bk,tola
integer, parameter::n=3, nrl=150,nal=300
double precision, dimension(n) ::xk,dk,gkm1,b,x
double precision, dimension(n,n) ::a,M
integer i,casrl
!----------------------------------------------------------------------------------------------
xk(1)=1
xk(2)=1.4
xk(3)=1.1
tola=1e-8
write(*,'(a60)',advance='no'),'Recherche lineaire:1->Armijo, 2->Goldestein, 3->Wolfe, 4->RLE :'
read *,casrl
!print *,"Q=",h1(xk)
!call choleski(xk,h1(xk),xk)

!---------------------------------Debut programme----------------------------------------------
open(1,FILE='newchol.dat')
open(2,FILE='x.dat')
open(3,FILE='y.dat')
al:do i=1,nal
	write (1,'(i3)') i
	write (1,'(1PG15.7E2 )'),tk
	write (1,'(A1)',advance='no') '('
	write (1,'(1PG15.7E2 )',advance='no'),xk(1)
	write (1,'(A1)',advance='no') ','
	write (1,'(1PG15.7E2 )',advance='no'),xk(2)
	write (1,'(A1)') ')'
	write (1,'(1PG15.7E2 )'),f(xk)
	write (1,'(1PG15.7E2 )'),norm(g(xk))
	write (2,'(1PG15.7E2 )'),xk(1)
	write (3,'(1PG15.7E2 )'),xk(2)
	print *,'Iteration num =>',i
	print *,'gn=',g(xk)
	print *,'------------------------------------------------'
	print *,'xn=',xk
	print *,'fn=',f(xk)
	print *,'norm=',norm(g1(xk))
!----------------------------------------------------------------------------------------------
	if(norm(g1(xk))<tola) exit al
	a=h1(xk)
	b=-g1(xk)
	call choleski(dk,a,b)
	rln:select case(casrl)
	case(1)
		call rl_armijo(tk,xk,dk)
	case(2)
		call rl_goldestein(tk,xk,dk)
	case(3)
		call rl_wolf(tk,xk,dk)
	case(4)
		call rl_exacte(tk,xk,dk)
	end select rln
		call rl_armijo(tk,xk,dk)
!print *,'tk=',tk
!read *,

	xk=xnew(tk,xk,dk)
enddo al
contains
!!!!!!!!!!!!!!!!!!!!!!!!!Les sous programmes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!la fonction!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION f(x)
	double precision f
	double precision,dimension(n) :: x
	f=(100*(x(1)**2-x(2))**2)+((x(1)-1)**2)
end function
!-------------------------------------Le gradien-----------------------------------------------
FUNCTION g(x)
	double precision,dimension(n) :: x,g
	g(1)=(400*x(1)*((x(1)**2)-x(2)))+(2*x(1)-2)
	g(2)=-200*x(1)**2+(200*x(2))
end function
!-------------------------------------Le hessien-----------------------------------------------
FUNCTION h(x)
	double precision,dimension(n) :: x
	double precision,dimension(n,n) :: h
	h(1,1)=(1200*x(1)**2)-(400*x(2))+2
	h(2,1)=-400*x(1)
	h(1,2)=-400*x(1)
	h(2,2)=200
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!Une fonction quadratique!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION f1(x)
	double precision f1
	double precision,dimension(n) :: x,u,b
	double precision,dimension(n,n) :: q
	q(1,1)=6.69
	q(2,1)=-3.842
	q(3,1)=-1.536
	q(1,2)=-3.842
	q(2,2)=5.54
	q(3,2)=-1.783
	q(1,3)=-1.536
	q(2,3)=-1.783
	q(3,3)=9.287
	b(1)=0
	b(2)=1
	b(3)=0
	u=matmul(q,x)
	f1=0.5*(dot_product(x,u))+dot_product(b,x)
end function
!!!!!!!!!!!!!!!!!!!!!!!!Le gradien!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION g1(x)
	double precision,dimension(n) :: x,g1,b
	double precision,dimension(n,n) :: q
	q(1,1)=6.69
	q(2,1)=-3.842
	q(3,1)=-1.536
	q(1,2)=-3.842
	q(2,2)=5.54
	q(3,2)=-1.783
	q(1,3)=-1.536
	q(2,3)=-1.783
	q(3,3)=9.287
	b(1)=0
	b(2)=1
	b(3)=0
	g1=matmul(q,x)+b
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!Le hessien!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION h1(x)
	double precision,dimension(n) :: x
	double precision,dimension(n,n) :: h1
	h1(1,1)=6.69
	h1(2,1)=-3.842
	h1(3,1)=-1.536
	h1(1,2)=-3.842
	h1(2,2)=5.54
	h1(3,2)=-1.783
	h1(1,3)=-1.536
	h1(2,3)=-1.783
	h1(3,3)=9.287
	!h1(1,1)=6.69
	!h1(2,1)=-3.842
	!h1(3,1)=-1.536
	!h1(1,2)=-3.842
	!h1(2,2)=5.54
	!h1(3,2)=-1.783
	!h1(1,3)=-1.536
	!h1(2,3)=-1.783
	!h1(3,3)=9.287
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------theta------------------------------------------------
FUNCTION theta(t,x,d)
	double precision theta,t
	double precision,dimension(n) :: x,d
	theta=f1(x+t*d)
end function
!-----------------------------------------theta'-----------------------------------------------
FUNCTION thetap(t,x,d)
	double precision thetap,t
	double precision,dimension(n) :: x,d
	thetap=dot_product(g1(x+t*d),d)
end function
!-------------------------------------------xnew-----------------------------------------------
FUNCTION xnew(t,x,d)
	double precision t
	double precision,dimension(n) :: xnew,x,d
	xnew=x+(t*d)
end function
!--------------------------------------------Norme----------------------------------------------
FUNCTION norm(x)
	double precision norm
	double precision,dimension(n) :: x
	norm=sqrt(dot_product(x,x))
end function
!----------------------------------RL Amijo----------------------------------------------------
subroutine rl_armijo(t,x,d)
	integer i
	double precision, parameter ::m=.3d0
	double precision t,t1
	double precision,dimension(n)::x,d
	t1=0.7       !le pas alfa1
	rl: do i=1,nrl
		t=(t1**(i-1))
		if (theta(t,x,d) <= (m*t*thetap(zero,x,d))+theta(zero,x,d)) exit rl
	enddo rl
10 end subroutine rl_armijo
!----------------------------------RL Goldestein----------------------------------------------------
subroutine rl_goldestein(t,x,d)
	double precision, parameter ::m1=.3d0,m2=.7d0,alf=10d0,lambda=.7d0
	double precision t,t1,td,tg
	double precision,dimension(n)::x,d
	integer i
	t=1
	td=0
	tg=0
	rl:do i=1,nrl
		if ((theta(t,x,d) <= (m1*t*thetap(zero,x,d))+theta(zero,x,d)).and.(theta(t,x,d) >= (m2*t*thetap(zero,x,d))&
		&+theta(zero,x,d))) exit rl
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
!!!!!!!!!!!!!!!!!!!!!!!!RL Wolf!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rl_wolf(t,x,d)
	double precision, parameter ::m1=.3d0,m2=.7d0,alf=10d0,lambda=.7d0
	double precision t,t1,td,tg
	double precision,dimension(n)::x,d
	integer i
	t=1
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
			t=td*lambda+((1-lambda)*tg) !t=(tg+td)/2
		end if
	enddo rl
end subroutine rl_wolf
!!!!!!!!!!!!!!!!!!!!!!!!RL exacte!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rl_exacte(t,x,d)
	double precision t
	double precision,dimension(n)::x,d,u
	u=matmul(h1(x),d)
	t=-dot_product(d,g1(x))/dot_product(d,u)
end subroutine rl_exacte
!------------------------------Choleski-----------------------------------------
subroutine choleski(x,a,b)
double precision som
double precision, dimension(n) ::b,x
double precision, dimension(n,n) ::a,l
integer i,j,k
!-------------------------Debut programme--------------------------------------------
!										Partie A
!---------------------Factorisation de Choleski--------------------------------------
l(:,:)=0.0
!!!!!!!!!!!!!!!!!!!!!!
!print *,"la matrice L:"
!write (*,fmt='(3F10.5)'),l
do i=1,n
	l(i,i)=sqrt(a(i,i)-somme(l(i,:)*l(i,:),i-1))
	do j=i+1,n
		l(j,i)=(a(i,j)-somme(l(i,:)*l(j,:),i-1))/(l(i,i))
	enddo
enddo
!print *,"la matrice triangulaire ==>",l
!print *,"la matrice LLt ==>",matmul(l,transpose(l))
!print *,"la matrice Q:"
!write (*,fmt='(3F10.5)'),a
!print *,"la matrice L:"
!write (*,fmt='(3F10.5)'),l
!read *,
!-----------------------------------------------------------------------------------
!											Partie B
!----------------------Resolution d'un systeme triangulaire inf----------------------
x(1)=b(1)/l(1,1)
do i=2,n
	som=0
	do j=1,i-1
		som=som+(l(i,j)*x(j))
	enddo
	x(i)=(b(i)-som)/l(i,i)
enddo
!print *,'x=',x
!-----------------------Resolution d'un systeme triangulaire sup----------------------
b=x
l = TRANSPOSE(l)
x(n)=b(n)/l(n,n)
do i=n-1,1,-1
	som=0
	do j=i+1,n
		som=som+(l(i,j)*x(j))
	enddo
	x(i)=(b(i)-som)/l(i,i)
enddo
!print *,'x===',x
!read *,
end subroutine choleski
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function somme(V,m)
double precision :: somme
double precision, dimension(n) :: V
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
