! Programme d'optmisation non lineaire sans contraintes du type L-BFGS
!avec recherche lineaire de: Goldestein, Armijo, Wolf
! (c) M.L. Sahari 18/08/2003
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!Partie declaration
implicit none
double precision tk,zero,tola,temps_debut,temps_fin
integer, parameter::m=5,nrl=100,nal=3000
double precision,dimension(:),allocatable::xk,xkm1,dk,s1
double precision, dimension(:,:),allocatable ::gamma,delta
integer i,j,casrl,numprob,n
!!!!!!!!!!!!!!!!!!!!Inititialisation
write(*,'(a96)',advance='no')'Probleme num:1->Dixon (PI), 2->Oren (PII), 3->Powel (PIII), 4->Rosenbrock (PIV), 5->Wood (PV):'
read *,numprob
write(*,'(a65)',advance='no')'Dimension ? (PI:n>2, PII:n>0, PIII:n=4*i, PIV:n=2*i, PV:n=4*i):'
read *,n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(xk(n),xkm1(n),dk(n),s1(n))
allocate(delta(n,m),gamma(n,m))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
select case(numprob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(1) !Dixon
do i=1,n
	xk(i)=0.5
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(2)	!Oren
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,n
	xk(i)=1
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(3)	!Powel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,n/4
	xk(2*i-3)=3
	xk(2*i-2)=-1
	xk(2*i-1)=0
	xk(4*i)=1
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(4)	!Rosenbrock
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,n
	xk(2*i)=1
	xk(2*i-1)=-1.2
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(5)	!Wood
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,n
	xk(i)=0
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endselect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tola=1e-5
s1(1:n)=1
write(*,'(a56)',advance='no')'Line Search Methods:1->Armijo, 2->Goldestein, 3->Wolf :'
read *,casrl
!!!!!!!!!!!!!!!!!!!!Create a data files
select case(casrl)
case(1)
	open(1,FILE='l_bfgs_arm.dat')
case(2)
	open(1,FILE='l_bfgs_glst.dat')
case(3)
	open(1,FILE='l_bfgs_wolf.dat')
end select
dk=-g(xk)
!!!!!!!!!!!!!!!!!!!!Debut du programme
CALL CPU_TIME (temps_debut)
do i=1,m
	write (1,'(i3)') i
	write (1,'(1PG15.7E2 )')tk
	write (1,'(1PG15.7E2 )')f(xk)
	write (1,'(1PG15.7E2 )')norm(g(xk))
	print *,'Iteration num =>',i
	print *,'------------------------------------------------'
	print *,'fn=',f(xk)
	print *,'|gk|=',norm(g(xk))
	if(norm(g(xk))<tola) exit
!!!!!!!!!!!!!!!!!!!!Recherche lineaire
	rln:select case(casrl)
	case(1)
		call rl_armijo(tk,xk,dk,n)
	case(2)
		call rl_goldestein(tk,xk,dk,n)
	case(3)
		call rl_wolf(tk,xk,dk,n)
	end select rln
!!!!!!!!!!!!!!!!!!!!calcul du nouveau point
	xkm1=xk
	xk=xnew(tk,xk,dk)
	delta(1:n,i)=xk-xkm1
	gamma(1:n,i)=g(xk)-(xkm1)
	call newdk(i,delta,gamma,s1,g(xk),dk,n)
enddo
do i=m+1,nal
	write (1,'(i3)') i
	write (1,'(1PG15.7E2 )')tk
	write (1,'(1PG15.7E2 )')f(xk)
	write (1,'(1PG15.7E2 )')norm(g(xk))
	print *,'Iteration num =>',i
	print *,'------------------------------------------------'
	print *,'fk=',f(xk)
	print *,'|gk|=',norm(g(xk))
	if(norm(g(xk))<tola) exit
	rl:select case(casrl)
	case(1)
		call rl_armijo(tk,xk,dk,n)
	case(2)
		call rl_goldestein(tk,xk,dk,n)
	case(3)
		call rl_wolf(tk,xk,dk,n)
	end select rl
	xkm1=xk
	xk=xnew(tk,xk,dk)
!!!	s1(1:n)=dot_product(gamma(:,m),delta(:,m))/dot_product(gamma(:,m),gamma(:,m))
	do j=1,m-1
		delta(:,j)=delta(:,j+1)
		gamma(:,j)=gamma(:,j+1)
	enddo
	delta(:,m)=xk-xkm1
	gamma(:,m)=g(xk)-g(xkm1)
	s1(1:n)=dot_product(gamma(:,m),delta(:,m))/dot_product(gamma(:,m),gamma(:,m))
	call newdk(m,delta,gamma,s1,g(xk),dk,n)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL CPU_TIME (temps_fin)
write (1,'(A19)',advance='no') 'le temps machine =>'
write (1,'(1PG15.7E2 )',advance='no') temps_fin - temps_debut,'seconde'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!Subroutine
!
!!!!!!!!!!!!!!!!!!!!functions
!
!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!P1 :Dixon generalisee
!!!										n>2
!!!!!!!!!!!!!!!!!!!!!
FUNCTION f1(x)
	double precision f1
	double precision,dimension(n) :: x
	integer i
	f1=0
	do i=2,n
		f1=f1+(i*((2*x(i)**2)-x(i-1))**2)
	enddo
	f1=f1+((x(1)-1)**2)
end function
!!!!!!!!!!!!!!!!!!!!Le gradien
FUNCTION g1(x)
	integer i
	double precision,dimension(n) :: x,g1
		g1(1)=(2*(x(1)-1))-(4*((2*x(2)**2)-x(1)))
		do i=2,n-1
			g1(i)=(8*i*x(i)*((2*x(i)**2)-x(i-1)))-(2*(i+1)*((2*x(i+1)**2)-x(i)))
		enddo
		g1(n)=(8*n*x(n)*((2*x(n)**2)-x(n-1)))
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!P2: Oren generalis�e
FUNCTION f2(x)
	double precision f2
	double precision,dimension(n) :: x
	integer i
	f2=0
	do i=1,n
		f2=f2+(i*(x(i)**2))
	enddo
	f2=f2**2
end function
!!!!!!!!!!!!!!!!!!!!Le gradien
FUNCTION g2(x)
	integer i
	double precision,dimension(n) :: x,g2
	do i=1,n
		g2(i)=4*(sqrt(f2(x)))*i*x(i)
	enddo
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!P3: Powel generalisee
!!!										n= 4*i, i=1,2,..
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION f3(x)
	double precision f3,e1,e2,e3,e4
	double precision,dimension(n) :: x
	integer i
	f3=0d0
	do i=1,n/4
		e1=x(4*i-3)+(10*x(4*i-2))
		e2=x(4*i-1)-x(4*i)
		e3=x(4*i-2)-(2*x(4*i-1))
		e4=x(4*i-3)-x(4*i)
		f3=f3+(e1**2)+(5*(e2**2))+(e3**4)+(10*(e4**4))
	enddo
end function
!!!!!!!!!!!!!!!!!!!!Le gradien
FUNCTION g3(x)
	double precision e1,e2,e3,e4
	double precision,dimension(n) :: x,g3
	integer i
	do i=1,n/4
		e1=x(4*i-3)+(10*x(4*i-2))
		e2=x(4*i-1)-x(4*i)
		e3=x(4*i-2)-(2*x(4*i-1))
		e4=x(4*i-3)-x(4*i)
		g3(4*i-3)=(2*e1)+(40*(e4**3))
		g3(4*i-2)=(20*e1)+(4*(e3**3))
		g3(4*i-1)=(10*e2)-(8*(e3**3))
		g3(4*i)=(-10*e2)-(40*(e4**3))
	enddo
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!P4: Rosenbrock generalisee
!!!								n= 2*i, i=1,2,..
!!!!!!!!!!!!!!!!!!!!
FUNCTION f4(x)
	double precision f4
	double precision,dimension(n) :: x
	integer i
	f4=0
	do i=1,n/2
		f4=f4+((100*(x(2*i)-x(2*i-1)**2)**2)+((1-x(2*i-1))**2))
	enddo
end function
!!!!!!!!!!!!!!!!!!!!Le gradien
FUNCTION g4(x)
	integer i
	double precision,dimension(n) :: x,g4
	do i=1,n/2
		g4(2*i-1)=(400*x(2*i-1)*((x(2*i-1)**2)-x(2*i)))+((2*x(2*i-1))-2)
		g4(2*i)=200*(x(2*i)-(x(2*i-1)**2))
	enddo
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!P5: Wood generalisee
!									n= 4*i, i=1,2,..
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION f5(x)
	double precision f5,inv,ten,e1,e2,e3,e4,e5,e6
	double precision,dimension(n) :: x
	integer i
	f5=0
	inv=.1d0
	ten = 10.d0
	do i=1,n/4
	e1=x(4*i-2)-(x(4*i-3)**2)
	e2=1-x(4*i-3)
	e3=x(4*i)-(x(4*i-1)**2)
	e4=1-x(4*i-1)
	e5=x(4*i-2)+x(4*i)-2
	e6=x(4*i-2)-x(4*i)
	f5=f5+((100*(e1**2))+(e2**2)+(90*(e3**2))+(e4**2)+(10*(e5**2))+(0.1*(e6**2)))
	enddo
end function
!!!!!!!!!!!!!!!!!!!!Le gradien!!!!!!!!!!!!!!!!!!!!
FUNCTION g5(x)
	integer i
	double precision inv,ten,e1,e2,e3,e4,e5,e6
	double precision,dimension(n) :: x,g5
	inv=.1d0
	ten = 10.d0
	do i=1,n/4
	e1=x(4*i-2)-(x(4*i-3)**2)
	e2=1-x(4*i-3)
	e3=x(4*i)-(x(4*i-1)**2)
	e4=1-x(4*i-1)
	e5=x(4*i-2)+x(4*i)-2
	e6=x(4*i-2)-x(4*i)
	g5(4*i-3)=(-40*ten*x(4*i-3)*e1)-(2*e2)
	g5(4*i-2)=(20*ten*e1)+(2*ten*e5)+(inv*2*e6)
	g5(4*i-1)=(-36*ten*x(4*i-1)*e3)-(2*e4)
	g5(4*i)=(18*ten*e3)+(2*ten*e5)-(inv*2*e6)
	enddo
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION f6(x)
	double precision f6,e1,e2,e3,e4,e5,e6,ten,inv
	double precision,dimension(n) :: x
	inv=.1d0
	ten = 10.d0
	e1=x(2)-(x(1)**2)
	e2=1-x(1)
	e3=x(4)-(x(3)**2)
	e4=1-x(3)
	e5=x(2)+x(4)-2
	e6=x(2)-x(4)
	f6=(100*(e1**2))+(e2**2)+(90*(e3**2))+(e4**2)+(10*(e5**2))+(0.1*(e6**2)) !((x(1)-2)**4)+((x(1)-2*x(2))**2) bazara
end function
!!!!!!!!!!!!!!!!!!!!Le gradien
FUNCTION g6(x)
	double precision e1,e2,e3,e4,e5,e6,ten,inv
	double precision,dimension(n) :: x,g6
	inv=.1d0
	ten = 10.d0
	e1=x(2)-(x(1)**2)
	e2=1-x(1)
	e3=x(4)-(x(3)**2)
	e4=1-x(3)
	e5=x(2)+x(4)-2
	e6=x(2)-x(4)
	g6(1)=(-40*ten*x(1)*e1)-(2*e2) !((x(1)**2)+(2*x(2)**2))!(4*(x(1)-2)**3)+(2*(x(1)-2*x(2)))
	g6(2)=(20*ten*e1)+(2*ten*e5)+(inv*2*e6) !((x(1)**2)+(2*x(2)**2))! -4*(x(1)-2*x(2))
	g6(3)=(-36*ten*x(3)*e3)-(2*e4)
	g6(4)=(18*ten*e3)+(2*ten*e5)-(inv*2*e6)
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								Affection des problemes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION f(x)
	double precision f
	double precision,dimension(n) :: x
	select case(numprob)
	case(1) !Dixon
		f=f1(x)
	case(2) !Oren
		f=f2(x)
	case(3) !Powel
		f=f3(x)
	case(4) !Rosenbrock
		f=f4(x)
	case(5) !Wood
		f=f5(x)
	case(6) !?
		f=f6(x)
	end select
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION g(x)
	double precision,dimension(n) :: x,g
	select case(numprob)
	case(1) !Dixon
		g=g1(x)
	case(2) !Oren
		g=g2(x)
	case(3) !Powel
		g=g3(x)
	case(4) !Rosenbrock
		g=g4(x)
	case(5) !Wood
		g=g5(x)
	case(6) !?
		g=g6(x)
	end select
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!theta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION theta(t,x,d)
	double precision theta,t
	double precision,dimension(n) :: x,d
	theta=f(x+t*d)
end function
!!!!!!!!!!!!!!!!!!!!theta'
FUNCTION thetap(t,x,d)
	double precision thetap,t
	double precision,dimension(n) :: x,d
	thetap=dot_product(g(x+t*d),d)
end function
!!!!!!!!!!!!!!!!!!!!xnew
FUNCTION xnew(t,x,d)
	double precision t
	double precision,dimension(n) :: xnew,x,d
	xnew=x+(t*d)
end function
!!!!!!!!!!!!!!!!!!!!Norme
FUNCTION norm(x)
	double precision norm
	double precision,dimension(n) :: x
	norm=sqrt(dot_product(x,x))
end function
!!!!!!!!!!!!!!!!!!!!RL Amijo
subroutine rl_armijo(t,x,d,n)
	integer i,n
	double precision, parameter ::m=.3d0
	double precision t,t1
	double precision,dimension(n)::x,d
	t1=0.7       !le pas alfa1
	rl: do i=1,nrl
		t=(t1**(i-1))
		if (theta(t,x,d) <= (m*t*thetap(zero,x,d))+theta(zero,x,d)) exit rl
	enddo rl
10 end subroutine rl_armijo
!!!!!!!!!!!!!!!!!!!!RL Goldestein
subroutine rl_goldestein(t,x,d,n)
	integer i,n
	double precision, parameter ::m1=.3d0,m2=.7d0,alf=10d0,lambda=.7d0
	double precision t,t1,td,tg
	double precision,dimension(n)::x,d
	t=1
	td=0
	tg=0
	rl:do i=1,nrl
		if ((theta(t,x,d) <= (m1*t*thetap(zero,x,d))+theta(zero,x,d)).and.(theta(t,x,d)&
		 >= (m2*t*thetap(zero,x,d))+theta(zero,x,d))) exit rl
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
!!!!!!!!!!!!!!!!!!!!RL Wolf
subroutine rl_wolf(t,x,d,n)
	integer i,n
	double precision, parameter ::m1=.3d0,m2=.7d0,alf=10d0,lambda=.7d0
	double precision t,t1,td,tg
	double precision,dimension(n)::x,d
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
!!!!!!!!!!!!!!!!!!!!Update L-BFGS
subroutine newdk(k,delta,gamma,s1,g,v,n)
	integer i,k,n
	double precision sigmma,c
	double precision,dimension(n,m)::delta,gamma
	double precision,dimension(n)::s1,g,v
	double precision,dimension(m)::sigma,ro
	v=-g
	do i=k,1,-1
		ro(i)=1/(dot_product(gamma(:,i),delta(:,i)))
		sigma(i)=ro(i)*(dot_product(delta(:,i),v))
		v=v-(sigma(i)*gamma(:,i))
	enddo
	v(:)=v(:)*s1(:)
	do i=1,k
		sigmma=ro(i)*(dot_product(gamma(:,i),v))
		v=v+((sigma(i)-sigmma)*delta(:,i))
	enddo
end subroutine newdk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end
