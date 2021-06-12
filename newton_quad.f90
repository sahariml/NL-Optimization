! programme d'optmisation non linéaire sans contraintes du type Newton
! (c) Sahari Med Lamine 08/2003
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Partie declaration
implicit none
double precision tk,zero,bk,tola
integer, parameter::n=3, nrl=150,nal=300
double precision, dimension(n) ::xk,dk,gkm1,b,x
double precision, dimension(n,n) ::a,M,l,q
integer i,casrl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
xk(1)=1
xk(2)=1.4
xk(3)=1.1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tola=1e-8
tk=1.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Debut programme
open(1,FILE='newchol.dat')
al:do i=1,nal
	write (1,'(i3)') i
	write (1,'(1PG15.7E2 )')tk
	write (1,'(A1)',advance='no') '('
	write (1,'(1PG15.7E2 )',advance='no')xk(1)
	write (1,'(A1)',advance='no') ','
	write (1,'(1PG15.7E2 )',advance='no')xk(2)
	write (1,'(A1)') ')'
	write (1,'(1PG15.7E2 )')f(xk)
	write (1,'(1PG15.7E2 )')norm(g(xk))
	print *,'Iteration num =>',i
	print *,'gk=',g(xk)
	print *,'norm(gk)=',norm(g(xk))
	print *,'------------------------------------------------'
	print *,'xk=',xk
	print *,'fk=',f(xk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(norm(g(xk))<tola) exit al
	call choleski(h(xk),l)
	call res_sys_tri(dk,l,-g(xk))
	xk=xnew(tk,xk,dk)
enddo al
contains
!!!!!!!!!!!!!!!!!!!!!!!!!Les sous programmes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!la fonction!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!Une fonction quadratique!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION f(x)
	double precision f
	double precision,dimension(n) :: x,u
	!double precision,dimension(n,n) :: q
	u=matmul(q,x)
	f=0.5*(dot_product(x,u))+dot_product(b,x)
end function
!!!!!!!!!!!!!!!!!!!!!!!!Le gradien!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION g(x)
	double precision,dimension(n) :: x,g
	!double precision,dimension(n,n) :: q
	g=matmul(q,x)+b
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!Le hessien!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION h(x)
	double precision,dimension(n) :: x
	double precision,dimension(n,n) :: h
	h=q
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!xnew
FUNCTION xnew(t,x,d)
	double precision t
	double precision,dimension(n) :: xnew,x,d
	xnew=x+(t*d)
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Norme
FUNCTION norm(x)
	double precision norm
	double precision,dimension(n) :: x
	norm=sqrt(dot_product(x,x))
end function
!!!!!!!!!!!!!!!!!!!!!!!!RL exacte!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rl_exacte(t,x,d)
	double precision t
	double precision,dimension(n)::x,d,u
	u=matmul(h(x),d)
	t=-dot_product(d,g(x))/dot_product(d,u)
end subroutine rl_exacte
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end
