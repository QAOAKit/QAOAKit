module BFGS

implicit none

contains

SUBROUTINE dfpmin(p,n,gtol,iter,fret,func,dfunc,graph_number)
!from Numerical Recipes in FORTRAN 77: The art of Scientific Programming
!https://www.goodreads.com/book/show/27815.Numerical_Recipes_in_FORTRAN_77
implicit none

INTEGER iter,n,NMAX,ITMAX,graph_number
Double precision fret,gtol,p(n),func,EPS,STPMX,TOLX
PARAMETER (NMAX=50,ITMAX=200,STPMX=6.d0,EPS=3.d-10,TOLX=4.d0*EPS) 
EXTERNAL dfunc,func
! USES dfunc,func,lnsrch

!assumes all values are of order unity

!Given a starting point p(1:n) that is a vector of length n, 
!the Broyden-Fletcher-Goldfarb-Shanno variant of Davidon-Fletcher-Powell 
!minimization is performed on a function func, using its gradient as calculated 
!by a routine dfunc. The convergence requirement on zeroing the gradient is 
!input as gtol. Returned quantities are p(1:n) (the location of the mini- mum), 
!iter (the number of iterations that were performed), 
!and fret (the minimum value of the function). The routine lnsrch is called to 
!perform approximate line minimizations. 
!Parameters: NMAX is the maximum anticipated value of n; 
!ITMAX is the maximum allowed number of iterations; 
!STPMX is the scaled maximum step length allowed in line searches; 
!TOLX is the convergence criterion on x values.

INTEGER i,its,j
LOGICAL check
Double precision den,fac,fad,fae,fp,stpmax,sum,sumdg,sumxi,temp,test, &
& dg(NMAX),g(NMAX),hdg(NMAX),hessin(NMAX,NMAX), pnew(NMAX),xi(NMAX)
!Calculate starting function value and gradient,

fp=func(n,p,graph_number)
call dfunc(n,p,g,func,graph_number)
sum=0.d0

!and initialize the inverse Hessian to the unit matrix.

do i=1,n
	do j=1,n
		hessin(i,j)=0.d0 
	enddo 
	hessin(i,i)=1.d0
	!Initial line direction.
	xi(i)=-g(i)
	sum=sum+p(i)**2
enddo 

stpmax=STPMX*max(sqrt(sum),float(n))

!Main loop over the iterations.
do its=1,ITMAX  
	iter=its
	!The new function evaluation occurs in lnsrch; save the function value in fp for the next
	!line search. It is usually safe to ignore the value of check. 
	call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check,func,graph_number)
	fp=fret
	do i=1,n 
		!Update the line direction, and the current point.
		xi(i)=pnew(i)-p(i) 
		p(i)=pnew(i)
	enddo 
	test=0.d0
	!Test for convergence on ∆x.
	do i=1,n
		temp=abs(xi(i))/max(abs(p(i)),1.d0)
		if(temp.gt.test)test=temp 
	enddo
	if(test.lt.TOLX)return 
	!Save the old gradient,
	do i=1,n
		dg(i)=g(i) 
	enddo
	!and get the new gradient.
	call dfunc(n,p,g,func,graph_number)
	!Test for convergence on zero gradient.
	test=0.d0
	den=max(fret,1.d0) 
	do i=1,n
		temp=abs(g(i))*max(abs(p(i)),1.d0)/den
		if(temp.gt.test)test=temp 
	enddo
	if(test.lt.gtol)return 
	!Compute difference of gradients,
	do i=1,n
		dg(i)=g(i)-dg(i) 
	enddo
	!and difference times current matrix.
	do i=1,n 
		hdg(i)=0.d0
		do j=1,n 
			hdg(i)=hdg(i)+hessin(i,j)*dg(j)
		enddo 
	enddo
	!Calculate dot products for the denominators.
	fac=0.d0
	fae=0.d0
	sumdg=0.d0
	sumxi=0.d0
 	do i=1,n 
 		fac=fac+dg(i)*xi(i) 
 		fae=fae+dg(i)*hdg(i) 
 		sumdg=sumdg+dg(i)**2 
 		sumxi=sumxi+xi(i)**2
	enddo 
	!Skip update if fac not sufficiently positive.
	if(fac.gt.sqrt(EPS*sumdg*sumxi))then
		fac=1.d0/fac
		fad=1.d0/fae
		!The vector that makes BFGS different from DFP:
		do i=1,n 
			dg(i)=fac*xi(i)-fad*hdg(i) 
		enddo
		!The BFGS updating formula: 
		do i=1,n 
			do j=i,n
				hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j) -fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j)
				hessin(j,i)=hessin(i,j) 
			enddo
		enddo 
	endif

	!Now calculate the next direction to go,
	do i=1,n  
		xi(i)=0.d0
		do j=1,n 
			xi(i)=xi(i)-hessin(i,j)*g(j)
		enddo 
	enddo
!and go back for another iteration. 
enddo 
print *, "too many iterations in dfpmin"
return
END subroutine dfpmin


SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func,graph_number) 
!from Numerical Recipes in FORTRAN 77: The art of Scientific Programming
!https://www.goodreads.com/book/show/27815.Numerical_Recipes_in_FORTRAN_77
implicit none
INTEGER n, graph_number
LOGICAL check
DOUBLE PRECISION f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX 
PARAMETER (ALF=1.d-4,TOLX=1.d-7)
EXTERNAL func
! USESfunc
!Given an n-dimensional point xold(1:n), the value of the function and gradient there, 
!fold and g(1:n), and a direction p(1:n), finds a new point x(1:n) along the direction p 
!from xold where the function func has decreased “sufficiently.” The new function value is
!returned in f. stpmax is an input quantity that limits the length of the steps so that 
!you do not try to evaluate the function in regions where it is undefined or subject to 
!overflow. p is usually the Newton direction. The output quantity check is false on a normal 
!exit. It is true when x is too close to xold. In a minimization algorithm, this usually 
!signals convergence and can be ignored. However, in a zero-finding algorithm the calling 
!program should check whether the convergence is spurious.

!Parameters: ALF ensures sufficient decrease in function value; TOLX is the convergence
!criterion on ∆x.
INTEGER i
Double Precision a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,test,tmplam
check=.false. 
sum=0.d0

do i=1,n
	sum=sum+p(i)*p(i) 
enddo
sum=sqrt(sum) 
!scale if attempted step is too big
if(sum.gt.stpmax) then
	do i=1,n 
		p(i)=p(i)*stpmax/sum
	enddo
endif
slope=0.d0 
do i=1,n
	slope=slope+g(i)*p(i) 
enddo
if (slope .ge. 0.d0) then
	print *, 'roundoff problem in lnsrch'
	print *, 'slope:',slope
	stop
endif
!Compute λmin.
test=0.d0 
do i=1,n
	temp=abs(p(i))/max(abs(xold(i)),1.d0)
	if(temp.gt.test)test=temp 
enddo
alamin=TOLX/test
!Always try full Newton step first.
alam=1.d0 
1 continue
	do i=1,n
		x(i)=xold(i)+alam*p(i) 
	enddo
	f=func(n,x,graph_number) 
	!Convergence on ∆x. For zero finding, the calling program should verify the convergence.
	if(alam.lt.alamin)then
		do i=1,n 
			x(i)=xold(i)
		enddo
		check=.true. 
		return
	!Sufficient function decrease.
	else if(f.le.fold+ALF*alam*slope)then 
		return
	!Backtrack. 
	else 
		!First time.
		if(alam.eq.1.d0)then
			tmplam=-slope/(2.d0*(f-fold-slope)) 
		!Subsequent backtracks.
		else
			rhs1=f-fold-alam*slope 
			rhs2=f2-fold-alam2*slope 
			a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2) 
			b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
			if(a.eq.0.d0)then
				tmplam=-slope/(2.d0*b) 
			else
				disc=b*b-3.d0*a*slope 
				if (disc .lt. 0.d0) then 
					tmplam=0.5d0*alam
				else if(b.le.0.d0)then 
					tmplam=(-b+sqrt(disc))/(3.d0*a)
				else 
					tmplam=-slope/(b+sqrt(disc))
				endif 
			endif
			if(tmplam.gt.0.5d0*alam)tmplam=0.5d0*alam
		endif
	endif
	alam2=alam
	f2=f
	alam=max(tmplam,0.1d0*alam)
goto 1
end subroutine lnsrch



subroutine dfunc(n,x,grad,func,graph_number)
!from Numerical Recipes in FORTRAN 77: The art of Scientific Programming
!https://www.goodreads.com/book/show/27815.Numerical_Recipes_in_FORTRAN_77
	implicit none

	integer :: n,i,graph_number
	double precision :: x(n), x_dif_plus(n), x_dif_min(n)
	double precision :: nominal_delta = 1.d-6
	double precision :: h,temp
	double precision :: grad(n)
	double precision :: func

	External func

	!set h exact to numerical precision
	!so h can be represented exactly in base 2
	temp =  x(1) + nominal_delta
 	h = temp-x(1)

 	do i = 1,n
  		x_dif_plus = x
 		x_dif_min = x
 		x_dif_plus(i) = x(i) + h
 		x_dif_min(i) = x(i) - h
 		grad(i) = (func(n,x_dif_plus,graph_number) - func(n,x_dif_min,graph_number))/(2.d0*h)
	enddo

end subroutine dfunc


end module BFGS