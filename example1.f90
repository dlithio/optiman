module user_functions
implicit none
double precision, allocatable, private :: uu(:)
double precision, allocatable, private :: b(:)
double precision, private :: alpha
contains

subroutine setup (par)
implicit none
double precision, intent(inout) :: par(:)
alpha = par(1)
allocate(uu(3*int(par(2))))
allocate(b(int(par(2))))
end subroutine setup

subroutine fcn ( n, x, fvec, iflag)
implicit none
integer, intent(in) :: n
double precision, intent(inout) :: x(n) ! really intent(in), want to avoid a mem copy though
double precision, intent(inout) :: fvec(n)
integer, intent(inout) :: iflag
integer :: j
iflag = 0
! A test manifold
!fvec(1) = x(1)
!fvec(2) = x(2)
!fvec(3) = -x(3)

! The lorenz manifold
!fvec(1)=10.d0*(x(2)-x(1));
!fvec(2)=28.d0*x(1)-x(2)-x(1) * x(3);
!fvec(3)=x(1)*x(2) - 8.d0/3.d0 * x(3);
!fvec = -1.d0*fvec

! The kse
call getBsingle(n,x,b)
do j=1,n
    fvec(j)=dfloat(-4*j**4)*x(j)+alpha*dfloat(j**2)*x(j) - b(j)
enddo
end subroutine fcn

subroutine getBsingle(n,u,b_out)
! Computes the bilinear term of the KSE
! Does not include the extra high modes.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT none
integer, intent(in) :: n
double precision, intent(inout) :: u(1:)
double precision, intent(inout) :: b_out(1:)
DOUBLE PRECISION :: mysum
integer :: k,kjsgn,j
uu = 0.d0
uu(1:n) = u
b_out = 0.d0
do  k=1,n
    mysum=0.d0
    do  j=1,n
        if (k .gt. j) then
            kjsgn=1
            mysum=mysum+j*u(j)*(uu(j+k)+kjsgn*uu(abs(k-j)))
        else if (k .lt. j) then
            kjsgn=-1
            mysum=mysum+j*u(j)*(uu(j+k)+kjsgn*uu(abs(k-j)))
        else
            kjsgn=0
            mysum=mysum+j*u(j)*uu(j+k)
        endif
    enddo
    b_out(k)=alpha*mysum/2
enddo
end subroutine getBsingle

end module
