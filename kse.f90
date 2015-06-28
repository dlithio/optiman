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
allocate(uu(3*int(par(36))))
allocate(b(int(par(36))))
end subroutine setup

subroutine fcn ( n, x, fvec)
implicit none
integer, intent(in) :: n
double precision, intent(inout) :: x(n) ! really intent(in), want to avoid a mem copy though
double precision, intent(inout) :: fvec(n)
integer :: j
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

subroutine get_jac(n,x,a)
implicit none
integer, intent(inout) :: n
double precision, intent(inout) :: x( 1: )
double precision, intent(inout) :: a( 1: , 1: )
double precision :: original_force(n),forward_state(n),forward_force(n)
integer :: i
double precision :: h
h = 1.d-4
a = 0.d0
call fcn ( n, x, original_force)
do i=1,n
    forward_state = x
    forward_state(i) = forward_state(i) + h
    call fcn ( n, forward_state, forward_force)
    forward_force = (forward_force - original_force)/h
    a(:,i) = forward_force
enddo
end subroutine get_jac

end module
