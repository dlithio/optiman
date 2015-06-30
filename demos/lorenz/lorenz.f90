module user_functions
implicit none
contains

subroutine setup (par)
implicit none
double precision, intent(inout) :: par(:)
end subroutine setup

subroutine fcn ( n, x, fvec)
implicit none
integer, intent(in) :: n
double precision, intent(inout) :: x(n) ! really intent(in), want to avoid a mem copy though
double precision, intent(inout) :: fvec(n)
! The lorenz manifold
fvec(1)=10.d0*(x(2)-x(1));
fvec(2)=28.d0*x(1)-x(2)-x(1) * x(3);
fvec(3)=x(1)*x(2) - 8.d0/3.d0 * x(3);
! We are calculating the stable manifold, so we reverse the vector field
fvec = -1.d0*fvec
end subroutine fcn

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
