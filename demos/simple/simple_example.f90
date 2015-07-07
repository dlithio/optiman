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
! A test manifold
fvec(1) = 2.d0*x(1)
fvec(2) = 2.d0*x(2)
fvec(3) = 2.d0*x(1)*dcos(x(1)) + 2.d0*x(2)*dcos(x(2)) - x(3) + dsin(x(1)) + dsin(x(2))
end subroutine fcn

subroutine get_jac(n,x,a)
implicit none
integer, intent(inout) :: n
double precision, intent(inout) :: x( 1: )
double precision, intent(inout) :: a( 1: , 1: )
double precision :: original_force(n),forward_state(n),forward_force(n)
integer :: i
a = 0.d0
a(1,1) = 2.d0
a(2,2) = 2.d0
a(3,3) = -1.d0
a(3,1) = 2.d0*dcos(x(1))-2.d0*x(1)*dsin(x(1))+dcos(x(1))
a(3,2) = 2.d0*dcos(x(2))-2.d0*x(2)*dsin(x(2))+dcos(x(2))
end subroutine get_jac

end module
