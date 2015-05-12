module user_functions
implicit none
contains

subroutine fcn ( n, x, fvec, iflag)
implicit none
integer, intent(in) :: n
double precision, intent(inout) :: x(n) ! really intent(in), want to avoid a mem copy though
double precision, intent(inout) :: fvec(n)
integer, intent(inout) :: iflag
iflag = 0
fvec(1) = x(1)
fvec(2) = x(2)
fvec(3) = -x(3)
end subroutine fcn
end module
