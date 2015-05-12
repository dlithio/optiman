subroutine fcn ( n, x, fvec, iflag)
implicit none
integer, intent(in) :: n
double precision, intent(inout) :: x(n) ! really intent(in), want to avoid a mem copy though
double precision, intent(inout) :: fvec(n)
integer, intent(inout) :: iflag

end subroutine fcn
