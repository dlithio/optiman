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

!fvec(1)=10.d0*(x(2)-x(1));
!fvec(2)=28.d0*x(1)-x(2)-x(1) * x(3);
!fvec(3)=x(1)*x(2) - 8.d0/3.d0 * x(3);

!fvec = -1.d0*fvec

end subroutine fcn
end module
