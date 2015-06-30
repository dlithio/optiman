module utility_mod
implicit none
integer, private :: num_errors
integer, private, parameter :: WP = selected_real_kind(15,307)
contains

subroutine utility_mod_init()
!-----------------------------------------------------------------------
! subroutine: utility_mod_init
! Sets the error count
!-----------------------------------------------------------------------
implicit none
num_errors = 0
end subroutine utility_mod_init

subroutine get_M(step_num_in,M)
!-----------------------------------------------------------------------
! subroutine: get_M
! Reads the file to get an M for a specific step
!-----------------------------------------------------------------------
implicit none
integer, intent(in) :: step_num_in
integer, intent(out) :: M
integer :: trash
integer :: counter
! The auto_ndim is the step_num-th line of the file auto_ndim_key
open(unit=200,file="M_key")
counter = 1
do while (counter .lt. step_num_in)
    read(200,*) trash
    counter = counter + 1 
enddo
read(200,*) M
close(200)
end subroutine get_M

subroutine get_N(step_num_in,N)
!-----------------------------------------------------------------------
! subroutine: get_N
! Reads the file to get an N for a specific step
!-----------------------------------------------------------------------
implicit none
integer, intent(in) :: step_num_in
integer, intent(out) :: N
integer :: trash
integer :: counter
! The auto_ndim is the step_num-th line of the file auto_ndim_key
open(unit=200,file="N_key")
counter = 1
do while (counter .lt. step_num_in)
    read(200,*) trash
    counter = counter + 1 
enddo
read(200,*) N
close(200)
end subroutine get_N

subroutine get_auto_ndim(step_num_in,auto_ndim)
!-----------------------------------------------------------------------
! subroutine: get_auto_ndim
! Reads the file to get an auto_ndim for a specific step
!-----------------------------------------------------------------------
implicit none
integer, intent(in) :: step_num_in
integer, intent(out) :: auto_ndim
integer :: trash
integer :: counter
! The auto_ndim is the step_num-th line of the file auto_ndim_key
open(unit=200,file="auto_ndim_key")
counter = 1
do while (counter .lt. step_num_in)
    read(200,*) trash
    counter = counter + 1 
enddo
read(200,*) auto_ndim
close(200)
end subroutine get_auto_ndim

subroutine check_dims_complex(array,first,second)
!-----------------------------------------------------------------------
! subroutine: check_dims_complex
! Checks the dimensions of a 2d complex array
!-----------------------------------------------------------------------
implicit none
complex(WP), intent(inout) :: array(:,:)
integer, intent(in) :: first
integer, intent(in) :: second
integer :: dims(2)
dims = shape(array)
if (dims(1) .ne. first) then
    num_errors = num_errors + 1
endif
if (dims(2) .ne. second) then
    num_errors = num_errors + 1
endif
end subroutine check_dims_complex

subroutine check_dims_real_1(array,first)
!-----------------------------------------------------------------------
! subroutine: check_dims_real_1
! Checks the dimensions of a 1d real array
!-----------------------------------------------------------------------
implicit none
real(WP), intent(inout) :: array(:)
integer, intent(in) :: first
integer :: dims(1)
dims = shape(array)
if (dims(1) .ne. first) then
    num_errors = num_errors + 1
endif
end subroutine check_dims_real_1

subroutine check_dims_real_2(array,first,second)
!-----------------------------------------------------------------------
! subroutine: check_dims_real_2
! Checks the dimensions of a 2d real array
!-----------------------------------------------------------------------
implicit none
real(WP), intent(inout) :: array(:,:)
integer, intent(in) :: first
integer, intent(in) :: second
integer :: dims(2)
dims = shape(array)
if (dims(1) .ne. first) then
    num_errors = num_errors + 1
endif
if (dims(2) .ne. second) then
    num_errors = num_errors + 1
endif
end subroutine check_dims_real_2

subroutine cleanup_utility_mod(err_return)
!-----------------------------------------------------------------------
! subroutine: cleanup_utility_mod
! Clears the number of errors
!-----------------------------------------------------------------------
implicit none
integer, intent(inout) :: err_return
err_return = err_return + num_errors
num_errors = 0
end subroutine cleanup_utility_mod

end module utility_mod
