module projector_mod_phys
use utility_mod
implicit none
integer, private, parameter :: WP = selected_real_kind(15,307)
integer, private, parameter :: kdp = selected_real_kind(15)
integer, private :: M
integer, private :: auto_ndim
integer, private :: step_num
integer, private :: num_errors
integer, private :: kxstar,kystar
integer, allocatable :: kx_projections(:)
integer, allocatable :: ky_projections(:)
contains

subroutine projector_mod_phys_init(step_num_in,m_return,n_return,auto_ndim_return,kxstar_in,kystar_in)
implicit none
integer, intent(in) :: kxstar_in, kystar_in
integer, intent(in) :: step_num_in
integer, intent(out) :: m_return
integer, intent(out) :: n_return
integer, intent(out) :: auto_ndim_return
integer :: counter
kxstar=kxstar_in
kystar=kystar_in
num_errors = 0
step_num = step_num_in
! Determine the auto_ndim
call get_auto_ndim(step_num_in,auto_ndim_return)
auto_ndim_return = auto_ndim_return/2
auto_ndim = auto_ndim_return
! Determine the M
call get_M(step_num_in,m_return)
M = m_return
! Determine the N
call get_N(step_num_in,n_return)
! Set up the kx_projections and ky_projections
allocate(kx_projections(auto_ndim_return))
allocate(ky_projections(auto_ndim_return))
open(unit=200,file="kx_projections")
open(unit=201,file="ky_projections")
do counter = 1,auto_ndim_return
    read(200,*) kx_projections(counter)
    read(201,*) ky_projections(counter)
enddo
close(200)
close(201)
end subroutine projector_mod_phys_init

subroutine project_to_1d(nse_array,one_d_array)
implicit none
complex(WP), intent(inout) :: nse_array( 0: , -M/2: )
real(WP), intent(inout) :: one_d_array( 1: )
integer :: dim_num
integer :: kx
integer :: ky
integer :: ktest
call check_dims_complex(nse_array,M/2+1,M+1)
call check_dims_real_1(one_d_array,auto_ndim)
one_d_array = 0.d0
do dim_num = 1,auto_ndim
    kx = kx_projections(dim_num)
    ky = ky_projections(dim_num)
    ktest = ky*kxstar-kx*kystar+1-kx-ky
    if (mod(ktest,2) .eq. 0) then
        one_d_array(dim_num) = dreal(nse_array(kx,ky))
        if (abs(dimag(nse_array(kx,ky))) .gt. 1.d-14) then
            write(*,*) "SYMMETRY PROBLEM", dimag(nse_array(kx,ky))
        endif
    else
        one_d_array(dim_num) = dimag(nse_array(kx,ky))
        if ((dreal(nse_array(kx,ky))) .gt. 1.d-14) then
            write(*,*) "SYMMETRY PROBLEM", dreal(nse_array(kx,ky))
        endif
    endif
enddo
end subroutine project_to_1d

subroutine project_to_jac(nse_array,jac,col_num)
implicit none
complex(WP), intent(inout) :: nse_array( 0: , -M/2: )
real(WP), intent(inout) :: jac( 1: , 1: )
integer, intent(in) :: col_num
integer :: dim_num
integer :: kx
integer :: ky
integer :: ktest
call check_dims_complex(nse_array,M/2+1,M+1)
call check_dims_real_2(jac,auto_ndim,auto_ndim)
jac(:,col_num) = 0.d0
do dim_num = 1,auto_ndim
    kx = kx_projections(dim_num)
    ky = ky_projections(dim_num)
    ktest = ky*kxstar-kx*kystar+1-kx-ky
    if (mod(ktest,2) .eq. 0) then
        jac(dim_num,col_num) = dreal(nse_array(kx,ky))
    else
        jac(dim_num,col_num) = dimag(nse_array(kx,ky))
    endif
enddo
end subroutine project_to_jac

subroutine project_to_2d(one_d_array,nse_array)
implicit none
complex(WP), intent(inout) :: nse_array( 0: , -M/2: )
real(WP), intent(inout) :: one_d_array( 1: )
integer :: dim_num
integer :: kx
integer :: ky
integer :: ktest
call check_dims_complex(nse_array,M/2+1,M+1)
call check_dims_real_1(one_d_array,auto_ndim)
nse_array = dcmplx(0.d0,0.d0)
do dim_num = 1,auto_ndim
    kx = kx_projections(dim_num)
    ky = ky_projections(dim_num)
    ktest = ky*kxstar-kx*kystar+1-kx-ky
    if (mod(ktest,2) .eq. 0) then
        nse_array(kx,ky) = dcmplx(one_d_array(dim_num),0.d0)
    else
        nse_array(kx,ky) = dcmplx(0.d0,one_d_array(dim_num))
    endif
enddo
kx=0
do ky=1,M/2
    nse_array(kx,-ky) = -dconjg(nse_array(kx,ky))
enddo
nse_array(0,0) = dcmplx(0.d0,0.d0)
end subroutine project_to_2d

subroutine cleanup_projector_mod_phys(err_return)
!-----------------------------------------------------------------------
! subroutine: cleanup_projector_mod
! Returns the numbers of errors while the mod has been running and
! resets the necessary scalars/deallocates the arrays
!-----------------------------------------------------------------------
implicit none
integer, intent(inout) :: err_return
err_return = err_return + num_errors
num_errors = 0
M = 0
auto_ndim = 0
step_num = 0
deallocate(kx_projections)
deallocate(ky_projections)
end subroutine cleanup_projector_mod_phys

end module projector_mod_phys
