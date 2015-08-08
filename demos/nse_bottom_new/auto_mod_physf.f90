module auto_mod_physf
use nse_mod_physf
use projector_mod_phys
use utility_mod
implicit none
integer, private, parameter :: WP = selected_real_kind(15,307)
real(WP), private :: par1,ar1,ai1,ar2,ai2
integer, private :: forced_modes,kx1,ky1,kx2,ky2
complex(WP), allocatable, private :: state_array_2d(:,:)
complex(WP), allocatable, private :: result_array_2d(:,:)
integer, private :: M
integer, private :: auto_ndim
real(WP), allocatable, private :: energy_mult(:)
real(WP), allocatable, private :: palinstrophy_mult(:)
real(WP), allocatable, private :: energy_dissip_mult(:)
real(WP), allocatable, private :: energy_multr(:)
real(WP), allocatable, private :: palinstrophy_multr(:)
real(WP), allocatable, private :: energy_dissip_multr(:)
real(WP), allocatable, private :: nl_1d(:)
real(WP), allocatable, private :: linear_1d(:)
real(WP), allocatable, private :: result_1d(:)
complex(WP), allocatable, private :: nl_2d(:,:)
complex(WP), allocatable, private :: linear_2d(:,:)
real(WP), private :: energy_of_force
contains

subroutine auto_mod_physf_init(m_input,auto_ndim_input,par1init,forced_modesinit,&
                        kx1init,ky1init,ar1init,ai1init,kx2init,ky2init,&
                        ar2init,ai2init)
implicit none
integer, intent(in) :: m_input
integer, intent(in) :: auto_ndim_input
real(WP), intent(in) :: par1init,ar1init,ai1init,ar2init,ai2init
integer, intent(in) :: forced_modesinit,kx1init,ky1init,kx2init,ky2init
integer :: mode_num
integer :: kx
integer :: ky
M = m_input
auto_ndim = auto_ndim_input
par1=par1init
ar1=ar1init
ai1=ai1init
ar2=ar2init
ai2=ai2init
forced_modes=forced_modesinit
kx1=kx1init
ky1=ky1init
kx2=kx2init
ky2=ky2init
allocate(state_array_2d(0:M/2,-M/2:M/2))
allocate(result_array_2d(0:M/2,-M/2:M/2))
allocate(nl_2d(0:M/2,-M/2:M/2))
allocate(linear_2d(0:M/2,-M/2:M/2))
allocate(nl_1d(auto_ndim))
allocate(linear_1d(auto_ndim))
allocate(energy_mult(auto_ndim))
allocate(palinstrophy_mult(auto_ndim))
allocate(energy_dissip_mult(auto_ndim))
allocate(energy_multr(auto_ndim))
allocate(palinstrophy_multr(auto_ndim))
allocate(energy_dissip_multr(auto_ndim))
allocate(result_1d(auto_ndim))
energy_mult = 0.d0
palinstrophy_mult = 0.d0
energy_dissip_mult = 0.d0
energy_multr = 0.d0
palinstrophy_multr = 0.d0
energy_dissip_multr = 0.d0
do mode_num=1,auto_ndim
    kx = kx_projections(mode_num)
    ky = ky_projections(mode_num)
    energy_mult(mode_num) = 1.d0/dble(kx**2+ky**2)
    palinstrophy_mult(mode_num) = dble(kx**2+ky**2)
    energy_dissip_mult(mode_num) = (dble(kx**2-ky**2)/dble(kx**2+ky**2))**2
enddo
end subroutine auto_mod_physf_init

subroutine auto_state_to_force(state,return_force,PAR_input)
implicit none
real(WP), intent(inout) :: state( 1: )
real(WP), intent(inout) :: return_force( 1: )
double precision, intent(inout) :: PAR_input(1:)
integer :: par_size
par_size = 36
return_force = 0.d0
call check_dims_real_1(PAR_input,par_size)
call check_dims_real_1(state,auto_ndim)
call check_dims_real_1(return_force,auto_ndim)
call project_to_2d(state,state_array_2d)
call nse_phys(result_array_2d,state_array_2d,PAR_input)
call project_to_1d(result_array_2d,return_force)
end subroutine auto_state_to_force

subroutine write_to_two_d(kx,ky,coefficent,two_d_array)
!-----------------------------------------------------------------------
! subroutine: write_two_d_to_one_d (replaces write_w_to_F)
! Given an unshifted mode, writes to the 2d array used in fortran
!
! Inputs:
! kx : the unshifted x-mode
! ky : the unshifted y-mode
! coefficent : the complex coefficient to write to that mode
! ndim : The dimesion of the one-dimesional array
! n : relates to the size of the two_d array (for NSE calculations)
!
! Outputs:
! two_d_array : the two-dimesional complex array that we're writing to
!               (for NSE calculations)
!-----------------------------------------------------------------------
implicit none
integer, intent(in) :: kx
integer, intent(in) :: ky
complex(WP), intent(in) :: coefficent
complex(WP), intent(inout) :: two_d_array( 0: , -M/2: )
call check_dims_complex(two_d_array,M/2+1,M+1)
if ((kx .gt. M/2) .or. (kx .lt. -M/2) .or. (ky .gt. M/2) .or. (ky .lt. -M/2)) then
    write(*,*) "M ",M
    write(*,*) "kx ",kx
    write(*,*) "ky ",ky
    write(*,*) "illegal write to 2d"
    stop
endif
if ((kx .eq. 0) .and. (ky .eq. 0)) then
    write(*,*) "illegal write to (0,0)"
    stop
endif
if (kx .gt. 0) then
    two_d_array(kx,ky) = coefficent
endif
if ((kx .eq. 0) .and. (ky .ne. 0)) then
    two_d_array(kx,ky) = coefficent
    two_d_array(-kx,-ky) = -dconjg(coefficent)
endif
if (kx .lt. 0) then
    two_d_array(-kx,-ky) = -dconjg(coefficent)
endif
end subroutine write_to_two_d

subroutine get_external_force_and_steady_state(external_force,steady_state,PAR_input)
!-----------------------------------------------------------------------
! subroutine: single_mode_external_force
! This subroutine is used for getting the desired external force array
! for a given problem and its steady state.
! 
! Outputs:
! external_force : a one-dimensional array that is the desired external 
!                  force (for use with AUTO)
! steady_state : a one-dimensional array that is the desired steady 
!                state (for use with AUTO)
!-----------------------------------------------------------------------
implicit none
integer, parameter :: real8   = selected_real_kind(15,307)
complex(WP), intent(inout) :: external_force( 0: , -M/2: )
complex(WP), intent(inout) :: steady_state( 0: , -M/2: )
integer :: num_modes
integer :: i
integer :: j
integer :: kx_1
integer :: ky_1
integer :: kx_high
integer :: ky_high
integer :: kx_low
integer :: ky_low
integer :: ktest
double precision :: force_r_1
double precision :: force_i_1
double precision :: alpha_r_high
double precision :: alpha_i_high
double precision :: alpha_r_low
double precision :: alpha_i_low
double precision :: k_high_times_k_low_times_diff
double precision :: v
complex(kind=real8), dimension(0:M/2,-M/2:M/2) :: external_force_two_d
complex(kind=real8), dimension(0:M/2,-M/2:M/2) :: external_force_two_d_temp
complex(kind=real8), dimension(0:M/2,-M/2:M/2) :: steady_state_two_d
complex(kind=real8), dimension(0:M/2,-M/2:M/2) :: steady_state_two_d_temp
complex(kind=real8), dimension(8) :: check_for_distinct
complex(kind=real8) :: coefficent
double precision, intent(inout) :: PAR_input(1:)
integer :: par_size
par_size = 36
call check_dims_real_1(PAR_input,par_size)
! Initialize Variables
call check_dims_complex(external_force,M/2+1,M+1)
call check_dims_complex(steady_state,M/2+1,M+1)
external_force = dcmplx(0.d0,0.d0)
steady_state = dcmplx(0.d0,0.d0)
num_modes = 0
i=0
j=0
kx_1 = 0
ky_1 = 0
kx_high = 0
ky_high = 0
kx_low = 0
ky_low = 0
force_r_1 = 0.d0
force_i_1 = 0.d0
alpha_r_high= 0.d0
alpha_i_high = 0.d0
alpha_r_low= 0.d0
alpha_i_low = 0.d0
k_high_times_k_low_times_diff = 0.d0
v=0.d0
external_force_two_d = dcmplx(0.d0,0.d0)
steady_state_two_d = dcmplx(0.d0,0.d0)
external_force_two_d_temp = dcmplx(0.d0,0.d0)
steady_state_two_d_temp = dcmplx(0.d0,0.d0)
coefficent = dcmplx(0.d0,0.d0)
check_for_distinct=dcmplx(0.d0,0.d0)
! Read in PAR
v = 1.d0
num_modes = forced_modes
kx_1 = kx1
ky_1 = ky1
force_r_1 = ar1
force_i_1 = ai1
! Use ktest to decide whether to write to imaginary or real
ktest = ky1*kx1-kx1*ky1+1-kx1-ky1
! If we want a single mode forced...
if (num_modes .eq. 1) then
    if (mod(ktest,2) .eq. 0) then
        call write_to_two_d(kx_1,ky_1, & 
                    dcmplx(PAR_input(1),0.d0), &
                    external_force_two_d)
        coefficent = dcmplx(PAR_input(1)/dble(kx_1**2+ky_1**2), &
                            0.d0)
        call write_to_two_d(kx_1,ky_1, & 
                            coefficent, &
                            steady_state_two_d)
    else
        call write_to_two_d(kx_1,ky_1, & 
                    dcmplx(0.d0,PAR_input(1)), &
                    external_force_two_d)
        coefficent = dcmplx(0.d0, &
                            PAR_input(1)/dble(kx_1**2+ky_1**2))
        call write_to_two_d(kx_1,ky_1, & 
                            coefficent, &
                            steady_state_two_d)
    endif
end if
! If we want a two mode stead state...
if (num_modes .eq. 2) then
    write(*,*) "can only do 1 mode forcing with the physcial case"
    stop
end if
! Copy results to output
external_force = external_force_two_d
steady_state = steady_state_two_d
end subroutine get_external_force_and_steady_state

subroutine get_extra_solution_measures(state,PAR_input)
implicit none
real(WP), intent(inout) :: state( 1: )
double precision, intent(inout) :: PAR_input(1:)
integer :: par_size
par_size = 36
call check_dims_real_1(PAR_input,par_size)
call check_dims_real_1(state,auto_ndim)
PAR_input(15) = sum(state*state*energy_mult)
PAR_input(16) = sum(state*state*palinstrophy_mult)
PAR_input(17) = sum(state*state*energy_dissip_mult)/1.d0
!PAR_input(18) = (2.d0*energy_of_force * ((PAR_input(1))**4))**0.5
PAR_input(18) = 0.d0
PAR_input(19) = state(1)
PAR_input(20) = 0
PAR_input(21) = 0
end subroutine get_extra_solution_measures

subroutine save_energy_of_force(force)
implicit none
complex(WP), intent(inout) :: force( 0: , -M/2: )
real(WP) :: one_d_array( 1:auto_ndim )
call check_dims_complex(force,M/2+1,M+1)
call project_to_1d(force,one_d_array)
energy_of_force = sum(one_d_array*one_d_array)
end subroutine save_energy_of_force

subroutine cleanup_auto_mod_physf(err_num)
implicit none
integer, intent(inout) :: err_num
! Return same number of errors that was given as input. Errors in this
! module stop the program
err_num = err_num
deallocate(state_array_2d)
deallocate(result_array_2d)
par1 = 0.d0
ar1 = 0.d0
ai1 = 0.d0
ar2 = 0.d0
ai2 = 0.d0
forced_modes = 0
kx1 = 0
kx2 = 0
ky1 = 0
ky2 = 0
M=0
auto_ndim =0
energy_of_force = 0.d0
deallocate(energy_mult)
deallocate(palinstrophy_mult)
deallocate(energy_dissip_mult)
deallocate(energy_multr)
deallocate(palinstrophy_multr)
deallocate(energy_dissip_multr)
deallocate(result_1d)
deallocate(nl_2d)
deallocate(linear_2d)
deallocate(nl_1d)
deallocate(linear_1d)
end subroutine cleanup_auto_mod_physf
end module auto_mod_physf
