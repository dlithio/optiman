module nse_mod_physf
use utility_mod
use projector_mod_phys
implicit none
integer, private, parameter :: WP = selected_real_kind(15,307)
integer, private :: M
integer, private :: kxstar,kystar
integer, private :: auto_ndim
complex(WP), private :: e(2)
double precision, allocatable, private :: beta(:,:,:,:)
complex(WP), private, allocatable :: omega_hat_c(:,:)
complex(WP), allocatable, private :: nl_e_omega(:,:)
complex(WP), allocatable, private :: nl_omega_e(:,:)
complex(WP), allocatable, private :: jac_temp(:,:)
complex(WP), allocatable, private :: forcing(:,:)
complex(WP), allocatable, private :: nl(:,:)
complex(WP), allocatable, private :: linear(:,:)
complex(WP), allocatable, private :: two_d_dummy(:,:)
real(WP), private :: viscosity
real(WP), private, allocatable :: dfdp_mult(:)


contains
subroutine nse_mod_physf_init(m_input,kxstar_in,kystar_in,auto_ndim_input,forcing_input)
implicit none
integer, intent(in) :: m_input
integer, intent(in) :: auto_ndim_input
integer, intent(in) :: kxstar_in, kystar_in
integer :: ix,iy,jx,jy,kx,ky
integer :: mode_num
complex(WP), intent(inout) :: forcing_input( 0: , -m_input/2: )
M = m_input
e(1) = dcmplx(1.d0,0.d0)
e(2) = dcmplx(0.d0,1.d0)
allocate(beta(-M/2:M/2,-M/2:M/2,-M/2:M/2,-M/2:M/2))
allocate(forcing(0:M/2,-M/2:M/2))
allocate(nl(0:M/2,-M/2:M/2))
allocate(linear(0:M/2,-M/2:M/2))
forcing = forcing_input
beta = 0.d0
auto_ndim = auto_ndim_input
do concurrent(ix=-M/2:M/2,iy=-M/2:M/2,jx=-M/2:M/2,jy=-M/2:M/2)
kx = ix+jx
ky = iy+jy
if (((ix**2+iy**2).ne. 0) .and. ((jx**2+jy**2).ne. 0) .and. (((ix+jx)**2+(iy+jy)**2).ne. 0)) then
beta(ix,iy,jx,jy) = dble((iy*jx - ix*jy)*(jy*ky + jx*kx))
beta(ix,iy,jx,jy) = beta(ix,iy,jx,jy)/&
                    dble(sqrt(dble(ix**2+iy**2))*sqrt(dble(jx**2+jy**2))*sqrt(dble(kx**2+ky**2)))
endif
enddo
allocate(omega_hat_c(-M/2:M/2,-M/2:M/2))
allocate(dfdp_mult(1:auto_ndim))
allocate(nl_e_omega(0:M/2,-M/2:M/2))
allocate(nl_omega_e(0:M/2,-M/2:M/2))
allocate(jac_temp(0:M/2,-M/2:M/2))
allocate(two_d_dummy(0:M/2,-M/2:M/2))
kxstar=kxstar_in
kystar=kystar_in
dfdp_mult = 0.d0
do mode_num=1,auto_ndim
    kx = kx_projections(mode_num)
    ky = ky_projections(mode_num)
    if ((kx .eq. kxstar) .and. (ky .eq. kystar)) then
        dfdp_mult(mode_num) = 1.d0
    endif
enddo
end subroutine

subroutine nse_phys(two_d_result,omega_hat,PAR_input)
implicit none
complex(WP), intent(inout) :: two_d_result( 0: , -M/2: )
complex(WP), intent(inout) :: omega_hat( 0: , -M/2: )
double precision, intent(inout) :: PAR_input(1:)
integer :: kx
integer :: ky
integer :: par_size
par_size = 36
call check_dims_real_1(PAR_input,par_size)
call update_forcing(PAR_input)
viscosity = 1.d0
call check_dims_complex(two_d_result,M/2+1,M+1)
call check_dims_complex(omega_hat,M/2+1,M+1)
two_d_result = dcmplx(0.d0,0.d0)
call get_nl_single(omega_hat,nl)
linear = dcmplx(0.d0,0.d0)
do concurrent (kx=0:M/2,ky=-M/2:M/2)
	linear(kx,ky) = omega_hat(kx,ky)*viscosity*(dble(kx**2+ky**2))
enddo
two_d_result = forcing - nl - linear
end subroutine nse_phys

subroutine update_forcing(PAR_input)
implicit none
double precision, intent(inout) :: PAR_input(1:)
integer :: ktest
integer :: par_size
par_size = 36
call check_dims_real_1(PAR_input,par_size)
ktest = kystar*kxstar-kxstar*kystar+1-kxstar-kystar
if (mod(ktest,2) .eq. 0) then
    forcing(kxstar,kystar) = dcmplx(PAR_input(1),0.d0)
else
    forcing(kxstar,kystar) = dcmplx(0.d0,PAR_input(1))
endif
end subroutine update_forcing 

subroutine get_nl_single(omega_hat,nl)
implicit none
complex(WP), intent(inout) :: omega_hat( 0: , -M/2: )
complex(WP), intent(inout) :: nl( 0: , -M/2: )
integer :: kx1
integer :: ky1
integer :: kx2
integer :: ky2
call check_dims_complex(omega_hat,M/2+1,M+1)
call check_dims_complex(nl,M/2+1,M+1)
nl = dcmplx(0.d0,0.d0)
call symmetrize(omega_hat)
call load_little_to_c(omega_hat,omega_hat_c)
do concurrent (kx1=1:M/2,ky1=-M/2:M/2,kx2=-kx1:min(M/2-kx1,M/2),ky2=max(-M/2-ky1,-M/2):min(M/2-ky1,M/2))
    nl(kx1+kx2,ky1+ky2) = nl(kx1+kx2,ky1+ky2) + &
                    omega_hat_c(kx1,ky1)*omega_hat_c(kx2,ky2)*beta(kx1,ky1,kx2,ky2) 
end do
do concurrent (kx1=-M/2:-1,ky1=-M/2:M/2,kx2=-kx1:min(M/2-kx1,M/2),ky2=max(-M/2-ky1,-M/2):min(M/2-ky1,M/2))
    nl(kx1+kx2,ky1+ky2) = nl(kx1+kx2,ky1+ky2) + &
                    omega_hat_c(kx1,ky1)*omega_hat_c(kx2,ky2)*beta(kx1,ky1,kx2,ky2) 
end do
kx1=0
do concurrent (ky1=-M/2:-1,kx2=-kx1:min(M/2-kx1,M/2),ky2=max(-M/2-ky1,-M/2):min(M/2-ky1,M/2))
    nl(kx1+kx2,ky1+ky2) = nl(kx1+kx2,ky1+ky2) + &
                    omega_hat_c(kx1,ky1)*omega_hat_c(kx2,ky2)*beta(kx1,ky1,kx2,ky2) 
end do
kx1=0
do concurrent (ky1=1:M/2,kx2=-kx1:min(M/2-kx1,M/2),ky2=max(-M/2-ky1,-M/2):min(M/2-ky1,M/2))
    nl(kx1+kx2,ky1+ky2) = nl(kx1+kx2,ky1+ky2) + &
                    omega_hat_c(kx1,ky1)*omega_hat_c(kx2,ky2)*beta(kx1,ky1,kx2,ky2) 
end do
nl = dcmplx(0.d0,1.d0) * nl
call symmetrize(nl)
end subroutine get_nl_single

subroutine get_nl_e_omega(input_e_hat,input_omega_hatc,result_array,kx_unit,ky_unit)
implicit none
complex(WP), intent(inout) :: input_e_hat
complex(WP), intent(inout) :: input_omega_hatc( -M/2: , -M/2: )
complex(WP), intent(inout) :: result_array( 0: , -M/2: )
integer, intent(inout) :: kx_unit
integer, intent(inout) :: ky_unit
integer :: kx1,ky1,kx2,ky2
call check_dims_complex(input_omega_hatc,M+1,M+1)
call check_dims_complex(result_array,M/2+1,M+1)
result_array = dcmplx(0.d0,0.d0)
kx1 = kx_unit
ky1 = ky_unit
do concurrent (kx2=-kx1:min(M/2-kx1,M/2),ky2=max(-M/2-ky1,-M/2):min(M/2-ky1,M/2))
    result_array(kx1+kx2,ky1+ky2) = result_array(kx1+kx2,ky1+ky2) + &
                    input_e_hat*omega_hat_c(kx2,ky2)*beta(kx1,ky1,kx2,ky2) 
end do
kx1 = -kx_unit
ky1 = -ky_unit
do concurrent (kx2=-kx1:min(M/2-kx1,M/2),ky2=max(-M/2-ky1,-M/2):min(M/2-ky1,M/2))
    result_array(kx1+kx2,ky1+ky2) = result_array(kx1+kx2,ky1+ky2) - &
                    dconjg(input_e_hat)*omega_hat_c(kx2,ky2)*beta(kx1,ky1,kx2,ky2) 
end do
result_array = dcmplx(0.d0,1.d0) * result_array
call symmetrize(result_array)
end subroutine get_nl_e_omega

subroutine get_nl_omega_e(input_e_hat,input_omega_hatc,result_array,kx_unit,ky_unit)
implicit none
complex(WP), intent(inout) :: input_e_hat
complex(WP), intent(inout) :: input_omega_hatc( -M/2: , -M/2: )
complex(WP), intent(inout) :: result_array( 0: , -M/2: )
integer, intent(inout) :: kx_unit
integer, intent(inout) :: ky_unit
integer :: kx1,ky1,kx2,ky2
call check_dims_complex(input_omega_hatc,M+1,M+1)
call check_dims_complex(result_array,M/2+1,M+1)
result_array = dcmplx(0.d0,0.d0)
kx1 = kx_unit
ky1 = ky_unit
do concurrent (kx2=-kx1:min(M/2-kx1,M/2),ky2=max(-M/2-ky1,-M/2):min(M/2-ky1,M/2))
    result_array(kx1+kx2,ky1+ky2) = result_array(kx1+kx2,ky1+ky2) + &
                    input_e_hat*omega_hat_c(kx2,ky2)*beta(kx2,ky2,kx1,ky1) 
end do
kx1 = -kx_unit
ky1 = -ky_unit
do concurrent (kx2=-kx1:min(M/2-kx1,M/2),ky2=max(-M/2-ky1,-M/2):min(M/2-ky1,M/2))
    result_array(kx1+kx2,ky1+ky2) = result_array(kx1+kx2,ky1+ky2) - &
                    dconjg(input_e_hat)*omega_hat_c(kx2,ky2)*beta(kx2,ky2,kx1,ky1) 
end do
result_array = dcmplx(0.d0,1.d0) * result_array
call symmetrize(result_array)
end subroutine get_nl_omega_e

subroutine get_jacobian_from_2d(omega_hat,PAR_input,jac)
implicit none
complex(WP), intent(inout) :: omega_hat( 0: , -M/2: )
double precision, intent(inout) :: PAR_input(1:)
real(WP), intent(inout) :: jac( 1: , 1: )
complex(WP) :: e_hat
integer :: kx_unit
integer :: ky_unit
integer :: mode_num
integer :: j
real(WP) :: viscosity
integer :: par_size
integer :: ktest
par_size = 36
jac = 0.d0
call check_dims_real_1(PAR_input,par_size)
call check_dims_complex(omega_hat,M/2+1,M+1)
call check_dims_real_2(jac,auto_ndim,auto_ndim)
call symmetrize(omega_hat)
call load_little_to_c(omega_hat,omega_hat_c)
viscosity = 1.d0
! We will loop through each mode that we project into 1d
do mode_num = 1, auto_ndim
    kx_unit = kx_projections(mode_num)
    ky_unit = ky_projections(mode_num)
    ! Now make that mode either 1+0i or 0+1i
    ktest = ky_unit*kxstar-kx_unit*kystar+1-kx_unit-ky_unit
    if (mod(ktest,2) .eq. 0) then
        j = 1
    else
        j = 2
    endif
    ! In this inner most loop, we're getting a column of 
    ! the jacobian
    e_hat = e(j)
    ! Get B(w,e_i) and B(e_i,w)
    call get_nl_e_omega(e_hat,omega_hat_c,nl_e_omega,kx_unit,ky_unit)
    call get_nl_omega_e(e_hat,omega_hat_c,nl_omega_e,kx_unit,ky_unit)
    ! Get df_i/de_i from the nl part (a column of the jacobian but in 2d)
    jac_temp = nl_e_omega(0:M/2,-M/2:M/2) + nl_omega_e(0:M/2,-M/2:M/2)
    ! Get df_i/de_i from the linear part
    jac_temp(kx_unit,ky_unit) = jac_temp(kx_unit,ky_unit) + &
        viscosity*dble(kx_unit**2+ky_unit**2)*e_hat
    if (kx_unit .eq. 0) then
        jac_temp(kx_unit,-ky_unit) = -dconjg(jac_temp(kx_unit,ky_unit))
    endif
    jac_temp = -jac_temp
    ! Project the result into the actual jacobian matrix
    call project_to_jac(jac_temp,jac,mode_num)
end do
end subroutine get_jacobian_from_2d

subroutine get_jacobian_from_1d(auto_state,PAR_input,jacobian_return)
implicit none
double precision, intent(inout) :: PAR_input(1:)
real(WP), intent(inout) :: auto_state( 1: )
real(WP), intent(inout) :: jacobian_return( 1: , 1: )
integer :: par_size
par_size = 36
jacobian_return = 0.d0
call check_dims_real_1(PAR_input,par_size)
call check_dims_real_1(auto_state,auto_ndim)
call check_dims_real_2(jacobian_return,auto_ndim,auto_ndim)
call project_to_2d(auto_state,two_d_dummy)
call get_jacobian_from_2d(two_d_dummy,PAR_input,jacobian_return)
end subroutine get_jacobian_from_1d

subroutine symmetrize(little)
implicit none
complex(WP), intent(inout) :: little( 0: , -M/2: )
integer :: ky
real(WP) :: biggest
call check_dims_complex(little,M/2+1,M+1)
biggest = maxval(abs(little))
do ky=1,M/2
    if (abs( little(0,-ky) + conjg(little(0,ky)) ) .gt. (biggest * 1.d-10)) then
        write(*,*) ky
        write(*,*) little(0,-ky)
        write(*,*) little(0,ky)
        write(*,*) (biggest * 1.d-10)
        write(*,*) "symmetry has been broken"
        stop
    endif
    little(0,-ky) = (-dconjg(little(0,ky)) + little(0,-ky))/2.d0
    little(0,ky) = -dconjg(little(0,-ky))
enddo
if (abs(little(0,0)) .gt. biggest*1.d-14) then
    write(*,*) "0 mean has been broken"
    stop
endif
little( 0, 0) = dcmplx(0.d0,0.d0)
end subroutine symmetrize

subroutine load_little_to_c(little,c_array)
implicit none
complex(WP), intent(inout) :: little( 0: , -M/2: )
complex(WP), intent(inout) :: c_array( -M/2: , -M/2: )
integer :: kx
integer :: ky
call check_dims_complex(little,M/2+1,M+1)
call check_dims_complex(c_array,M+1,M+1)
c_array = dcmplx(0.d0,0.d0)
c_array( 0:M/2 , -M/2:M/2 ) = little( 0:M/2, -M/2:M/2 ) 
do kx = -M/2,-1
    do ky = -M/2,M/2
        c_array(kx,ky) = -dconjg(little(-kx,-ky))
    enddo
enddo
end subroutine load_little_to_c

subroutine get_dfdp_from_1d(auto_state,PAR_input,dfdp_return)
implicit none
double precision, intent(inout) :: PAR_input(1:)
real(WP), intent(inout) :: auto_state( 1: )
real(WP), intent(inout) :: dfdp_return( 1: , 1: )
integer :: par_size
dfdp_return = 0.d0
par_size = 36
call check_dims_real_1(PAR_input,par_size)
call check_dims_real_1(auto_state,auto_ndim)
call check_dims_real_2(dfdp_return,auto_ndim,par_size)
dfdp_return(:,1) = dfdp_mult
end subroutine get_dfdp_from_1d

end module nse_mod_physf
