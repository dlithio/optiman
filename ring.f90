module ring
use FFTW3
implicit none
double precision, private :: pi=4.d0*datan(1.d0)
double precision, private :: max_dist
integer, private :: sdiff_switch
integer, private :: t_switch
integer, private :: ts_switch
integer, private :: integral_switch
integer :: big
double precision, allocatable, private :: points(:,:)
double precision, allocatable, private :: big_points(:,:)
double precision, allocatable, private :: t(:,:)
double precision, allocatable, private :: ts(:,:)
double precision, allocatable, private :: big_t(:,:)
double precision, allocatable, private :: f(:,:)
double precision, allocatable, private :: fideal(:,:)
double precision, allocatable, private :: sdiff(:)
double precision, allocatable, private :: dist_diff(:)
double precision, allocatable, private :: big_sdiff(:)
double precision, allocatable, private :: s(:)
double precision, allocatable, private :: big_s(:)
double precision, allocatable, private :: f_dot_t(:)
double precision, allocatable, private :: f_dot_ts(:)
double precision, allocatable, private :: first_integral(:)
double precision, allocatable, private :: second_integral(:)
double precision, allocatable, private :: phi(:)
double precision, allocatable, private :: position_vec(:)
! FFTW arrays
real(C_DOUBLE), pointer, private :: phys(:)
complex(C_DOUBLE_COMPLEX), pointer, private :: coef(:)
type(C_PTR), private :: planr2c, planc2r, fftw_data
complex(kind=8), allocatable, private :: fftw_derivative(:)
complex(kind=8), allocatable, private :: fftw_integral(:)
contains

subroutine set_switches(sdiff_switch_input,t_switch_input,ts_switch_input,integral_switch_input)
implicit none
integer, intent(in) :: sdiff_switch_input
integer, intent(in) :: t_switch_input
integer, intent(in) :: ts_switch_input
integer, intent(in) :: integral_switch_input
sdiff_switch = sdiff_switch_input
t_switch = t_switch_input
ts_switch = ts_switch_input
integral_switch = integral_switch_input
big = 1
end subroutine set_switches

subroutine set_when_to_adapt(radius,npoints,distance_percentage)
implicit none
double precision, intent(in) :: distance_percentage
double precision, intent(in) :: radius
integer, intent(in) :: npoints
max_dist = 2*pi*radius/dble(npoints)*distance_percentage
write(*,*) max_dist
end subroutine set_when_to_adapt

subroutine allocate_arrays(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
integer :: i
allocate(points(ndim,npoints))
allocate(big_points(ndim,(1-big):(npoints+big)))
allocate(t(ndim,npoints))
allocate(ts(ndim,npoints))
allocate(big_t(ndim,(1-big):(npoints+big)))
allocate(f(ndim,npoints))
allocate(fideal(ndim,npoints))
allocate(sdiff(npoints))
allocate(dist_diff(npoints))
allocate(big_sdiff((1-big):(npoints+big)))
allocate(s(npoints))
allocate(big_s((1-big):(npoints+big)))
allocate(f_dot_t(npoints))
allocate(f_dot_ts(npoints))
allocate(first_integral(npoints))
allocate(second_integral(npoints))
allocate(phi(npoints))
allocate(position_vec(npoints))
! Allocate the fftw working arrays and set the plan
fftw_data = fftw_alloc_complex(int((npoints/2+1), C_SIZE_T))
call c_f_pointer(fftw_data, phys, [2*(npoints/2+1)])
call c_f_pointer(fftw_data, coef, [npoints/2+1])
planr2c = fftw_plan_dft_r2c_1d(npoints,phys,coef,FFTW_ESTIMATE)
planc2r= fftw_plan_dft_c2r_1d(npoints,coef,phys,FFTW_ESTIMATE)
! Set the derivative
allocate(fftw_derivative(npoints/2+1))
do concurrent (i=1:npoints/2+1)
fftw_derivative(i) = dcmplx(0.d0,dble(i-1))
end do
! Set the integral
allocate(fftw_integral(npoints/2+1))
fftw_integral(1) = dcmplx(0.d0,0.d0)
do concurrent (i=2:npoints/2+1)
fftw_integral(i) = dcmplx(0.d0,1.d0/dble(i-1))
end do
open(unit=100,file="output")
! Set up the position vector
do i=1,npoints
position_vec(i) = dble(i)/dble(npoints)
enddo 
end subroutine allocate_arrays

subroutine set_initial_points(eigvec1,eigvec2,eigval1,eigval2,fixed_point,radius,npoints)
implicit none
double precision, intent(inout) :: eigvec1(:)
double precision, intent(inout) :: eigvec2(:)
double precision, intent(in) :: eigval1
double precision, intent(in) :: eigval2
double precision, intent(inout) :: fixed_point(:)
double precision, intent(in) :: radius
integer, intent(in) :: npoints
integer :: i
do concurrent (i=1:npoints)
points(:,i) = fixed_point + radius * &
              (eigval1*dcos(dble(i-1)/dble(npoints)*2.d0*pi)*eigvec1 &
               + eigval2*dsin(dble(i-1)/dble(npoints)*2.d0*pi)*eigvec2)
!write(*,*) '**',i,points(:,i),i,'**'
end do
!stop
!points(:,1) =  (/ 1.99969539031d0, 0.03490481287d0, 0.d0 /)
end subroutine set_initial_points

subroutine points_to_tangent(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
if (t_switch .eq. 1) then
call points_to_tangent1(ndim,npoints)
endif
if (t_switch .eq. 2) then
call points_to_tangent2(ndim,npoints)
endif
end subroutine points_to_tangent

subroutine points_to_tangent1(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
integer :: i
do i=1,ndim
    phys = 0.d0
    phys(1:npoints) = points(i,:)
    call fftw_execute_dft_r2c(planr2c, phys, coef)
    coef = coef*fftw_derivative
    call fftw_execute_dft_c2r(planc2r, coef, phys)
    t(i,:) = phys(1:npoints)/dble(npoints)
end do
call normc(t,npoints)
!do concurrent (i=1:ndim)
!    write(*,*) i,'AAA',points(i,:),'AAA',t(i,:)
!enddo
end subroutine points_to_tangent1

subroutine points_to_tangent2(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
integer :: i
big_points(:,1:npoints) = points
big_points(:,(1-big):0) = points(:,(npoints-big+1):npoints)
big_points(:,(npoints+1):(npoints+big)) = points(:,1:(1+big-1))
big_sdiff(1:npoints) = sdiff
big_sdiff((1-big):0) = sdiff((npoints-big+1):npoints)
big_sdiff((npoints+1):(npoints+big)) = sdiff(1:(1+big-1))
do i=1,npoints
    t(:,i) = 0.5d0*(big_points(:,i)-big_points(:,i-1))/big_sdiff(i) + &
            0.5d0*(big_points(:,i+1)-big_points(:,i))/big_sdiff(i+1)
end do
call normc(t,npoints)
!do concurrent (i=1:ndim)
!    write(*,*) i,'AAA',points(i,:),'AAA',t(i,:)
!enddo
end subroutine points_to_tangent2

subroutine tangent_to_ts(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
if (ts_switch .eq. 1) then
call tangent_to_ts1(ndim,npoints)
endif
if (ts_switch .eq. 2) then
call tangent_to_ts2(ndim,npoints)
endif
end subroutine tangent_to_ts

subroutine tangent_to_ts1(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
integer :: i
do i=1,ndim
    phys = 0.d0
    phys(1:npoints) = t(i,:)
    call fftw_execute_dft_r2c(planr2c, phys, coef)
    coef = coef*fftw_derivative
    call fftw_execute_dft_c2r(planc2r, coef, phys)
    ts(i,:) = phys(1:npoints)/dble(npoints)
end do
!do concurrent (i=1:ndim)
!    write(*,*) i,'AAA',points(i,:),'AAA',ts(i,:)
!enddo
end subroutine tangent_to_ts1

subroutine tangent_to_ts2(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
integer :: i
big_t(:,1:npoints) = t
big_t(:,(1-big):0) = t(:,(npoints-big+1):npoints)
big_t(:,(npoints+1):(npoints+big)) = t(:,1:(1+big-1))
big_sdiff(1:npoints) = sdiff
big_sdiff((1-big):0) = sdiff((npoints-big+1):npoints)
big_sdiff((npoints+1):(npoints+big)) = sdiff(1:(1+big-1))
do i=1,npoints
    ts(:,i) = 0.5d0*(big_t(:,i)-big_t(:,i-1))/big_sdiff(i) + &
            0.5d0*(big_t(:,i+1)-big_t(:,i))/big_sdiff(i+1)
end do
!do concurrent (i=1:ndim)
!    write(*,*) i,'AAA',points(i,:),'AAA',t(i,:)
!enddo
end subroutine tangent_to_ts2

subroutine points_to_f(ndim,npoints)
use user_functions
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
integer :: i
integer :: iflag
iflag = 0
do i=1,npoints
call fcn ( ndim, points(:,i), f(:,i), iflag)
end do
call normc(f,npoints)
!do i=1,npoints
!write(*,*) i,f(:,i)
!end do
end subroutine points_to_f

subroutine normc(v,npoints)
implicit none
double precision, intent(inout) :: v(:,:)
integer, intent(in) :: npoints
integer :: i
do concurrent (i=1:npoints)
    v(:,i) = v(:,i)/(dsqrt(sum(v(:,i)*v(:,i))))
enddo
end subroutine normc

subroutine find_sdiff(npoints)
implicit none
integer, intent(in) :: npoints
if (sdiff_switch .eq. 1) then
call find_sdiff1(npoints)
endif
if (sdiff_switch .eq. 2) then
call find_sdiff2(npoints)
endif
end subroutine find_sdiff

subroutine find_sdiff1(npoints)
implicit none
integer, intent(in) :: npoints
integer :: i
sdiff(1) = 2.d0*pi*position_vec(1)
do concurrent(i=2:npoints)
sdiff(i) = 2.d0*pi*(position_vec(i) - position_vec(i-1))
enddo
!write(*,*) position_vec
!write(*,*) sdiff
!stop
end subroutine find_sdiff1

subroutine find_sdiff2(npoints)
implicit none
integer, intent(in) :: npoints
integer :: i
sdiff(1) = sqrt(sum((points(:,1)-points(:,npoints))*(points(:,1)-points(:,npoints))))
do concurrent(i=2:npoints)
sdiff(i) = sqrt(sum((points(:,i)-points(:,i-1))*(points(:,i)-points(:,i-1))))
enddo
sdiff = sdiff/sum(sdiff)*2.d0*pi
end subroutine find_sdiff2

subroutine find_distance(npoints)
implicit none
integer, intent(in) :: npoints
integer :: i
dist_diff(1) = sqrt(sum((points(:,1)-points(:,npoints))*(points(:,1)-points(:,npoints))))
do concurrent(i=2:npoints)
dist_diff(i) = sqrt(sum((points(:,i)-points(:,i-1))*(points(:,i)-points(:,i-1))))
enddo
end subroutine find_distance

subroutine check_points_far(npoints)
implicit none
integer, intent(in) :: npoints
integer :: i
call find_distance(npoints)
write(*,*) ANY(dist_diff .gt. max_dist)
end subroutine check_points_far

subroutine find_s(npoints,diff_vec,return_vec)
implicit none
integer, intent(in) :: npoints
double precision, intent(inout) :: diff_vec(1:npoints)
double precision, intent(inout) :: return_vec(1:npoints)
integer :: i
return_vec(1) = diff_vec(1)
do i=2,npoints
return_vec(i) = return_vec(i-1) + diff_vec(i)
enddo
!write(*,*) 'sdiff',diff_vec
!write(*,*) 's',return_vec
end subroutine find_s

subroutine dot(v1,v2,vr,npoints)
implicit none
double precision, intent(inout) :: v1(:,:)
double precision, intent(inout) :: v2(:,:)
double precision, intent(inout) :: vr(:)
integer, intent(in) :: npoints
integer :: i
do concurrent(i=1:npoints)
vr(i) = sum(v1(:,i)*v2(:,i))
!write(*,*) '**',i,vr(i),v1(:,i),v2(:,i),i,'***'
enddo
end subroutine dot

subroutine find_first_integral(npoints)
implicit none
integer, intent(in) :: npoints
if (integral_switch .eq. 1) then
call find_first_integral1(npoints)
endif
if (integral_switch .eq. 2) then
call find_first_integral2(npoints)
endif
end subroutine find_first_integral

subroutine find_first_integral1(npoints)
implicit none
integer, intent(in) :: npoints
phys = 0.d0
phys(1:npoints) = f_dot_ts
call fftw_execute_dft_r2c(planr2c, phys, coef)
! Integrate the mean
first_integral = dreal(coef(1))/12.d0*s
coef = coef*fftw_integral
call fftw_execute_dft_c2r(planc2r, coef, phys)
first_integral = first_integral + phys(1:npoints)/dble(npoints)
!write(*,*) 'f_dot_ts',f_dot_ts
!write(*,*) 'first_integral',first_integral
end subroutine find_first_integral1

subroutine find_first_integral2(npoints)
implicit none
integer, intent(in) :: npoints
integer :: i
first_integral(1) = sdiff(1)*(0.5d0*f_dot_ts(npoints) + 0.5d0*f_dot_ts(1))
do i=2,npoints
first_integral(i) = first_integral(i-1)+sdiff(i)*(0.5d0*f_dot_ts(i) + 0.5d0*f_dot_ts(i-1))
enddo
end subroutine find_first_integral2

subroutine find_second_integral(npoints)
implicit none
integer, intent(in) :: npoints
second_integral = s/(2.d0*pi)*first_integral(npoints)
!write(*,*) 'f_dot_ts',f_dot_ts
!write(*,*) 'second_integral',second_integral
end subroutine find_second_integral

subroutine find_fideal(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
integer :: i
call points_to_f(ndim,npoints)
call find_sdiff(npoints)
call find_s(npoints,sdiff,s)
call points_to_tangent(ndim,npoints)
call tangent_to_ts(ndim,npoints)
call dot(ts,f,f_dot_ts,npoints)
call find_first_integral(npoints)
call find_second_integral(npoints)
phi = first_integral - second_integral
phi = phi - sum(phi)/npoints
call dot(t,f,f_dot_t,npoints)
do i=1,npoints
    fideal(:,i) = f(:,i) + (phi(i) - f_dot_t(i)) * t(:,i)
!    write(*,*) '**',i,fideal(:,i),i,'**'
enddo
end subroutine find_fideal

subroutine timestep(dt)
implicit none
double precision, intent(in) :: dt
points = points + dt*fideal
end subroutine timestep

subroutine write_output(ringnum,npoints)
implicit none
integer, intent(in) :: ringnum
integer, intent(in) :: npoints
integer :: i
do i=1,npoints
write(100,*) dble(ringnum),points(:,i)
enddo
end subroutine write_output

subroutine deallocate_arrays()
implicit none
!integer :: i
!do concurrent (i=1:12)
!write(*,*) '**',i,points(:,i),i,'**'
!end do
deallocate(points)
deallocate(big_points)
deallocate(t)
deallocate(big_t)
deallocate(ts)
deallocate(f)
deallocate(fideal)
deallocate(sdiff)
deallocate(dist_diff)
deallocate(big_sdiff)
deallocate(s)
deallocate(big_s)
deallocate(f_dot_t)
deallocate(f_dot_ts)
deallocate(first_integral)
deallocate(second_integral)
deallocate(phi)
call fftw_destroy_plan(planr2c)
call fftw_destroy_plan(planc2r)
call fftw_free(fftw_data)
nullify(phys)
nullify(coef)
deallocate(fftw_derivative)
deallocate(fftw_integral)
deallocate(position_vec)
close(100)
end subroutine deallocate_arrays

end module ring
