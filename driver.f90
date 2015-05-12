program driver
use ring
implicit none
integer :: ndim
integer :: npoints, npoints_start
integer :: sdiff_switch
integer :: t_switch
double precision, allocatable :: eigvec1(:)
double precision, allocatable :: eigvec2(:)
double precision :: eigval1
double precision :: eigval2
double precision, allocatable :: fixed_point(:)
double precision :: radius

npoints_start = 12
ndim = 3
radius = 2.d0
sdiff_switch = 1
t_switch = 1

npoints = npoints_start
call set_switches(sdiff_switch,t_switch)
call allocate_arrays(ndim,npoints)

allocate(eigvec1(ndim))
allocate(eigvec2(ndim))
eigvec1 = (/ 1.d0, 0.d0, 0.d0 /)
eigvec2 = (/ 0.d0, 1.d0, 0.d0 /)
eigval1 = 1.d0
eigval2 = 1.d0
allocate(fixed_point(ndim))
fixed_point = (/ 0.d0, 0.d0, 0.d0 /)
call set_initial_points(eigvec1,eigvec2,eigval1,eigval2,fixed_point,radius,npoints)

call points_to_f(ndim,npoints)
call find_sdiff(npoints)
call points_to_tangent(ndim,npoints)
call tangent_to_ts(ndim,npoints)


call deallocate_arrays()



end program driver
