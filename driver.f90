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
double precision :: dt
integer :: steps_per_save,stepnum
integer :: saved_rings,ringnum

npoints_start = 12
ndim = 3
radius = 2.d0
sdiff_switch = 1
t_switch = 1
dt = 1.d-2
steps_per_save = 1;
saved_rings = 100;

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

call find_fideal(ndim,npoints)
call deallocate_arrays()
call allocate_arrays(ndim,npoints)
call set_initial_points(eigvec1,eigvec2,eigval1,eigval2,fixed_point,radius,npoints)

do ringnum=2,saved_rings
    do stepnum=1,steps_per_save
        call find_fideal(ndim,npoints)
        call timestep(dt)
    enddo
enddo


call deallocate_arrays()



end program driver
