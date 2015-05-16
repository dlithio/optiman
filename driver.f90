program driver
use ring
use user_functions
implicit none
integer :: ndim
integer :: npoints, npoints_start
integer :: sdiff_switch
integer :: t_switch
integer :: f_switch
integer :: ts_switch
integer :: integral_switch
double precision, allocatable :: eigvec1(:)
double precision, allocatable :: eigvec2(:)
double precision, allocatable :: par(:)
double precision :: eigval1
double precision :: eigval2
double precision, allocatable :: fixed_point(:)
double precision :: radius
double precision :: dt
double precision :: distance_percentagefar
double precision :: distance_percentageclose
integer :: steps_per_save,stepnum
integer :: saved_rings,ringnum
logical :: something_wrong

allocate(par(36))
npoints_start = 256
ndim = 8
radius = 1.d-1
sdiff_switch = 2
t_switch = 2
ts_switch = 2
integral_switch = 2
dt = 1.d-3
steps_per_save = 2500;
saved_rings = 40;
distance_percentagefar = 3.d0;
distance_percentageclose = 0.50d-30;
f_switch = 3;
par(1) = 33.0d0
par(2) = dble(ndim)
call setup(par)

npoints = npoints_start
call set_switches(sdiff_switch,t_switch,ts_switch,integral_switch,f_switch)
call allocate_arrays(ndim,npoints)

allocate(eigvec1(ndim))
allocate(eigvec2(ndim))
!eigvec1 = (/ 1.d0, 0.d0, 0.d0 /)
!eigvec2 = (/ 0.d0, 1.d0, 0.d0 /)
!eigvec1 = (/ 0.61482d0, -0.78867d0, 0.d0 /)
!eigvec2 = (/ 0.d0, 0.d0, 1.d0 /)
eigvec1 = 0.d0
eigvec1(1) = .2659d0
eigvec1(3) = -.88239d0
eigvec1(5) = .2651d0
eigvec1(7) = -.04816d0
eigvec2 = 0.d0
eigvec2(1) = .27920d0
eigvec2(3) = 0.0d0
eigvec2(5) = .01184d0
eigvec2(7) = -.00151d0
eigval1 = 1.d0
eigval2 = 1.d0
allocate(fixed_point(ndim))
fixed_point = 0.d0
fixed_point = (/ 0.d0,   5.7274d00,   0.d0,&
        -1.965d00,  0.d0,   2.7422d-01,&
         0.d0,  -3.2382d-02 /)

call set_initial_points(eigvec1,eigvec2,eigval1,eigval2,fixed_point,radius,npoints)

call find_f(ndim,npoints)
open(unit=101,file="header")
write(101,*) ndim
write(101,*) saved_rings
write(101,*) steps_per_save
write(101,*) dt
close(101)
call deallocate_arrays()
call allocate_arrays(ndim,npoints)
call set_initial_points(eigvec1,eigvec2,eigval1,eigval2,fixed_point,radius,npoints)
ringnum=1
call set_when_to_adapt(radius,npoints,distance_percentagefar,distance_percentageclose)
call write_output(ringnum,npoints)

do ringnum=2,saved_rings
    do stepnum=1,steps_per_save
        call find_f(ndim,npoints)
        call timestep(dt)
        call check_new_ring(npoints,something_wrong)
        do while (something_wrong)
            call fix_old_ring(ndim,npoints)
            call find_f(ndim,npoints)
            call timestep(dt) 
            call check_new_ring(npoints,something_wrong)
        enddo
        call accept_new_ring()
    enddo
    call write_output(ringnum,npoints)
enddo

call deallocate_arrays()

end program driver
