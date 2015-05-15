program driver
use ring
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
integer :: numchecks

npoints_start = 128
ndim = 3
radius = 15.d0
sdiff_switch = 2
t_switch = 2
ts_switch = 2
integral_switch = 2
dt = 1.d-2
steps_per_save = 1000;
saved_rings = 15;
distance_percentagefar = 3.d0;
distance_percentageclose = 0.25d0;
f_switch = 3;

npoints = npoints_start
call set_switches(sdiff_switch,t_switch,ts_switch,integral_switch,f_switch)
call allocate_arrays(ndim,npoints)

allocate(eigvec1(ndim))
allocate(eigvec2(ndim))
eigvec1 = (/ 1.d0, 0.d0, 0.d0 /)
eigvec2 = (/ 0.d0, 1.d0, 0.d0 /)
eigvec1 = (/ 0.61482d0, -0.78867d0, 0.d0 /)
eigvec2 = (/ 0.d0, 0.d0, 1.d0 /)
eigval1 = 1.d0
eigval2 = 1.d0
allocate(fixed_point(ndim))
fixed_point = (/ 0.d0, 0.d0, 0.d0 /)
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
        numchecks = 0
        do while (something_wrong)
            call fix_old_ring(ndim,npoints)
            call find_f(ndim,npoints)
            call timestep(dt) 
            call check_new_ring(npoints,something_wrong)
            numchecks = numchecks + 1
            if (numchecks .gt. 1000) then
                write(*,*) 'hi'
            endif
        enddo
        call accept_new_ring()
    enddo
    call write_output(ringnum,npoints)
enddo

call deallocate_arrays()

end program driver
