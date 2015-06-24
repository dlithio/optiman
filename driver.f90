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
integer :: i
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
npoints_start = 64
ndim = 112
radius = 1.d-1
sdiff_switch = 2
t_switch = 2
ts_switch = 2
integral_switch = 2
dt = 1.d-3
steps_per_save = 300;
saved_rings = 50;
distance_percentagefar = 3.d0;
distance_percentageclose = 0.1d0;
f_switch = 1;
open(100,file="par")
do i=1,36
read(100,*) PAR(i)
enddo
close(100)
call setup(par)

npoints = npoints_start
call set_switches(sdiff_switch,t_switch,ts_switch,integral_switch,f_switch)
call allocate_arrays(ndim,npoints)

allocate(eigvec1(ndim))
allocate(eigvec2(ndim))
open(100,file="eigvec1")
do i=1,ndim
read(100,*) eigvec1(i)
enddo
close(100)
open(100,file="eigvec2")
do i=1,ndim
read(100,*) eigvec2(i)
enddo
close(100)
allocate(fixed_point(ndim))
open(100,file="85")
do i=1,ndim
read(100,*) fixed_point(i)
enddo
close(100)

eigval1 = 1.d0
eigval2 = 1.d0
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
    write(*,*) dble(ringnum)/dble(saved_rings)*100.d0,"%"
enddo

call deallocate_arrays()
close(217)

end program driver
