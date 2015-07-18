program driver
use ring
use user_functions
implicit none
integer :: ndim
integer :: npoints, npoints_start
integer :: sdiff_switch
integer :: f_switch
integer :: integral_switch
integer :: order_of_accuracy
integer :: i
integer :: fix_trys
double precision, allocatable :: eigvec1(:)
double precision, allocatable :: eigvec2(:)
double precision :: par(36)
complex*16 :: eigval1
complex*16 :: eigval2
double precision :: real_part,imag_part
double precision, allocatable :: fixed_point(:)
double precision :: radius
double precision :: dt
double precision :: distance_percentagefar
double precision :: distance_percentageclose
integer :: steps_per_save,stepnum
integer :: saved_rings,ringnum
integer :: initial_points_switch
integer :: tstep_switch
logical :: something_wrong
namelist /optiman_input/ npoints_start,radius,dt,steps_per_save,saved_rings,distance_percentagefar,distance_percentageclose,f_switch,order_of_accuracy,initial_points_switch,tstep_switch

open(100,file="optiman_input",delim='APOSTROPHE')
read(100,nml=optiman_input)
close(100)

sdiff_switch = 2
integral_switch = 2

open(100,file="par")
do i=1,36
read(100,*) PAR(i)
enddo
close(100)
ndim = int(par(36))
call setup(par)

npoints = npoints_start
call init_ring_mod(sdiff_switch,f_switch,order_of_accuracy,initial_points_switch,tstep_switch,ndim,npoints)

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
open(100,file="eigval1")
read(100,*) real_part
read(100,*) imag_part
eigval1 = complex(real_part,imag_part)
close(100)
open(100,file="eigval2")
read(100,*) real_part
read(100,*) imag_part
eigval2 = complex(real_part,imag_part)
close(100)
allocate(fixed_point(ndim))
open(100,file="fixed_point")
do i=1,ndim
read(100,*) fixed_point(i)
enddo
close(100)

call set_initial_points(eigvec1,eigvec2,eigval1,eigval2,fixed_point,radius,ndim,npoints)

!call find_f(ndim,npoints)
open(unit=101,file="header")
write(101,*) ndim
write(101,*) saved_rings
write(101,*) steps_per_save
write(101,*) dt
close(101)
!call deallocate_arrays()
!call allocate_arrays(ndim,npoints)
!call set_initial_points(eigvec1,eigvec2,eigval1,eigval2,fixed_point,radius,npoints)
ringnum=1
call set_when_to_adapt(radius,npoints,distance_percentagefar,distance_percentageclose)
call find_f(ndim,npoints)
call write_output(ringnum,npoints,ndim)

do ringnum=2,saved_rings
    do stepnum=1,steps_per_save
        call timestep(ndim,npoints,dt)
        call check_new_ring(npoints,something_wrong)
        fix_trys = 0
        do while (something_wrong)
            call fix_old_ring(ndim,npoints)
            call timestep(ndim,npoints,dt) 
            call check_new_ring(npoints,something_wrong)
            if (npoints .gt. 100000) then
                write(*,*) "too many points ",npoints
                stop
            endif
            fix_trys = fix_trys + 1
        enddo
        call accept_new_ring()
    enddo
    call write_output(ringnum,npoints,ndim)
    call progress(ringnum,saved_rings,npoints) 
enddo

call write_interp_info()
call deallocate_arrays()
close(217)
close(218)
close(219)

end program driver

subroutine progress(j,max_step,active_points)
  implicit none
  integer(kind=4)::j,k,max_step,num_stars,active_points
  double precision :: percent_done
  character(len=56)::bar="?????? active points, overall manifold ???% |          |"
!?????? active points, overall manifold ???% |          |
  percent_done = int(dble(j)/dble(max_step)*100.d0)
  num_stars = floor(percent_done/10.d0)
  write(unit=bar(1:6),fmt="(i6)") active_points
  write(unit=bar(40:42),fmt="(i3)") int(percent_done)
  do k=1, num_stars
    bar(45+k:45+k)="*"
  enddo
!  write(*,*) bar
  ! print the progress bar.
!  write(unit=6,fmt="(a1,a17)",advance="no") char(13), bar
write(unit=6,fmt="(a1,a56,$)") char(13), bar
  if (num_stars/=10) then
flush(unit=6)
  else
    write(unit=6,fmt=*)
  endif
  return
end subroutine progress
