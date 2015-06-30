module ring
implicit none
double precision, private :: pi=4.d0*datan(1.d0)
double precision, private :: max_dist
double precision, private :: min_dist
integer, private :: sdiff_switch
integer, private :: t_switch
integer, private :: ts_switch
integer, private :: integral_switch
integer, private :: f_switch
integer, private :: interp_switch
integer :: big
double precision, allocatable, private :: points(:,:)
double precision, allocatable, private :: big_points(:,:)
double precision, allocatable, private :: t(:,:)
double precision, allocatable, private :: ts(:,:)
double precision, allocatable, private :: normct(:,:)
double precision, allocatable, private :: normcts(:,:)
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
double precision, allocatable, private :: points_new(:,:)
double precision, allocatable, private :: position_vec_new(:)
double precision, allocatable, private :: f_change_scalar(:)
double precision, allocatable, private :: f_change_scalart(:)
logical, private :: first_run = .TRUE.
contains

subroutine set_switches(sdiff_switch_input,t_switch_input,ts_switch_input,integral_switch_input,f_switch_input,interp_switch_input)
implicit none
integer, intent(in) :: sdiff_switch_input
integer, intent(in) :: t_switch_input
integer, intent(in) :: ts_switch_input
integer, intent(in) :: integral_switch_input
integer, intent(in) :: f_switch_input
integer, intent(in) :: interp_switch_input
sdiff_switch=sdiff_switch_input
t_switch=t_switch_input
ts_switch=ts_switch_input
integral_switch=integral_switch_input
f_switch=f_switch_input
interp_switch=interp_switch_input
big = 2
call system( 'rm output' )
call system( 'rm fdot' )
call system( 'rm t_angle' )
open(unit=217,file="output",access='stream')
open(unit=218,file="fdot",access='stream')
open(unit=219,file="t_angle",access='stream')
end subroutine set_switches

subroutine set_when_to_adapt(radius,npoints,distance_percentagefar,distance_percentageclose)
implicit none
double precision, intent(in) :: distance_percentagefar
double precision, intent(in) :: distance_percentageclose
double precision, intent(in) :: radius
integer, intent(in) :: npoints
max_dist = 2*pi*radius/dble(npoints)*distance_percentagefar
min_dist = 2*pi*radius/dble(npoints)*distance_percentageclose
!write(*,*) max_dist
end subroutine set_when_to_adapt

subroutine allocate_arrays(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
integer :: i
allocate(points(ndim,npoints))
if (first_run) then 
allocate(points_new(ndim,npoints))
endif
allocate(big_points(ndim,(1-big):(npoints+big)))
allocate(t(ndim,npoints))
allocate(ts(ndim,npoints))
allocate(normct(ndim,npoints))
allocate(normcts(ndim,npoints))
allocate(big_t(ndim,(1-big):(npoints+big)))
allocate(f(ndim,npoints))
allocate(fideal(ndim,npoints))
allocate(sdiff(npoints))
allocate(dist_diff(npoints))
allocate(big_sdiff((1-big):(npoints+big)))
allocate(s(npoints))
allocate(f_change_scalar(npoints))
allocate(f_change_scalart(npoints))
allocate(big_s((1-big):(npoints+big)))
allocate(f_dot_t(npoints))
allocate(f_dot_ts(npoints))
allocate(first_integral(npoints))
allocate(second_integral(npoints))
allocate(phi(npoints))
allocate(position_vec(npoints))
if (first_run) then 
allocate(position_vec_new(npoints))
endif
! Set up the position vector
do i=1,npoints
position_vec(i) = dble(i)/dble(npoints)
enddo 
first_run = .FALSE.
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
end subroutine set_initial_points

subroutine points_to_tangent(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
if (t_switch .eq. 2) then
call points_to_tangent2(ndim,npoints)
endif
if (t_switch .eq. 3) then
call points_to_tangent3(ndim,npoints)
endif
end subroutine points_to_tangent

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
end subroutine points_to_tangent2

subroutine points_to_tangent3(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
double precision :: cw(2,2*big+1)
integer :: i,ii
big_points(:,1:npoints) = points
big_points(:,(1-big):0) = points(:,(npoints-big+1):npoints)
big_points(:,(npoints+1):(npoints+big)) = points(:,1:(1+big-1))
big_sdiff(1:npoints) = sdiff
big_sdiff((1-big):0) = sdiff((npoints-big+1):npoints)
big_sdiff((npoints+1):(npoints+big)) = sdiff(1:(1+big-1))
call find_s(npoints+2*big,big_sdiff,big_s)
t = 0.d0
do i=1,npoints
    call weights(big_s(i),big_s((i-big):(i+big)),1,cw)
    do ii=-big,big
        t(:,i) = t(:,i) + big_points(:,i+ii)*cw(2,ii+big+1)
    enddo
end do
call normc(t,npoints)
!do concurrent (i=1:ndim)
!    write(*,*) i,'AAA',points(i,:),'AAA',t(i,:)
!enddo
end subroutine points_to_tangent3

subroutine weights(z,x,m,c)
! Calculates FD weights. The parameters are:
!  z   location where approximations are to be accurate,
!  x   vector with x-coordinates for grid points,
!  m   highest derivative that we want to find weights for
!  c   array size m+1,lentgh(x) containing (as output) in 
!      successive rows the weights for derivatives 0,1,...,m.
implicit none
double precision, intent(in) :: z
integer, intent(in) :: m
double precision, intent(in) :: x(2*big+1) 
double precision, intent(inout) :: c(m+1,2*big+1)
double precision :: c1,c2,c3,c4,c5
integer :: i,n,mn,j,k
n=2.d0*big+1
c = 0.d0
c1=1.d0
c4=x(1)-z
c(1,1)=1.d0
do i=2,n
   mn=min(i,m+1)
   c2=1.d0
   c5=c4
   c4=x(i)-z
   do j=1,(i-1)
      c3=x(i)-x(j)
      c2=c2*c3
      if (j .eq. (i-1)) then
         c(2:mn,i)=c1*((/(k, k=1,mn-1, 1)/)*c(1:mn-1,i-1)-c5*c(2:mn,i-1))/c2
         c(1,i)=-c1*c5*c(1,i-1)/c2
      end if
      c(2:mn,j)=(c4*c(2:mn,j)-(/(k, k=1,mn-1, 1)/)*c(1:mn-1,j))/c3
      c(1,j)=c4*c(1,j)/c3
    end do
    c1=c2
end do
end subroutine

subroutine weights2(z,x,m,c)
! Calculates FD weights. The parameters are:
!  z   location where approximations are to be accurate,
!  x   vector with x-coordinates for grid points,
!  m   highest derivative that we want to find weights for
!  c   array size m+1,lentgh(x) containing (as output) in 
!      successive rows the weights for derivatives 0,1,...,m.
implicit none
double precision, intent(in) :: z
integer, intent(in) :: m
double precision, intent(in) :: x(2*big) 
double precision, intent(inout) :: c(m+1,2*big)
double precision :: c1,c2,c3,c4,c5
integer :: i,n,mn,j,k
n=2.d0*big
c = 0.d0
c1=1.d0
c4=x(1)-z
c(1,1)=1.d0
do i=2,n
   mn=min(i,m+1)
   c2=1.d0
   c5=c4
   c4=x(i)-z
   do j=1,(i-1)
      c3=x(i)-x(j)
      c2=c2*c3
      if (j .eq. (i-1)) then
         c(2:mn,i)=c1*((/(k, k=1,mn-1, 1)/)*c(1:mn-1,i-1)-c5*c(2:mn,i-1))/c2
         c(1,i)=-c1*c5*c(1,i-1)/c2
      end if
      c(2:mn,j)=(c4*c(2:mn,j)-(/(k, k=1,mn-1, 1)/)*c(1:mn-1,j))/c3
      c(1,j)=c4*c(1,j)/c3
    end do
    c1=c2
end do
end subroutine

subroutine tangent_to_ts(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
if (ts_switch .eq. 2) then
call tangent_to_ts2(ndim,npoints)
endif
if (ts_switch .eq. 3) then
call tangent_to_ts3(ndim,npoints)
endif
end subroutine tangent_to_ts

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

subroutine tangent_to_ts3(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
double precision :: cw(2,2*big+1)
integer :: i,ii
big_t(:,1:npoints) = t
big_t(:,(1-big):0) = t(:,(npoints-big+1):npoints)
big_t(:,(npoints+1):(npoints+big)) = t(:,1:(1+big-1))
big_sdiff(1:npoints) = sdiff
big_sdiff((1-big):0) = sdiff((npoints-big+1):npoints)
big_sdiff((npoints+1):(npoints+big)) = sdiff(1:(1+big-1))
call find_s(npoints+2*big,big_sdiff,big_s)
ts = 0.d0
do i=1,npoints
    call weights(big_s(i),big_s((i-big):(i+big)),1,cw)
    do ii=-big,big
        ts(:,i) = ts(:,i) + big_t(:,i+ii)*cw(2,ii+big+1)
    enddo
end do
!do concurrent (i=1:ndim)
!    write(*,*) i,'AAA',points(i,:),'AAA',t(i,:)
!enddo
end subroutine tangent_to_ts3

subroutine points_to_f(ndim,npoints)
use user_functions
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
integer :: i
integer :: iflag
iflag = 0
do i=1,npoints
call fcn ( ndim, points(:,i), f(:,i))
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

subroutine find_distance(vin,vout,npoints)
implicit none
double precision, intent(inout) :: vin(:,:)
double precision, intent(inout) :: vout(:)
integer, intent(in) :: npoints
integer :: i
vout(1) = sqrt(sum((vin(:,1)-vin(:,npoints))*(vin(:,1)-vin(:,npoints))))
do concurrent(i=2:npoints)
vout(i) = sqrt(sum((vin(:,i)-vin(:,i-1))*(vin(:,i)-vin(:,i-1))))
enddo
end subroutine find_distance

subroutine check_points_far(npoints,something_wrong)
implicit none
integer, intent(in) :: npoints
logical, intent(inout) :: something_wrong
call find_distance(points_new,dist_diff,npoints)
something_wrong = any(dist_diff .gt. max_dist)
end subroutine check_points_far

subroutine check_points_close(npoints,something_wrong)
implicit none
integer, intent(in) :: npoints
logical, intent(inout) :: something_wrong
something_wrong = any(dist_diff .lt. min_dist)
end subroutine check_points_close

subroutine check_new_ring(npoints,something_wrong)
implicit none
integer, intent(in) :: npoints
logical, intent(inout) :: something_wrong
logical :: something_wrong_far
logical :: something_wrong_close
call check_points_far(npoints,something_wrong_far)
call check_points_close(npoints,something_wrong_close)
something_wrong = (something_wrong_far .OR. something_wrong_close)
end subroutine check_new_ring

subroutine accept_new_ring()
implicit none
points = points_new
end subroutine accept_new_ring

subroutine fix_old_ring(ndim,npointsold)
implicit none
integer, intent(in) :: ndim
integer, intent(inout) :: npointsold
integer :: new_npoints
integer :: old_point
integer :: new_point
logical :: made_change
made_change = .false.
if (any(dist_diff .gt. max_dist)) then
    made_change = .true.
    ! First interpolate where the points are too far away
    deallocate(points_new)
    deallocate(position_vec_new)
    new_npoints = npointsold + count(dist_diff .gt. max_dist)
    allocate(points_new(ndim,new_npoints))
    allocate(position_vec_new(new_npoints))
    new_point = 1
    old_point = 1
    big_points(:,1:npointsold) = points
    big_points(:,(1-big):0) = points(:,(npointsold-big+1):npointsold)
    big_points(:,(npointsold+1):(npointsold+big)) = points(:,1:(1+big-1))
    big_sdiff(1:npointsold) = sdiff
    big_sdiff((1-big):0) = sdiff((npointsold-big+1):npointsold)
    big_sdiff((npointsold+1):(npointsold+big)) = sdiff(1:(1+big-1))
    call find_s(npointsold+2*big,big_sdiff,big_s)
    do old_point=1,npointsold
        if (dist_diff(old_point) .gt. max_dist) then
            call interpolate(points,points_new,old_point,new_point,npointsold)
            call interpolatepos(position_vec,position_vec_new,old_point,new_point,npointsold)
            new_point = new_point + 1
        endif
        points_new(:,new_point) = points(:,old_point)
        position_vec_new(new_point) = position_vec(old_point)
        new_point = new_point + 1
    enddo
    call deallocate_arrays()
    call allocate_arrays(ndim,new_npoints)
    npointsold = new_npoints
    points = points_new
    position_vec = position_vec_new
    call find_distance(points,dist_diff,new_npoints)
endif
if (any(dist_diff .lt. min_dist) .and. (.not. made_change)) then
    ! Now remove points that are too close
    deallocate(points_new)
    deallocate(position_vec_new)
    new_npoints = npointsold - count(dist_diff .lt. min_dist)
    allocate(points_new(ndim,new_npoints))
    allocate(position_vec_new(new_npoints))
    new_point = 1
    do old_point=1,npointsold
        if (dist_diff(old_point) .ge. min_dist) then
            points_new(:,new_point) = points(:,old_point)
            position_vec_new(new_point) = position_vec(old_point)
            new_point = new_point + 1
        endif
    enddo
    call deallocate_arrays()
    call allocate_arrays(ndim,new_npoints)
    npointsold = new_npoints
    points = points_new
    position_vec = position_vec_new
endif
end subroutine fix_old_ring

subroutine interpolate1(oldarray,newarray,old_point,new_point,npointsold)
implicit none
double precision, intent(inout) :: oldarray(:,:)
double precision, intent(inout) :: newarray(:,:)
integer, intent(in) :: old_point
integer, intent(in) :: new_point
integer, intent(inout) :: npointsold
if (old_point .eq. 1) then
newarray(:,new_point) = (oldarray(:,1) + oldarray(:,npointsold))/2.d0
endif
if (old_point .ne. 1) then
newarray(:,new_point) = (oldarray(:,old_point) + oldarray(:,old_point-1))/2.d0
endif
end subroutine

subroutine interpolate2(newarray,old_point,new_point)
implicit none
double precision, intent(inout) :: newarray(:,:)
integer, intent(in) :: old_point
integer, intent(in) :: new_point
!double precision :: where_to_interp
double precision :: cw2(1,2*big)
integer :: ii
cw2 = 0.d0
!where_to_interp = (big_s(old_point) + big_s(old_point - 1))/2.d0
call weights2((big_s(old_point) + big_s(old_point - 1))/2.d0,big_s((old_point-big):(old_point+big-1)),0,cw2)
newarray(:,new_point) = 0.d0
do ii=-big,(big-1)
    newarray(:,new_point) = newarray(:,new_point) + big_points(:,old_point+ii)*cw2(1,ii+big+1)
enddo
end subroutine

subroutine interpolate(oldarray,newarray,old_point,new_point,npointsold)
implicit none
double precision, intent(inout) :: oldarray(:,:)
double precision, intent(inout) :: newarray(:,:)
integer, intent(in) :: old_point
integer, intent(in) :: new_point
integer, intent(inout) :: npointsold
if (interp_switch .eq. 1) then
call interpolate1(oldarray,newarray,old_point,new_point,npointsold)
endif
if (interp_switch .eq. 2) then
call interpolate2(newarray,old_point,new_point)
endif
end subroutine

subroutine interpolatepos(oldarray,newarray,old_point,new_point,npointsold)
implicit none
double precision, intent(inout) :: oldarray(:)
double precision, intent(inout) :: newarray(:)
integer, intent(in) :: old_point
integer, intent(in) :: new_point
integer, intent(inout) :: npointsold
if (old_point .eq. 1) then
newarray(new_point) = (oldarray(1))/2.d0
endif
if (old_point .ne. 1) then
newarray(new_point) = (oldarray(old_point) + oldarray(old_point-1))/2.d0
endif
end subroutine

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
integer :: i
first_integral(1) = sdiff(1)*(0.5d0*f_dot_ts(npoints) + 0.5d0*f_dot_ts(1))
do i=2,npoints
first_integral(i) = first_integral(i-1)+sdiff(i)*(0.5d0*f_dot_ts(i) + 0.5d0*f_dot_ts(i-1))
enddo
end subroutine find_first_integral

subroutine find_second_integral(npoints)
implicit none
integer, intent(in) :: npoints
second_integral = s/(2.d0*pi)*first_integral(npoints)
!write(*,*) 'f_dot_ts',f_dot_ts
!write(*,*) 'second_integral',second_integral
end subroutine find_second_integral

subroutine find_f(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
if (f_switch .eq. 1) then
call find_fideal(ndim,npoints)
endif
if (f_switch .eq. 2) then
call find_fzero(ndim,npoints)
endif
if (f_switch .eq. 3) then
call find_forig(ndim,npoints)
endif
end subroutine find_f

subroutine find_forig(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
integer :: i
call points_to_f(ndim,npoints)
call find_sdiff(npoints)
call find_s(npoints,sdiff,s)
fideal = f
end subroutine find_forig

subroutine find_fzero(ndim,npoints)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: npoints
integer :: i
call points_to_f(ndim,npoints)
call find_sdiff(npoints)
call find_s(npoints,sdiff,s)
call points_to_tangent(ndim,npoints)
call dot(t,f,f_dot_t,npoints)
do i=1,npoints
    fideal(:,i) = f(:,i) - f_dot_t(i) * t(:,i)
!    write(*,*) '**',i,fideal(:,i),i,'**'
enddo
call normc(fideal,npoints)
end subroutine find_fzero

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
!! My ad-hoc changes
!f_change_scalart = 0.d0
!do i=1,npoints
!    if ((dabs(f_dot_t(i)) .ge. 0.2d0) .and. (dabs(f_dot_t(i)) .le. 0.8d0)) then
!        !f_change_scalart(i) = (phi(i) - f_dot_t(i))*(0.5d0-dabs(f_dot_t(i)))/0.5d0
!        f_change_scalart(i) = (phi(i) - f_dot_t(i))
!    endif
!    if ((dabs(f_dot_t(i)) .ge. 0.1d0) .and. (dabs(f_dot_t(i)) .lt. 0.2d0)) then
!        !f_change_scalart(i) = (phi(i) - f_dot_t(i))*(0.5d0-dabs(f_dot_t(i)))/0.5d0
!        f_change_scalart(i) = (phi(i) - f_dot_t(i))*(0.2d0-dabs(f_dot_t(i)))/0.1d0
!    endif
!    if ((dabs(f_dot_t(i)) .gt. 0.8d0) .and. (dabs(f_dot_t(i)) .le. 0.9d0)) then
!        !f_change_scalart(i) = (phi(i) - f_dot_t(i))*(0.5d0-dabs(f_dot_t(i)))/0.5d0
!        f_change_scalart(i) = (phi(i) - f_dot_t(i))*(0.9d0-dabs(f_dot_t(i)))/0.1d0
!    endif
!!    if (dabs(f_dot_t(i)) .gt. 0.5d0) then
!!        f_change_scalart(i) = 0.d0
!!    endif
!enddo
!!f_change_scalar(1) = (sum(f_change_scalart(1:2))+f_change_scalart(npoints))/3.d0
!!f_change_scalar(npoints) = (sum(f_change_scalart(npoints-1:npoints))+f_change_scalart(1))/3.d0
!!do i=2,npoints-1
!!    f_change_scalar(i) = sum(f_change_scalart(i-1:i+1))/3.d0
!!enddo
!f_change_scalar = f_change_scalart
!! And now make the changes
!do i=1,npoints
!    fideal(:,i) = f(:,i) + f_change_scalar(i) * t(:,i)
!!    write(*,*) '**',i,fideal(:,i),i,'**'
!enddo
! My ad-hoc changes2
!f_change_scalart = 0.d0
!do i=1,npoints
!    if (dabs(f_dot_t(i)) .le. 0.9d0) then
!        fideal(:,i) = f(:,i) + (phi(i) - f_dot_t(i)) * t(:,i)
!    endif
!    if ((dabs(f_dot_t(i)) .gt. 0.9d0)) then
!        fideal(:,i) = 0.d0
!    endif
!enddo
!f_change_scalar(1) = (sum(f_change_scalart(1:2))+f_change_scalart(npoints))/3.d0
!f_change_scalar(npoints) = (sum(f_change_scalart(npoints-1:npoints))+f_change_scalart(1))/3.d0
!do i=2,npoints-1
!    f_change_scalar(i) = sum(f_change_scalart(i-1:i+1))/3.d0
!enddo
!f_change_scalar = f_change_scalart
!! And now make the changes
!do i=1,npoints
!    fideal(:,i) = f(:,i) + f_change_scalar(i) * t(:,i)
!!    write(*,*) '**',i,fideal(:,i),i,'**'
!enddo
! Another important change (?) with the adhoc changes.
! call normc(fideal,npoints)
normct(:,1:(npoints-1)) = t(:,2:npoints)
normct(:,npoints) = t(:,1)
do i=1,npoints
    !if ((sum(t(:,i)*normct(:,i))**0.5 .gt. 0.999d0) .and. (dabs(f_dot_t(i)) .le. 0.99d0)) then
        fideal(:,i) = f(:,i) + (phi(i) - f_dot_t(i)) * t(:,i)
    !endif
    !if (sum(t(:,i)*normct(:,i))**0.5 .le. 0.999d0) then
    !    fideal(:,i) = f(:,i)
    !endif
    !if ((dabs(f_dot_t(i)) .gt. 0.99d0)) then
    !    fideal(:,i) = 0.d0
    !endif
enddo
end subroutine find_fideal

subroutine timestep(dt)
implicit none
double precision, intent(in) :: dt
points_new = points + dt*fideal
end subroutine timestep

subroutine write_output(ringnum,npoints)
implicit none
integer, intent(in) :: ringnum
integer, intent(in) :: npoints
integer :: i
normct(:,1:(npoints-1)) = t(:,2:npoints)
normct(:,npoints) = t(:,1)
do i=1,npoints
write(217) dble(ringnum),points(:,i)
!write(218) dble(ringnum),f_dot_t(i)
!write(219) dble(ringnum),sum(t(:,i)*normct(:,i))**0.5
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
deallocate(normct)
deallocate(normcts)
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
deallocate(position_vec)
deallocate(f_change_scalar)
deallocate(f_change_scalart)
end subroutine deallocate_arrays

end module ring
