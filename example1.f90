module user_functions
use projector_mod_phys
use auto_mod_physf
use nse_mod_physf
use utility_mod
use status_mod
implicit none
integer, parameter, private :: NDIM=112
double precision, private :: modpar(36)
!double precision, allocatable, private :: uu(:)
!double precision, allocatable, private :: b(:)
!double precision, private :: alpha
contains

subroutine setup (par)

implicit none
double precision, intent(inout) :: par(:)
!alpha = par(1)
!allocate(uu(3*int(par(2))))
!allocate(b(int(par(2))))
integer m, n, NDIM_check
double precision U(NDIM), T, F(NDIM)
integer, parameter :: WP = selected_real_kind(15,307)
complex(WP), allocatable :: steady_state(:,:)
complex(WP), allocatable :: forcing(:,:)
integer :: step_num
integer :: nmodes
integer :: err_num
integer :: kx1,ky1,kx2,ky2
real(WP) :: par1,ar1,ai1,ar2,ai2
logical :: my_status
namelist /old_par/ par1,step_num,nmodes,kx1,ky1,ar1,ai1,kx2,ky2,ar2,ai2

modpar = par
open(100,file="old_par",delim='APOSTROPHE')
read(100,nml=old_par)
close(100)
PAR = 0.d0
err_num = 0
PAR(1) = par1
call projector_mod_phys_init(step_num,m,n,NDIM_check,kx1,ky1)
if (NDIM .ne. NDIM_check) then
    write(*,*) "ndim mismatch"
    stop
endif
call utility_mod_init()
call auto_mod_physf_init(m,NDIM,par1,nmodes,&
                            kx1,ky1,ar1,ai1,kx2,ky2,&
                            ar2,ai2)
allocate(steady_state(0:m/2,-m/2:m/2))
allocate(forcing(0:m/2,-m/2:m/2))
! Get the steady state, ignore the set forcing
call get_external_force_and_steady_state(forcing,steady_state,par)
call project_to_1d(steady_state,U)
! Reinit the nse_mod with that new force
call nse_mod_physf_init(m,kx1,ky1,NDIM,forcing)
F = 0.d0
call save_energy_of_force(forcing)
deallocate(steady_state)
deallocate(forcing)
call auto_state_to_force(U,F,PAR)
my_status = .TRUE.
call set_status(my_status)
call system( 'rm output' )
open(unit=217,file="output",access='stream')

end subroutine setup

subroutine fcn ( n, x, fvec, iflag)
implicit none
integer, intent(in) :: n
double precision, intent(inout) :: x(n) ! really intent(in), want to avoid a mem copy though
double precision, intent(inout) :: fvec(n)
integer, intent(inout) :: iflag
integer :: j
iflag = 0
! A test manifold
!fvec(1) = x(1)
!fvec(2) = x(2)
!fvec(3) = -x(3)

! The lorenz manifold
!fvec(1)=10.d0*(x(2)-x(1));
!fvec(2)=28.d0*x(1)-x(2)-x(1) * x(3);
!fvec(3)=x(1)*x(2) - 8.d0/3.d0 * x(3);
!fvec = -1.d0*fvec

!! The kse
!call getBsingle(n,x,b)
!do j=1,n
!    fvec(j)=dfloat(-4*j**4)*x(j)+alpha*dfloat(j**2)*x(j) - b(j)
!enddo
call auto_state_to_force(x,fvec,modpar)
end subroutine fcn

!subroutine getBsingle(n,u,b_out)
!! Computes the bilinear term of the KSE
!! Does not include the extra high modes.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!IMPLICIT none
!integer, intent(in) :: n
!double precision, intent(inout) :: u(1:)
!double precision, intent(inout) :: b_out(1:)
!DOUBLE PRECISION :: mysum
!integer :: k,kjsgn,j
!uu = 0.d0
!uu(1:n) = u
!b_out = 0.d0
!do  k=1,n
!    mysum=0.d0
!    do  j=1,n
!        if (k .gt. j) then
!            kjsgn=1
!            mysum=mysum+j*u(j)*(uu(j+k)+kjsgn*uu(abs(k-j)))
!        else if (k .lt. j) then
!            kjsgn=-1
!            mysum=mysum+j*u(j)*(uu(j+k)+kjsgn*uu(abs(k-j)))
!        else
!            kjsgn=0
!            mysum=mysum+j*u(j)*uu(j+k)
!        endif
!    enddo
!    b_out(k)=alpha*mysum/2
!enddo
!end subroutine getBsingle

end module
