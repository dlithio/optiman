! Reads in from the following files:
! 1) fixed_point_input
! 2) par (note that par(36) should be ndim. Also, par should be 36 long)
! 3) initial_guess

program fixed_point
use user_functions
implicit none
double precision, allocatable :: a(:,:),b(:,:),f(:),x(:)
double precision, allocatable :: wr(:),wi(:),vl(:,:),vr(:,:),work(:)
character :: trans,JOBVL,JOBVR
integer :: M,N,NRHS,LDA,LDB,INFO
integer :: LDVL,LDVR,LWORK
integer, allocatable :: IPIV(:)
integer :: i,user1,user2
integer :: max_iterations,current_iteration
double precision :: tolerance,par(36)
external :: DGETRF,DGETRS,DGEEV
namelist /fixed_point_input/ tolerance,max_iterations

open(100,file="fixed_point_input",delim='APOSTROPHE')
read(100,nml=fixed_point_input)
close(100)
open(100,file="par")
do i=1,36
read(100,*) PAR(i)
enddo
close(100)
N = par(36)
call setup(par)
TRANS = 'N'
M = N
NRHS = 1
LDA = N
LDB = N
allocate(a(LDA,N))
allocate(b(LDB,NRHS))
allocate(IPIV(N))
allocate(x(N))
allocate(f(N))

open(100,file="initial_guess")
do i=1,N
read(100,*) x(i)
enddo
close(100)
call fcn ( n, x, f)

current_iteration = 0
do while (sum(f*f)**0.5 .gt. tolerance)
    b(:,1) = f
    call get_jac(n,x,a)
    call DGETRF( M, N, a, LDA, IPIV, INFO )
    call DGETRS( TRANS, N, NRHS, a, LDA, IPIV, b, LDB, INFO )
    x = x - b(:,1)
    call fcn ( n, x, f)
    current_iteration = current_iteration + 1
    if (current_iteration .gt. max_iterations) then
        write(*,*) "there have been ",current_iteration," newton iterations"
        write(*,*) "the current x is ",x
        write(*,*) "the current f is ",f
        write(*,*) "the current fnorm is ",sum(f*f)**0.5
        write(*,*) "Please adjust parameters and run again to get convergence"
        stop
    endif
enddo

write(*,*) "In ",current_iteration," newton iterations we found a fixed point at"
write(*,*) ""
write(*,*) x
write(*,*) ""
write(*,*) "the norm of the vector field at this point is ",sum(f*f)**0.5
write(*,*) "it's has been saved in the file fixed_point for use by optiman."
write(*,*) ""
open(100,file="fixed_point")
do i=1,n
write(100,*) x(i)
enddo
close(100)

! Now we find the eigenvectors and values
call get_jac(n,x,a)
JOBVL = 'N'
JOBVR = 'V'
LDVL = n
LDVR = n
lwork = 4*n
allocate(wr(n))
allocate(wi(n))
allocate(VL(LDVL,N))
allocate(VR(LDVR,N))
allocate(work(lwork))
call DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

write(*,*) "Here we list out all the eigenvalues of the jacobian at that fixed point"
write(*,*) "that have positive real parts. If you do not see at least 2 positive eigenvalues,"
write(*,*) "or see the ones you want, you should manually reverse the vector field simply"
write(*,*) "by making the output of fcn negative"
write(*,*) "If you select complex eigenvalues, you must select them in their conjugate pairs"

do i=1,n
    if (wr(i) .gt. 0.d0) then
        write(*,*)
        write(*,*) "i ",i
        write(*,*) "real part of eigenvalue ",wr(i)
        write(*,*) "imaginary part of eigenvalue ",wi(i)
        ! write(*,*) "vr ",vr(:,i)
    end if
enddo

write(*,*) ""
write(*,*) "Select the value of the first eigenvalue/vector you would like"
write(*,*) "to use"
read(*,*) user1
write(*,*) "Select the value of the first eigenvalue/vector you would like"
write(*,*) "to use"
read(*,*) user2
write(*,*) "your selections have been saved to eigval1,eigval2,eigvec1,eigvec2"
write(*,*) "for use by optiman"
open(100,file="eigval1")
write(100,*) wr(user1)
write(100,*) wi(user1)
close(100)
open(100,file="eigval2")
write(100,*) wr(user2)
write(100,*) wi(user2)
close(100)
open(100,file="eigvec1")
do i=1,n
write(100,*) vr(i,user1)
enddo
close(100)
open(100,file="eigvec2")
do i=1,n
write(100,*) vr(i,user2)
enddo
close(100)

deallocate(a)
deallocate(b)
deallocate(IPIV)
deallocate(x)
deallocate(f)
deallocate(wr)
deallocate(wi)
deallocate(VL)
deallocate(VR)
deallocate(work)

end program
