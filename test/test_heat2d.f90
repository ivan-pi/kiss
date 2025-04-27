
module heat2d
  use kiss, only: kiss_linop
  implicit none
  private

  public :: apply_laplacian
  public :: stencil_op

  integer, parameter :: dp = kind(1.0d0)

  type, extends(kiss_linop) :: stencil_op
    integer :: Nx, Ny
    real(dp) :: dx, dy, coeff
  contains
    procedure :: apply => apply_stencil
  end type

  type, extends(kiss_linop) :: jacobi_pc_op
    integer :: Nx, Ny
    real(dp) :: dx, dy, coeff
  contains
    procedure :: apply => apply_jacobi
  end type

contains

  subroutine apply_stencil(op,x,y)
      class(stencil_op), intent(in) :: op
      real(dp), intent(in) :: x(op%n)
      real(dp), intent(out) :: y(op%n)

      if (op%N /= op%Nx*op%Nx) error stop "shape mismatch"
      call apply_laplacian(y,x,op%coeff,op%dx,op%dy,op%Nx,op%Ny)

  end subroutine

  ! (I + coef * Laplacian) * u
  subroutine apply_laplacian(Ru, u, coeff, dx, dy, Nx, Ny)
    real(dp), intent(inout) :: Ru(Nx,Ny)
    real(dp), intent(in)  :: u(Nx,Ny)
    real(dp), intent(in)  :: coeff, dx, dy
    integer, intent(in)   :: Nx, Ny
    integer :: i, j
    real(dp) :: rdx2, rdy2

    rdx2 = 1.0_dp / dx**2
    rdy2 = 1.0_dp / dy**2

    ! copy BC's (not efficient)
    Ru = u

    do j = 2, Ny-1
      do i = 2, Nx-1
        Ru(i,j) = u(i,j) + coeff * ( &
        rdx2 * (u(i+1,j) - 2.0_dp*u(i,j) + u(i-1,j)) + &
        rdy2 * (u(i,j+1) - 2.0_dp*u(i,j) + u(i,j-1)) )
      end do
    end do

  end subroutine

  ! Jacobi preconditioner for the (I + coef *Laplacian) operator
  subroutine apply_jacobi(op,x,y)
      class(jacobi_pc_op), intent(in) :: op
      real(dp), intent(in) :: x(op%n)
      real(dp), intent(out) :: y(op%n)

      real(dp) :: invm, m
      integer :: i, j

      m = 1 - op%coeff * (-2.0_dp/op%dx**2 - 2.0_dp/op%dy**2)
      invm = 1.0_dp/m

      do j = 2, op%Ny-1
        do i = 2, op%Nx-1
          y((j-1)*op%Nx + i) = invm * x((j-1)*op%Nx + i)
        end do
      end do

  end subroutine


  ! Reaction-diffusion stencil in a hollow cylindrical catalyst
  !
  ! This problem is taken from a former university assignment
  ! under supervision of Prof. Janez Levec.
  !
  subroutine apply_crd(Ru,u,r,dr,dy,k,De,Nr,Ny)
    integer, intent(in) :: Nr, Ny
    real(dp), intent(in) :: r(Nr), dr, dy, k, De
    real(dp), intent(out) :: Ru(Nr,Ny)
    real(dp), intent(in) :: u(Nr,Ny)

    real(dp) :: rdr2, rdy2
    integer :: i, j

    rdr2 = 1.0_dp/dr**2
    rdy2 = 1.0_dp/dy**2

    do j = 2, Ny-1
      do i = 2, Nr-1
        Ru(i,j) = &
          rdr2*(u(i+1,j) - 2.0_dp*u(i,i) + u(i-1,j)) + &
          0.5_dp*dr*(u(i+1,j) - u(i-1,2))/r(i) + &
          rdy2*(u(i,j+1) - 2.0_dp*u(i,i) + u(i,j+1)) - &
          - (k/De) * u(i,j)
      end do
    end do

  end subroutine

end module

program test_heat2d

  use kiss_bmp, only: kiss_print_bmp
  use heat2d, only: stencil_op, apply_laplacian
  use kiss, only: cg_linop

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: Nx = 100, Ny = 100
  integer, parameter :: Nt = 100
  integer, parameter :: N = Nx * Ny
  real(dp), parameter :: Lx = 1.0_dp, Ly = 1.0_dp
  real(dp), parameter :: alpha = 0.01_dp, dt = 0.1_dp
  real(dp), parameter :: tol = 1.0e-8_dp
  integer, parameter :: maxiter = 1000

  real(dp) :: dx, dy, Lcoef, Rcoef
  real(dp) :: u(N), u_next(N), rhs(N)
  integer :: step, num_steps, info

  type(stencil_op) :: euler, cn_R, cn_L

  dx = Lx / (Nx - 1)
  dy = Ly / (Ny - 1)

  euler = stencil_op(N,Nx,Ny,dx,dy,coeff=alpha*dt)
  cn_R = stencil_op(N,Nx,Ny,dx,dy,coeff=alpha*dt*0.5)
  cn_L = stencil_op(N,Nx,Ny,dx,dy,coeff=-alpha*dt*0.5)

  call initialize_field(u,nx,nx,dx,dy)

  call kiss_print_bmp( &
    reshape(real(u),shape=[nx,ny]), &
    minv=0.0, &
    maxv=1.0, &
    filename="heat0.bmp")

  block
    real(dp) :: x(Nx), y(Ny)
    integer :: i, unit
    x = [((i-1)*dx,i=1,Nx)]
    y = [((i-1)*dy,i=1,Ny)]
    open(newunit=unit,file="heat0.txt")
    call print_nonuniform_matrix(unit,x,y,&
        reshape(u,shape=[Nx,Ny]))
    close(unit)
  end block

  ! Copy BC's
  u_next = u

  Rcoef = alpha * dt

  num_steps = 100
  do step = 0, num_steps

!    call apply_laplacian(&
!      u_next,u,Rcoef,dx,dy,Nx,Ny)
!    call euler%apply(u,u_next)

    call cn_R%apply(u,u_next)
    u = u_next

    print *, "starting cg"
    call cg_linop(A=cn_L,b=u_next,x=u,&
            callback=print_residual_norm,info=info)
    if (info /= 0) error stop "cg did not converge"
    print *, "cg done"
    !u = u_next

    if (mod(step,1000) == 0) then
      print *, step
      print *, maxval(u_next), minval(u_next)
      print *, maxval(u), maxval(u_next)
    end if

    !u = u_next

  end do

  call kiss_print_bmp( &
    reshape(real(u),shape=[nx,ny]), &
    minv=0.0, &
    maxv=1.0, &
    filename="heat_end.bmp")

  block
    real(dp) :: x(Nx), y(Ny)
    integer :: i, unit
    x = [((i-1)*dx,i=1,Nx)]
    y = [((i-1)*dy,i=1,Ny)]
    open(newunit=unit,file="heat_end.txt")
    call print_nonuniform_matrix(unit,x,y,&
        reshape(u,shape=[Nx,Ny]))
    close(unit)
  end block

contains

    subroutine print_residual_norm(x)
        real(dp), intent(in) :: x(:)
        real(dp) :: ax(size(x))
        call cn_L%apply(x,ax)
        print *, "  residual norm = ", norm2(u_next - ax)
    end subroutine


  subroutine initialize_field(u, Nx, Ny, dx, dy)
    real(dp), intent(out), target :: u(Nx*Ny)
    integer, intent(in) :: Nx, Ny
    real(dp), intent(in) :: dx, dy

    real(dp), pointer, contiguous :: u2d(:,:)
    integer :: i, j
    real(dp) :: xval, yval

    u2d(1:Nx,1:Ny) => u
    u2d = 0.0_dp

    ! Left boundary (x = 0): varies along y from 0 to 1
    do j = 1, Ny
      u2d(1,j) = Ny - j
    end do
    u2d(1,:) = u2d(1,:) / (Ny - 1)

    ! Bottom boundary (y = 0): varies along x from 0 to 1
    do i = 1, Nx
      u2d(i,1) = Nx - i
    end do
    u2d(:,1) = u2d(:,1) / (Nx - 1)

    ! Other boundaries (top and right) remain 0

  end subroutine

  ! ASCII output format suitable for gnuplot
  !
  ! Plot the resulting data with using
  !
  !    plot "data" nonuniform matrix with image
  !
  subroutine print_nonuniform_matrix(unit,x,y,z,fmt)
  integer, intent(in) :: unit
  real(dp), intent(in) :: x(:), y(:)
  real(dp), intent(in) :: z(:,:)      ! shape(z) := (NY,NX)
  character(len=*), intent(in), optional :: fmt

  character(len=32) :: cx, cy
  character(len=:), allocatable :: fmt1, fmt2, dp_fmt
  integer :: row

  if (size(z,1) /= size(y) .or. &
  size(z,2) /= size(x)) then
  error stop "[plot_nonuniform_matrix] error: incompatible array dimensions"
  end if

  write(cx,'(I0)') size(x)
  write(cy,'(I0)') size(y)

  dp_fmt = 'G0'
  if (present(fmt)) then
  dp_fmt = trim(fmt)
  end if

  fmt1 = '(I0,'//trim(cx)//'(1X,'//dp_FMT//'))'
  !fmt1 = "(I0,4(1X,G0))"

  fmt2 = '('//dp_FMT//','//trim(cx)//'(1X,'//dp_FMT//'))'

  print *, fmt1
  print *, fmt2

  write(unit,fmt1) size(x)+1, x
  do row = 1, size(y)
  write(unit,fmt2) y(row),z(row,:)
  end do

  ! Add empty line so we can continue plotting if needed
  write(unit,*)

  close(unit)

  end subroutine


  ! Legacy VTK format for structured grid
  !
  ! Routine is borrowed from the repo of Neil Carlson
  ! under MIT License:
  !
  !    https://github.com/nncarlson/index-map
  !
  subroutine vtk_plot(filename)
    character(len=*), intent(in) :: filename
    integer :: iunit
    open(newunit=iunit,file=trim(filename),action="write")
    write(iunit,'("# vtk DataFile Version 3.0")')
    write(iunit,'("test_heat2d")')
    write(iunit,'("ASCII")')
    write(iunit,'("DATASET STRUCTURED_POINTS")')
    write(iunit,'("DIMENSIONS",3(1x,i0))') Nx+1, Ny+1, 1
    write(iunit,'("ORIGIN 0 0 0")')
    write(iunit,'("SPACING",3(1x,g0))') dx, dy, 1
    write(iunit,'("POINT_DATA ",i0)') (Nx*Ny)**2
    write(iunit,'("SCALARS u float 1")')
    write(iunit,'("LOOKUP_TABLE default")')
    write(iunit,'(G0)') u
    close(iunit)
  end subroutine

end program
