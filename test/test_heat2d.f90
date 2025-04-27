program test_heat2d

  use kiss_bmp, only: kiss_print_bmp

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: Nx = 500, Ny = 500
  integer, parameter :: Nt = 100
  integer, parameter :: N = Nx * Ny
  real(dp), parameter :: Lx = 1.0_dp, Ly = 1.0_dp, alpha = 1.0_dp, dt = 0.001_dp
  real(dp), parameter :: tol = 1.0e-8_dp
  integer, parameter :: maxiter = 1000
  real(dp) :: dx, dy
  real(dp) :: u(N), u_next(N), rhs(N)

  dx = Lx / (Nx - 1)
  dy = Ly / (Ny - 1)

  call initialize_field(u,nx,nx,dx,dy)

  call kiss_print_bmp( &
    reshape(real(u),shape=[nx,ny]), &
    minv=0.0, &
    maxv=1.0, &
    filename="heat0.bmp")

contains

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
        yval = dy * (j - 1) / (dy * (Ny - 1))
        u2d(1,j) = yval
      end do

      ! Bottom boundary (y = 0): varies along x from 0 to 1
      do i = 1, Nx
        xval = dx * (i - 1) / (dx * (Nx - 1))
        u2d(i,1) = xval
      end do

      ! Other boundaries (top and right) remain 0

    end subroutine

end program
