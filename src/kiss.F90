module kiss_strat
    implicit none
    public

    integer, parameter :: wp = kind(1.0d0)
    private :: wp

contains

    function default_norm2(n,x) result(res)
        integer, intent(in) :: n
        real(wp), intent(in) :: x(n)
        real(wp) :: res
#if USE_BLAS
        real(kind(1.0d0)) :: dnrm2
        external :: dnrm2
        res = dnrm2(n,x,1)
#else
        res = norm2(x)
#endif
    end function

    function default_dotprod(n,x,y) result(res)
        integer, intent(in) :: n
        real(wp), intent(in) :: x(n), y(n)
        real(wp) :: res
#if USE_BLAS
        real(kind(1.0d0)) :: ddot
        external :: ddot
        res = ddot(n,x,1,y,1)
#else
        res = dot_product(x,y)
#endif
    end function

    ! FIXME: hacky way to include preconditioner
    subroutine eye(n,p,y)
        integer, intent(in) :: n
        real(wp), intent(in) :: p(n)
        real(wp), intent(out) :: y(n)
#if USE_BLAS
        external :: dcopy
        call dcopy(n,p,1,y,1)
#else
        y = p
#endif
    end subroutine

end module


module kiss

    ! The Intel Compiler appears to have a bug here.
    ! The problem is similar to the one reported here:
    !    https://community.intel.com/t5/Intel-Fortran-Compiler/Problem-with-procedure-pointers-using-both-ifx-and-ifort/m-p/1442308
    ! Renaming the routines appears to be a workaround
    use kiss_strat, times_eye =>  eye, &
        dnorm2 => default_norm2, dotprod => default_dotprod

    implicit none
    private

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: wp = dp

    public :: wp
    public :: it_callback

    public :: bicg_dense
    public :: bicgstab_dense
    public :: cg_dense
    public :: cgs_dense
    public :: qmr_dense
    public :: tfqmr_dense

    ! public callback
    abstract interface
        subroutine it_callback(x)
            import wp
            real(wp), intent(in) :: x(:)
        end subroutine
    end interface

    interface
        module subroutine bicgstab_dense(A,b,x,rtol,atol,maxiter,callback,info)
            real(wp), intent(in) :: A(:,:)
            real(wp), intent(in) :: b(:)
            real(wp), intent(inout) :: x(:)
            real(wp), intent(in), optional :: rtol, atol
            integer, intent(in), optional :: maxiter
            procedure(it_callback), optional :: callback
            integer, intent(out), optional :: info
        end subroutine
        module subroutine bicg_dense(A,b,x,rtol,atol,maxiter,callback,info)
            real(wp), intent(in) :: A(:,:)
            real(wp), intent(in) :: b(:)
            real(wp), intent(inout) :: x(:)
            real(wp), intent(in), optional :: rtol, atol
            integer, intent(in), optional :: maxiter
            procedure(it_callback), optional :: callback
            integer, intent(out), optional :: info
        end subroutine
        module subroutine cg_dense(A,b,x,rtol,atol,maxiter,callback,info)
            real(wp), intent(in) :: A(:,:)
            real(wp), intent(in) :: b(:)
            real(wp), intent(inout) :: x(:)
            real(wp), intent(in), optional :: rtol, atol
            integer, intent(in), optional :: maxiter
            procedure(it_callback), optional :: callback
            integer, intent(out), optional :: info
        end subroutine
        module subroutine cgs_dense(A,b,x,rtol,atol,maxiter,callback,info)
            real(wp), intent(in) :: A(:,:)
            real(wp), intent(in) :: b(:)
            real(wp), intent(inout) :: x(:)
            real(wp), intent(in), optional :: rtol, atol
            integer, intent(in), optional :: maxiter
            procedure(it_callback), optional :: callback
            integer, intent(out), optional :: info
        end subroutine
        module subroutine qmr_dense(A,b,x,rtol,atol,maxiter,callback,info)
            real(wp), intent(in) :: A(:,:)
            real(wp), intent(in) :: b(:)
            real(wp), intent(inout) :: x(:)
            real(wp), intent(in), optional :: rtol, atol
            integer, intent(in), optional :: maxiter
            procedure(it_callback), optional :: callback
            integer, intent(out), optional :: info
        end subroutine
        module subroutine tfqmr_dense(A,b,x,rtol,atol,maxiter,callback,info)
            real(wp), intent(in) :: A(:,:)
            real(wp), intent(in) :: b(:)
            real(wp), intent(inout) :: x(:)
            real(wp), intent(in), optional :: rtol, atol
            integer, intent(in), optional :: maxiter
            procedure(it_callback), optional :: callback
            integer, intent(out), optional :: info
        end subroutine
    end interface


    abstract interface
        function norm_fun(n,x) result(res)
            import wp
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            real(wp) :: res
        end function
        function dotprod_fun(n,x,y) result(res)
            import wp
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n), y(n)
            real(wp) :: res
        end function
        subroutine matvec_fun(n,x,y)
            import wp
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            real(wp), intent(out) :: y(n)
        end subroutine
    end interface


    real(wp), parameter :: DEFAULT_RTOL = 1.0e-5_wp

#ifdef __INTEL_COMPILER
    type :: bicgstab_workspace(n)
        integer, len :: n
        real(wp) :: x0(n), r(n), p(n), v(n)
        real(wp) :: rtilde(n), s(n), phat(n), shat(n), t(n)
#else
    type :: bicgstab_workspace
        integer :: n
        real(wp), allocatable :: x0(:), r(:), p(:), v(:)
        real(wp), allocatable :: rtilde(:), s(:), phat(:), shat(:), t(:)
#endif
        real(wp) :: rtol = DEFAULT_RTOL
        real(wp) :: atol = 0

        procedure(matvec_fun), pointer, nopass :: matvec => null() ! A * x

        ! Strategies, might be useful to replace with custom implementations
        procedure(norm_fun), pointer, nopass :: norm => dnorm2
        procedure(dotprod_fun), pointer, nopass :: dotprod => dotprod
        procedure(matvec_fun), pointer, nopass :: psolve => times_eye ! Left preconditioner

        procedure(matvec_fun), pointer, nopass :: rmatvec => null() ! A^T * x
        procedure(matvec_fun), pointer, nopass :: rpsolve => times_eye ! Right preconditioner
    end type


#ifdef __INTEL_COMPILER
    type :: bicg_workspace(n)
        integer, len :: n
        real(wp), dimension(n) :: p, ptilde
        real(wp), dimension(n) :: r, rtilde
        real(wp), dimension(n) :: z, ztilde
        real(wp), dimension(n) :: q, qtilde
#else
    type :: bicg_workspace
        integer :: n
        real(wp), allocatable :: p(:), ptilde(:)
        real(wp), allocatable :: r(:), rtilde(:)
        real(wp), allocatable :: z(:), ztilde(:)
        real(wp), allocatable :: q(:), qtilde(:)
#endif
        real(wp) :: rtol = DEFAULT_RTOL
        real(wp) :: atol = 0

        procedure(matvec_fun), pointer, nopass :: matvec => null() ! A * x
        procedure(matvec_fun), pointer, nopass :: rmatvec => null() ! A^T * x

        ! Strategies, might be useful to replace with custom implementations
        procedure(norm_fun), pointer, nopass :: norm => dnorm2
        procedure(dotprod_fun), pointer, nopass :: dotprod => dotprod

        procedure(matvec_fun), pointer, nopass :: psolve => times_eye ! Left preconditioner
        procedure(matvec_fun), pointer, nopass :: rpsolve => times_eye ! Right preconditioner
    end type


#ifdef __INTEL_COMPILER
    type :: cg_workspace(n)
        integer, len :: n
        real(wp), dimension(n) :: r, p, z, q
#else
    type :: cg_workspace
        integer :: n
        real(wp), allocatable, dimension(:) :: r, p, z, q
#endif
        real(wp) :: rtol = DEFAULT_RTOL
        real(wp) :: atol = 0

        procedure(matvec_fun), pointer, nopass :: matvec => null() ! A * x

        ! Strategies, might be useful to replace with custom implementations
        procedure(norm_fun), pointer, nopass :: norm => dnorm2
        procedure(dotprod_fun), pointer, nopass :: dotprod => dotprod

        procedure(matvec_fun), pointer, nopass :: psolve => times_eye ! Left preconditioner
    end type


#ifdef __INTEL_COMPILER
    type :: cgs_workspace(n)
        integer, len :: n
        real(wp), dimension(n) :: r, rtilde, q
        real(wp), dimension(n) :: p, phat, u, uhat, vhat
#else
    type :: cgs_workspace
        integer :: n
        real(wp), allocatable, dimension(:) :: r, rtilde, q
        real(wp), allocatable, dimension(:) :: p, phat, u, uhat, vhat
#endif
        real(wp) :: rtol = DEFAULT_RTOL
        real(wp) :: atol = 0

        procedure(matvec_fun), pointer, nopass :: matvec => null() ! A * x

        ! Strategies, might be useful to replace with custom implementations
        procedure(norm_fun), pointer, nopass :: norm => dnorm2
        procedure(dotprod_fun), pointer, nopass :: dotprod => dotprod

        procedure(matvec_fun), pointer, nopass :: psolve => times_eye ! Left preconditioner
    end type

    type :: linop
        procedure(matvec_fun), pointer, nopass :: matvec => times_eye
        procedure(matvec_fun), pointer, nopass :: rmatvec => times_eye
    end type


#ifdef __INTEL_COMPILER
    type :: qmr_workspace(n)
        integer, len :: n
        real(wp), dimension(n) :: r, q, d, s
        real(wp), dimension(n) :: y, ytilde
        real(wp), dimension(n) :: v, vtilde
        real(wp), dimension(n) :: w, wtilde
        real(wp), dimension(n) :: z, ztilde
        real(wp), dimension(n) :: p, ptilde
#else
    type :: qmr_workspace
        integer :: n
        real(wp), allocatable, dimension(:) :: r, q, d, s
        real(wp), allocatable, dimension(:) :: y, ytilde
        real(wp), allocatable, dimension(:) :: v, vtilde
        real(wp), allocatable, dimension(:) :: w, wtilde
        real(wp), allocatable, dimension(:) :: z, ztilde
        real(wp), allocatable, dimension(:) :: p, ptilde
#endif
        real(wp) :: rtol = DEFAULT_RTOL
        real(wp) :: atol = 0

        procedure(matvec_fun), pointer, nopass :: matvec => null() ! A * x
        procedure(matvec_fun), pointer, nopass :: rmatvec => null() ! A^T * x

        ! Strategies, might be useful to replace with custom implementations
        procedure(norm_fun), pointer, nopass :: norm => dnorm2
        procedure(dotprod_fun), pointer, nopass :: dotprod => dotprod

        type(linop) :: m1, m2
    end type


#ifdef __INTEL_COMPILER
    type :: tfqmr_workspace(n)
        integer, len :: n
        real(wp), dimension(n) :: r, u, w, rstar, v, uhat, unext, d, z
#else
    type :: tfqmr_workspace
        integer :: n
        real(wp), allocatable, dimension(:) :: r, u, w, rstar, v, uhat, unext, d, z
#endif
        real(wp) :: rtol = DEFAULT_RTOL
        real(wp) :: atol = 0

        procedure(matvec_fun), pointer, nopass :: matvec => null() ! A * x

        ! Strategies, might be useful to replace with custom implementations
        procedure(norm_fun), pointer, nopass :: norm => dnorm2
        procedure(dotprod_fun), pointer, nopass :: dotprod => dotprod

        procedure(matvec_fun), pointer, nopass :: psolve => times_eye ! Left preconditioner
    end type


    ! Scalar constants
    real(wp), parameter :: zero = 0
    real(wp), parameter :: one  = 1

end module
