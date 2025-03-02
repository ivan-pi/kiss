! kiss_mkl.F90 -- Procedural interface to oneMKL RCI ISS
!
! This module provides a procedural interface around the
! reverse-communication-based solvers provided in MKL. To use
! these, the caller has to implement a matvec subroutine, which
! is passed as a procedure pointer.
!
! Only double precision is available.
!
module kiss_mkl
implicit none
private

public :: dp
public :: rci_iss_cache
public :: matvec_fun
public :: pcg, new_cg_cache
public :: fgmres, new_fgmres_cache

integer, parameter :: dp = kind(1.0d0)

abstract interface
    subroutine matvec_fun(n,x,y)
        import dp
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: y(n)
    end subroutine
end interface

type :: rci_iss_cache
    integer :: n
    integer :: ipar(128)
    real(dp) :: dpar(128)
    real(dp), allocatable :: tmp(:)
    integer :: nrestart = 0
    procedure(matvec_fun), pointer, nopass :: matvec => null()
    procedure(matvec_fun), pointer, nopass :: psolve => null()
end type

include "mkl_rci.fi"

real(dp), parameter :: DEFAULT_RTOL = 1.0d-6
real(dp), parameter :: DEFAULT_ATOL = 0

contains

    function new_cg_cache(n,pc) result(cache)
        integer, intent(in) :: n
        logical, intent(in) :: pc
        type(rci_iss_cache) :: cache

        integer :: sz

        cache%n = n
        if (pc) then
            sz = n*4
        else
            sz = n*3
        end if
        allocate(cache%tmp(sz))

    end function

    subroutine pcg(cache,x,b,info,rtol,atol,maxiter,silent)
        type(rci_iss_cache), intent(inout) :: cache
        real(dp), intent(inout) :: x(cache%n)
        real(dp), intent(in) :: b(cache%n)
        integer, intent(out) :: info
        real(dp), intent(in), optional :: rtol, atol
        integer, intent(in), optional :: maxiter
        logical, intent(in), optional :: silent

        integer :: n, rci_request, itercount

        associate(n=>cache%n,ipar=>cache%ipar,dpar=>cache%dpar,tmp=>cache%tmp)

        call dcg_init(n,x,b,rci_request,ipar,dpar,tmp)
        if (rci_request /= 0) error stop "Failure in dcg_init"

        ! Silence warnings and errors
        if (present(silent)) then
            if (silent) ipar(6:7) = 0
        end if

        if (present(rtol)) dpar(1) = rtol
        if (present(atol)) dpar(2) = atol

        if (present(maxiter)) ipar(5) = maxiter

        ! Let MKL do the convergence test
        ipar(8) = 1
        ipar(9) = 1
        ipar(10) = 0

        ! Preconditioning
        if (associated(cache%psolve)) ipar(11) = 1

        ! Check zero norm of generated vector
        ipar(12) = 1

        call dcg_check(n,x,b,rci_request,ipar,dpar,tmp)
        if (rci_request /= 0) error stop "Failure in dcg_check"

        ! FIXME: if ipar(6) == 0 or ipar(7) == 0, we must handles
        !        errors ourselves

        revcom: do
            call dcg(n,x,b,rci_request,ipar,dpar,tmp)
            if (rci_request == 0) then
                info = 0
                return
            end if
            select case(rci_request)
            case(1)
                ! Matvec
                call cache%matvec(n,tmp(1),tmp(n+1))
                !call dcg_get(n,x,b,rci_request,ipar,dpar,tmp,itercount)
            case(2)
                ! Stopping test
                error stop "not supposed to be here (ipar(10) = 0)"
            case(3)
                ! Preconditioner
                call cache%psolve(n,tmp(2*n+1),tmp(3*n+1))
            case default
                error stop "in iterate loop, unknown request"
            end select
        end do revcom

        call dcg_get(n,x,b,rci_request,ipar,dpar,tmp,itercount)
        info = itercount

        end associate
    end subroutine pcg


    function new_fgmres_cache(n,nrestart) result(cache)
        integer, intent(in) :: n
        integer, intent(in), optional :: nrestart
        type(rci_iss_cache) :: cache

        integer :: sz

        cache%n = n

        if (present(nrestart)) then
            cache%nrestart = nrestart
        else
            cache%nrestart = min(150, n)
        end if

        associate(nrestart => cache%nrestart)
            sz = (2*nrestart+1)*n + nrestart*(nrestart+9)/2 + 1
        end associate

        allocate(cache%tmp(sz))

    end function

    subroutine fgmres(cache,x,b,info,rtol,atol,maxiter,silent)
        type(rci_iss_cache), intent(inout), target :: cache
        real(dp), intent(inout) :: x(cache%n)
        real(dp), intent(in) :: b(cache%n)
        integer, intent(out) :: info
        real(dp), intent(in), optional :: rtol, atol
        integer, intent(in), optional :: maxiter
        logical, intent(in), optional :: silent

        integer :: n, rci_request, itercount

        associate(n=>cache%n,ipar=>cache%ipar,dpar=>cache%dpar,tmp=>cache%tmp)

        call dfgmres_init(n,x,b,rci_request,ipar,dpar,tmp)
        if (rci_request /= 0) then
            info = rci_request
            return
        end if

        ! Silence warnings and errors
        if (present(silent)) then
            if (silent) ipar(6:7) = 0
        end if

        if (present(rtol)) dpar(1) = rtol
        if (present(atol)) dpar(2) = atol

        if (present(maxiter)) ipar(5) = maxiter

        ! Let MKL do the convergence test
        ipar(8) = 1
        ipar(9) = 1
        ipar(10) = 0

        ! Preconditioning
        if (associated(cache%psolve)) ipar(11) = 1

        ! Zero-norm check
        ipar(12) = 1

        ! Number of restarts allowed
        ipar(15) = cache%nrestart

        print *, "ipar(13) =", ipar(13)

        call dfgmres_check(n,x,b,rci_request,ipar,dpar,tmp)
         if (rci_request /= 0) then
            info = rci_request
            return
        end if

        revcom: do
            call dfgmres(n,x,b,rci_request,ipar,dpar,tmp)
            !print *, rci_request, dpar(5), dpar(7)
            if (rci_request == 0) then
                call dfgmres_get(n,x,b,rci_request,ipar,dpar,tmp,itercount)
                info = 0
                return
            end if
            select case(rci_request)
            case(1)
                ! multiply the vector with the matrix
                call cache%matvec(n,tmp(ipar(22)),tmp(ipar(23)))
            case(2)
                error stop "stopping test not available; disabled by ipar(10) = 0"
            case(3)
                ! apply the preconditioner
                call cache%psolve(n,tmp(ipar(22)),tmp(ipar(23)))
            case(4)
                ! check if the norm of the current orthogonal vector is zero
                error stop "zero norm test not available; disabled by ipar(12) = 0"
            case default
                error stop "fgmres: unknown rci_request"
            end select
        end do revcom

        call dfgmres_get(n,x,b,rci_request,ipar,dpar,tmp,itercount)
        info = itercount

        end associate
    end subroutine fgmres


    ! Here, B should contain the ILU-based preconditioner of the
    ! original matrix A
    ! function psolve_csr(n,x,y,b,ib,jb)
    !     integer, intent(in) :: n
    !     real(wp), intent(in) :: x(n)
    !     real(wp), intent(out) :: y(n)
    !     real(wp), intent(inout) :: trvec(n)
    !     real(wp), intent(in) :: a(*)
    !     integer, intent(in) :: ia(n+1), ja(*)

    !     call mkl_dcsrtrsv('L','N','U',n,b,ib,jb,x,trvec)
    !     call mkl_dcsrtrsv('U','N','N',n,b,ib,jb,trvec,y)

    ! end function

end module kiss_mkl
