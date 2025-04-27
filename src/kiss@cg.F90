submodule (kiss) cg

contains

    module subroutine cg_linop(A,b,x,rtol,atol,maxiter,callback,info)
        class(kiss_linop), intent(in) :: A
        real(wp), intent(in) :: b(:)
        real(wp), intent(inout) :: x(:)
        real(wp), intent(in), optional :: rtol, atol
        integer, intent(in), optional :: maxiter
        procedure(it_callback), optional :: callback
        integer, intent(out), optional :: info

        integer :: info_, maxiter_, n
#ifdef __INTEL_COMPILER
        type(cg_workspace(n=:)), allocatable :: it
#else
        type(cg_workspace) :: it
#endif

        n = size(x)

#ifdef __INTEL_COMPILER
        allocate ( cg_workspace(n=n) :: it)
#else
        it%n = n
        allocate(it%r(n), it%p(n), it%z(n), it%q(n))
#endif
        it%matvec => matvec_linop

        maxiter_ = 10*n

        if (present(rtol)) it%rtol = rtol
        if (present(atol)) it%atol = atol
        if (present(maxiter)) maxiter_ = maxiter

        call cg_wrk(it,b,x,maxiter_,info_,callback)

        if (present(info)) info = info_

    contains

        subroutine matvec_linop(n,x,y)
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            real(wp), intent(out) :: y(n)
            call A%apply(x,y)
        end subroutine

    end subroutine

    module subroutine cg_dense(A,b,x,rtol,atol,maxiter,callback,info)
        real(wp), intent(in) :: A(:,:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(inout) :: x(:)
        real(wp), intent(in), optional :: rtol, atol
        integer, intent(in), optional :: maxiter
        procedure(it_callback), optional :: callback
        integer, intent(out), optional :: info

        integer :: info_, maxiter_, n
#ifdef __INTEL_COMPILER
        type(cg_workspace(n=:)), allocatable :: it
#else
        type(cg_workspace) :: it
#endif

        n = size(x)

#ifdef __INTEL_COMPILER
        allocate ( cg_workspace(n=n) :: it)
#else
        it%n = n
        allocate(it%r(n), it%p(n), it%z(n), it%q(n))
#endif
        it%matvec => matvec_dense

        maxiter_ = 10*n

        if (present(rtol)) it%rtol = rtol
        if (present(atol)) it%atol = atol
        if (present(maxiter)) maxiter_ = maxiter

        call cg_wrk(it,b,x,maxiter_,info_,callback)

        if (present(info)) info = info_

    contains

#include "matvec_dense.fi"

    end subroutine


    subroutine cg_wrk(it,b,x,maxiter,info,callback)
#ifdef __INTEL_COMPILER
        type(cg_workspace(n=*)), intent(inout) :: it
#else
        type(cg_workspace), intent(inout) :: it
#endif
        real(wp), intent(in) :: b(it%n)
        real(wp), intent(inout) :: x(it%n)
        integer, intent(in) :: maxiter
        integer, intent(out) :: info
        procedure(it_callback), optional :: callback

        ! Local scalars
        real(wp) :: bnrm2, atol, rho_cur, rho_prev, alpha, beta, rv
        integer :: n, iter

        real(wp), parameter :: rhotol = epsilon(x)**2

        n = it%n
        bnrm2 = it%norm(n,b)

        if (bnrm2 == zero) then
            x = b
            info = 0
            return
        end if

        associate(r => it%r, p => it%p, z => it%z, q => it%q)

        atol = max(it%atol,it%rtol*bnrm2)
        call it%matvec(it%n,x,r)
        r = b - r

        do iter = 0, maxiter-1

            if (it%norm(it%n,r) < atol) then
                info = 0
                return
            end if

            call it%psolve(it%n,r,z)
            rho_cur = it%dotprod(it%n,r,z)

            if (iter > 0) then
                beta = rho_cur / rho_prev
                p = p*beta
                p = p + z
            else
                p = z
            end if

            call it%matvec(it%n,p,q)
            alpha = rho_cur / it%dotprod(it%n,p,q)
            x = x + alpha*p
            r = r - alpha*q
            rho_prev = rho_cur

            if (present(callback)) call callback(x)
        end do

        info = iter

        end associate

    end subroutine

end submodule
