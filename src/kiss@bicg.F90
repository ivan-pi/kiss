submodule (kiss) bicg

contains

    module subroutine bicg_dense(A,b,x,rtol,atol,maxiter,callback,info)
            real(wp), intent(in) :: A(:,:)
            real(wp), intent(in) :: b(:)
            real(wp), intent(inout) :: x(:)
            real(wp), intent(in), optional :: rtol, atol
            integer, intent(in), optional :: maxiter
            procedure(it_callback), optional :: callback
            integer, intent(out), optional :: info

        integer :: info_, maxiter_, n
#ifdef __INTEL_COMPILER
        type(bicg_workspace(n=:)), allocatable :: it
#else
        type(bicg_workspace) :: it
#endif

        n = size(x)

#ifdef __INTEL_COMPILER
        allocate ( bicg_workspace(n=n) :: it)
#else
        it%n = n
        allocate(it%p(n), it%ptilde(n))
        allocate(it%r(n), it%rtilde(n))
        allocate(it%z(n), it%ztilde(n))
        allocate(it%q(n), it%qtilde(n))
#endif
        it%matvec => matvec_dense
        it%rmatvec => rmatvec_dense

        maxiter_ = 10*n

        if (present(rtol)) it%rtol = rtol
        if (present(atol)) it%atol = atol
        if (present(maxiter)) maxiter_ = maxiter

        call bicg_wrk(it,b,x,maxiter_,info_,callback)

        if (present(info)) info = info_

    contains

#include "matvec_dense.fi"
#include "rmatvec_dense.fi"

    end subroutine

    subroutine bicg_wrk(it,b,x,maxiter,info,callback)
#ifdef __INTEL_COMPILER
        type(bicg_workspace(n=*)), intent(inout) :: it
#else
        type(bicg_workspace), intent(inout) :: it
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
        atol = max(it%atol,it%rtol*bnrm2)

        if (bnrm2 == zero) then
            x = b
            info = 0
            return
        end if

        associate(p => it%p, ptilde => it%ptilde, &
                  r => it%r, rtilde => it%rtilde, &
                  z => it%z, ztilde => it%ztilde, &
                  q => it%q, qtilde => it%qtilde)

        call it%matvec(n,x,r)
        r = b - r
        rtilde = r

        do iter = 0, maxiter-1

            if (it%norm(n,r) < atol) then
                info = 0
                return
            end if

            call it%psolve(n,r,z)
            call it%rpsolve(n,rtilde,ztilde)
            rho_cur = it%dotprod(n,rtilde,z)

            if (abs(rho_cur) < rhotol) then
                info = -10
                return
            end if

            if (iter > 0) then
                beta = rho_cur / rho_prev
                p = p*beta
                p = p + z
                ptilde = ptilde * beta ! conj
                ptilde = ptilde + ztilde
            else
                p = z
                ptilde = ztilde
            end if

            call it%matvec(n,p,q)
            call it%rmatvec(n,ptilde,qtilde)
            rv = it%dotprod(n,ptilde,q)

            if (rv == zero) then
                info = -11
                return
            end if

            alpha = rho_cur / rv

            x = x + alpha*p
            r = r - alpha*q
            rtilde = rtilde - alpha * qtilde    ! conj
            rho_prev = rho_cur

            if (present(callback)) call callback(x)

        end do

        info = iter

        end associate

    end subroutine

end submodule
