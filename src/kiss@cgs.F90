submodule (kiss) cgs

contains

    module subroutine cgs_dense(A,b,x,rtol,atol,maxiter,callback,info)
        real(wp), intent(in) :: A(:,:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(inout) :: x(:)
        real(wp), intent(in), optional :: rtol, atol
        integer, intent(in), optional :: maxiter
        procedure(it_callback), optional :: callback
        integer, intent(out), optional :: info

        integer :: info_, maxiter_, n
#ifdef __INTEL_COMPILER
        type(cgs_workspace(n=:)), allocatable :: it
#else
        type(cgs_workspace) :: it
#endif

        n = size(x)

#ifdef __INTEL_COMPILER
        allocate ( cgs_workspace(n=n) :: it)
#else
        it%n = n
        allocate(it%r(n),it%rtilde(n),it%q(n),it%p(n),it%phat(n),&
                 it%u(n),it%uhat(n),it%vhat(n))
#endif
        it%matvec => matvec_dense

        maxiter_ = 10*n

        if (present(rtol)) it%rtol = rtol
        if (present(atol)) it%atol = atol
        if (present(maxiter)) maxiter_ = maxiter

        call cgs_wrk(it,b,x,maxiter_,info_,callback)

        if (present(info)) info = info_

    contains

#include "matvec_dense.fi"

    end subroutine


    subroutine cgs_wrk(it,b,x,maxiter,info,callback)
#ifdef __INTEL_COMPILER
        type(cgs_workspace(n=*)), intent(inout) :: it
#else
        type(cgs_workspace), intent(inout) :: it
#endif
        real(wp), intent(in) :: b(it%n)
        real(wp), intent(inout) :: x(it%n)
        integer, intent(in) :: maxiter
        integer, intent(out) :: info
        procedure(it_callback), optional :: callback

        ! Local scalars
        real(wp) :: bnrm2, atol, omegatol
        real(wp) :: rho_cur, rho_prev, alpha, beta, rnorm, rv
        integer :: iter

        real(wp), parameter :: rhotol = epsilon(x)**2

        bnrm2 = it%norm(it%n,b)
        if (bnrm2 == zero) then
            x = b
            info = 0
            return
        end if

        atol = max(it%atol,it%rtol*bnrm2)

        associate(r=>it%r, rtilde=>it%rtilde,q=>it%q,p=>it%p, phat=>it%phat,&
                  u=>it%u, uhat=>it%uhat, vhat=>it%vhat)

        call it%matvec(it%n,x,r)
        r = b - r
        rtilde = r

        do iter = 0, maxiter - 1

            rnorm = it%norm(it%n,r)
            if (rnorm < atol) then
                info = 0
                return
            end if

            rho_cur = it%dotprod(it%n,rtilde,r)
            if (abs(rho_cur) < rhotol) then
                info = -10
                return
            end if

            if (iter > 0) then
                beta = rho_cur / rho_prev

                u = r + beta*q

                p = beta*p + q
                p = beta*p + u
            else
                p = r
                u = r
            end if

            call it%psolve(it%n,p,phat)
            call it%matvec(it%n,phat,vhat)
            rv = it%dotprod(it%n,rtilde,vhat)

            if (rv == zero) then
                info = -11
                return
            end if

            alpha = rho_cur / rv
            q = u - alpha*vhat
            u = u + q
            call it%psolve(it%n,u,uhat)
            x = x + alpha*uhat

            ! Due to numerical error build-up the actual residual is computed
            ! instead of the following two lines that were in the original
            ! FORTRAN templates, still using a single matvec.

            call it%matvec(it%n,x,r)
            r = b - r

            rho_prev = rho_cur

            if (present(callback)) call callback(x)

        end do

        info = iter

        end associate
    end subroutine

end submodule
