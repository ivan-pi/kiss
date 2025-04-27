submodule (kiss) bicgstab

contains

    module subroutine bicgstab_linop(A,b,x,rtol,atol,maxiter,callback,info)
        class(kiss_linop), intent(in) :: A
        real(wp), intent(in) :: b(:)
        real(wp), intent(inout) :: x(:)
        real(wp), intent(in), optional :: rtol, atol
        integer, intent(in), optional :: maxiter
        procedure(it_callback), optional :: callback
        integer, intent(out), optional :: info

        integer :: info_, maxiter_, n
#ifdef __INTEL_COMPILER
        type(bicgstab_workspace(n=:)), allocatable :: it
#else
        type(bicgstab_workspace) :: it
#endif
        n = size(x)

#ifdef __INTEL_COMPILER
        allocate ( bicgstab_workspace(n=n) :: it)
#else
        it%n = n
        allocate(it%x0(n), it%r(n), it%p(n), it%v(n))
        allocate(it%rtilde(n), it%s(n), it%phat(n), it%shat(n), it%t(n))
#endif
        it%matvec => matvec

        maxiter_ = 10*n

        if (present(rtol)) it%rtol = rtol
        if (present(atol)) it%atol = atol
        if (present(maxiter)) maxiter_ = maxiter

        call bicgstab_wrk(it,b,x,maxiter_,info_,callback)

        if (present(info)) info = info_

    contains

        subroutine matvec(n,x,y)
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            real(wp), intent(out) :: y(n)
            call A%apply(x,y)
        end subroutine

    end subroutine

    module subroutine bicgstab_dense(A,b,x,rtol,atol,maxiter,callback,info)
        real(wp), intent(in) :: A(:,:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(inout) :: x(:)
        real(wp), intent(in), optional :: rtol, atol
        integer, intent(in), optional :: maxiter
        procedure(it_callback), optional :: callback
        integer, intent(out), optional :: info

        integer :: info_, maxiter_, n
#ifdef __INTEL_COMPILER
        type(bicgstab_workspace(n=:)), allocatable :: it
#else
        type(bicgstab_workspace) :: it
#endif

        n = size(x)

#ifdef __INTEL_COMPILER
        allocate ( bicgstab_workspace(n=n) :: it)
#else
        it%n = n
        allocate(it%x0(n), it%r(n), it%p(n), it%v(n))
        allocate(it%rtilde(n), it%s(n), it%phat(n), it%shat(n), it%t(n))
#endif
        it%matvec => matvec_dense

        maxiter_ = 10*n

        if (present(rtol)) it%rtol = rtol
        if (present(atol)) it%atol = atol
        if (present(maxiter)) maxiter_ = maxiter

        call bicgstab_wrk(it,b,x,maxiter_,info_,callback)

        if (present(info)) info = info_

    contains

#include "matvec_dense.fi"

    end subroutine


    subroutine bicgstab_wrk(it,b,x,maxiter,info,callback)
#ifdef __INTEL_COMPILER
        type(bicgstab_workspace(n=*)), intent(inout) :: it
#else
        type(bicgstab_workspace), intent(inout) :: it
#endif
        real(wp), intent(in) :: b(it%n)
        real(wp), intent(inout) :: x(it%n)
        integer, intent(in) :: maxiter
        integer, intent(out) :: info

        procedure(it_callback), optional :: callback

        real(wp) :: bnrm2, atol, rhotol, omegatol
        real(wp) :: beta, rho, rho_prev, alpha, omega, rv
        integer :: iter

#if USE_BLAS
        external :: dcopy, dscal, daxpy
#endif

        bnrm2 = it%norm(it%n,b)

        if (bnrm2 == zero) then
#if USE_BLAS
            call dcopy(it%n,b,1,x,1)
#else
            x = b
#endif
            info = 0
            return
        end if

        atol = max(it%atol,it%rtol*bnrm2)

        rhotol = epsilon(x)**2
        omegatol = rhotol

        associate(r => it%r, p => it%p, v => it%v, &
                  rtilde => it%rtilde, s => it%s, phat => it%phat, &
                  shat => it%shat, t => it%t)

            call it%matvec(it%n,x,r)
#if USE_BLAS
            call dscal(it%n,-one,r,1)
            call daxpy(it%n,one,b,1,r,1)
            call dcopy(it%n,r,1,rtilde,1)
#else
            r = b - r
            rtilde = r
#endif
            do iter = 0, maxiter-1

                if (it%norm(it%n,r) < it%atol) then
                    info = 0
                    return
                end if

                rho = it%dotprod(it%n,rtilde,r)
                if (abs(rho) < rhotol) then
                    info = -10                      ! rho breakdown
                    return
                end if

                if (iter > 0) then
                    if (abs(omega) < omegatol) then
                        info = -11                  ! omega breakdown
                        return
                    end if
                    beta = (rho / rho_prev) * (alpha / omega)
                    p = p - omega*v
                    p = p * beta
                    p = p + r
                else
                    ! First spin
#if USE_BLAS
                    call dcopy(it%n,r,1,p,1)
#else
                    p = r
#endif
                end if

            call it%psolve(it%n,p,phat)
            call it%matvec(it%n,phat,v)
            rv = it%dotprod(it%n,rtilde,v)

            if (rv == zero) then
                info = -11
                return
            end if

            alpha = rho / rv
#if USE_BLAS
            call daxpy(it%n,-alpha,v,1,r,1)
            call dcopy(it%n,r,1,s,1)
#else
            r = r - alpha*v
            s = r
#endif

            if (it%norm(it%n,s) < atol) then
                x = x + alpha*phat
                info = 0
                return
            end if

            call it%psolve(it%n,s,shat)
            call it%matvec(it%n,shat,t)
            omega = it%dotprod(it%n,t,s) / it%dotprod(it%n,t,t)

#if USE_BLAS
            call daxpy(it%n,alpha,phat,1,x,1)
            call daxpy(it%n,omega,shat,1,x,1)
            call daxpy(it%n,-omega,t,1,r,1)
#else
            x = x + alpha*phat
            x = x + omega*shat
            r = r - omega*t
#endif
            rho_prev = rho

            if (present(callback)) call callback(x)

        end do

        info = iter

    end associate

    end subroutine

end submodule