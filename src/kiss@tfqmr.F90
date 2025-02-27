submodule (kiss) tfqmr

contains

    module subroutine tfqmr_dense(A,b,x,rtol,atol,maxiter,callback,info)
        real(wp), intent(in) :: A(:,:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(inout) :: x(:)
        real(wp), intent(in), optional :: rtol, atol
        integer, intent(in), optional :: maxiter
        procedure(it_callback), optional :: callback
        integer, intent(out), optional :: info

        integer :: info_, maxiter_, n
#ifdef __INTEL_COMPILER
        type(tfqmr_workspace(n=:)), allocatable :: it
#else
        type(tfqmr_workspace) :: it
#endif

        n = size(x)

#ifdef __INTEL_COMPILER
        allocate ( tfqmr_workspace(n=n) :: it)
#else
        it%n = n
        allocate(it%r(n),it%u(n),it%w(n),it%rstar(n),it%v(n),it%uhat(n),&
                 it%unext(n),it%d(n),it%z(n))
#endif
        it%matvec => matvec_dense

        maxiter_ = min(10000, n * 10)

        if (present(rtol)) it%rtol = rtol
        if (present(atol)) it%atol = atol
        if (present(maxiter)) maxiter_ = maxiter

        call tfqmr_wrk(it,b,x,maxiter_,info_,callback)

        if (present(info)) info = info_

    contains

        subroutine matvec_dense(n,x,y)
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            real(wp), intent(out) :: y(n)
#if USE_BLAS
            external :: dgemv
            real(kind(1.0d0)), parameter :: alpha = 1.0d0, beta = 0.0d0
            call dgemv('N',n,n,alpha,A,size(A,1),x,1,beta,y,1)
#else
            y = matmul(A, x)
#endif
        end subroutine

    end subroutine


    subroutine tfqmr_wrk(it,b,x,maxiter,info,callback)
#ifdef __INTEL_COMPILER
        type(tfqmr_workspace(n=*)), intent(inout) :: it
#else
        type(tfqmr_workspace), intent(inout) :: it
#endif
        real(wp), intent(in) :: b(it%n)
        real(wp), intent(inout) :: x(it%n)
        integer, intent(in) :: maxiter
        integer, intent(out) :: info
        procedure(it_callback), optional :: callback

        ! Local scalars
        real(wp) :: bnrm2, atol
        real(wp) :: rho, rho_last, alpha, beta, rnorm, rv
        real(wp) :: c, d, theta, eta, r0norm, tau, rstar, vtrstar
        integer :: n, iter
        logical :: even

        real(wp), allocatable :: utmp(:)

        real(wp), parameter :: rhotol = epsilon(x)**2

        bnrm2 = it%norm(it%n,b)
        if (bnrm2 == zero) then
            x = b
            info = 0
            return
        end if

        n = it%n

        atol = max(it%atol,it%rtol*bnrm2)
        allocate(utmp(n))

        associate(r=>it%r,u=>it%u,w=>it%w,rstar=>it%rstar,v=>it%v,&
                  uhat=>it%uhat,unext=>it%unext,d=>it%d,z=>it%z)

        call it%matvec(it%n,x,r)
        r = b - r

        u = r
        w = r
        rstar = r

        call it%matvec(n,r,uhat)
        call it%psolve(n,uhat,v)
        uhat = v

        d     = zero
        theta = zero
        eta   = zero

        rho = sum(rstar*r)
        rho_last = rho
        r0norm = sqrt(rho)
        tau = r0norm

        if (r0norm == zero) then
            info = 0
            return
        end if

        ! Iterative solver loop
        do iter = 0, maxiter-1

            even = mod(iter, 2) == 0

            if (even) then
                vtrstar = it%dotprod(n,rstar,v)
                if (vtrstar == zero) then
                    info = -1
                    return
                end if
                alpha = rho / vtrstar
                unext = u - alpha * v               ! [1]-(5.6)
            end if

            w = w - alpha * uhat                    ! [1]-(5.8)
            d = u + (theta**2 / alpha) * eta * d    ! [1]-(5.5)
            theta = it%norm(n,w) / tau              ! [1]-(5.2)
            c = sqrt(one/(one + theta**2))
            tau = tau * (theta * c)
            eta = (c**2) * alpha                    ! [1]-(5.4)
            call it%psolve(n, d, z)
            x = x + eta * z

            if (present(callback)) call callback(x)

            ! Convergence criterion
            if (tau * sqrt(real(iter+1, wp)) < atol) then
                info = 0
                return
            end if

            if (.not. even) then
                rho = it%dotprod(n,rstar,w)  ! [1]-(5.7)
                beta = rho / rho_last
                u = w + beta * u
                v = beta * uhat + (beta**2) * v
                call it%matvec(n,u,utmp)
                call it%psolve(n,utmp,uhat)
                v = v + uhat
            else
                call it%matvec(n,unext,utmp)
                call it%psolve(n,utmp,uhat)
                u = unext
                rho_last = rho
            end if
        end do

        info = iter

        end associate
    end subroutine

end submodule
