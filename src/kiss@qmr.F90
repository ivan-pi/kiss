submodule (kiss) qmr

contains

    module subroutine qmr_dense(A,b,x,rtol,atol,maxiter,callback,info)
        real(wp), intent(in) :: A(:,:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(inout) :: x(:)
        real(wp), intent(in), optional :: rtol, atol
        integer, intent(in), optional :: maxiter
        procedure(it_callback), optional :: callback
        integer, intent(out), optional :: info

        integer :: info_, maxiter_, n
#ifdef __INTEL_COMPILER
        type(qmr_workspace(n=:)), allocatable :: it
#else
        type(qmr_workspace) :: it
#endif

        n = size(x)
#ifdef __INTEL_COMPILER
        allocate ( qmr_workspace(n=n) :: it)
#else
        it%n = n
        allocate(it%r(n), it%q(n), it%d(n), it%s(n))
        allocate(it%y(n), it%ytilde(n))
        allocate(it%v(n), it%vtilde(n))
        allocate(it%w(n), it%wtilde(n))
        allocate(it%z(n), it%ztilde(n))
        allocate(it%p(n), it%ptilde(n))
#endif
        it%matvec => matvec_dense
        it%rmatvec => rmatvec_dense

        ! it%M1 = ...
        ! it%M2 = ...

        maxiter_ = 10*n

        if (present(rtol)) it%rtol = rtol
        if (present(atol)) it%atol = atol
        if (present(maxiter)) maxiter_ = maxiter

        call qmr_wrk(it,b,x,maxiter_,info_,callback)

        if (present(info)) info = info_

    contains

        subroutine matvec_dense(n,x,y)
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            real(wp), intent(out) :: y(n)
#if USE_BLAS
            external :: dgemv
            real(dp), parameter :: alpha = 1.0_dp, beta = 0.0_dp
            call dgemv('N',n,n,alpha,A,size(A,1),x,1,beta,y,1)
#else
            y = matmul(A, x)
#endif
        end subroutine

        subroutine rmatvec_dense(n,x,y)
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            real(wp), intent(out) :: y(n)
#if USE_BLAS
            external :: dgemv
            real(dp), parameter :: alpha = 1.0_dp, beta = 0.0_dp
            call dgemv('T',n,n,alpha,A,size(A,1),x,1,beta,y,1)
#else
            y = matmul(transpose(A), x)
#endif
        end subroutine
    end subroutine


    subroutine qmr_wrk(it,b,x,maxiter,info,callback)
#ifdef __INTEL_COMPILER
        type(qmr_workspace(n=*)), intent(inout) :: it
#else
        type(qmr_workspace), intent(inout) :: it
#endif
        real(wp), intent(in) :: b(it%n)
        real(wp), intent(inout) :: x(it%n)
        integer, intent(in) :: maxiter
        integer, intent(out) :: info
        procedure(it_callback), optional :: callback

        ! Local scalars
        real(wp) :: bnrm2, atol, gamma, eta, theta, eps, delta, beta
        real(wp) :: rho, rho_prev, xi, gamma_prev, theta_prev
        integer :: n, iter
        real(wp), allocatable :: tmp(:)

        real(wp), parameter :: rhotol   = epsilon(x)**2
        real(wp), parameter :: betatol  = rhotol
        real(wp), parameter :: gammatol = rhotol
        real(wp), parameter :: deltatol = rhotol
        real(wp), parameter :: epstol   = rhotol
        real(wp), parameter :: xitol    = rhotol

        n = it%n
        bnrm2 = it%norm(n,b)
        if (bnrm2 == zero) then
            info = 0
            return
        end if
        atol = max(it%atol,it%rtol*bnrm2)
        allocate(tmp(n)) ! FIXME

        associate(M1 => it%m1, M2 => it%m2, &
                  r=>it%r, q=>it%q, d=>it%d, s=>it%s, &
                  y=>it%y, ytilde=>it%ytilde, &
                  v=>it%v, vtilde=>it%vtilde, &
                  w=>it%w, wtilde=>it%wtilde, &
                  z=>it%z, ztilde=>it%ztilde, &
                  p=>it%p, ptilde=>it%ptilde)

        call it%matvec(n,x,r)
        r = b - r

        vtilde = r
        call M1%matvec(n,vtilde,y)
        rho = it%norm(n,y)
        wtilde = r
        call M2%rmatvec(n,wtilde,z)
        xi = it%norm(n,z)

        gamma =  one
        eta   = -one
        theta =  zero

        do iter = 0, maxiter - 1

            if (it%norm(n,r) < atol) then
                info = 0
                return
            end if

            if (abs(rho) < rhotol) then
                info = -10
                return
            end if

            if (abs(xi) < xitol) then
                info = -15
                return
            end if

            v = vtilde
            v = v*(one/rho)
            y = y*(one/rho)
            w = wtilde
            w = w*(one/xi)
            z = z*(one/xi)
            delta = it%dotprod(n,z,y)

            if (abs(delta) < deltatol) then
                info = -13
                return
            end if

            call M2%matvec(n,y,ytilde)
            call M1%rmatvec(n,z,ztilde)

            if (iter > 0) then
                ytilde = ytilde - (xi * (delta / eps)) * p
                p = ytilde
                ztilde = ztilde - (rho * (delta / eps)) * q
                q = ztilde
            else
                p = ytilde
                q = ztilde
            end if

            call it%matvec(n,p,ptilde)
            eps = it%dotprod(n,q,ptilde)
            if (abs(eps) < epstol) then
                info = -14
                return
            end if

            beta = eps / delta
            if (abs(beta) < betatol) then
                info = -11
                return
            end if

            vtilde = ptilde - beta*v
            call M1%matvec(n,vtilde,y)

            rho_prev = rho
            rho = it%norm(n,y)
            wtilde = w * (-beta)
            call it%rmatvec(n,q,tmp)  ! A @ q
            wtilde = wtilde + tmp
            call M2%rmatvec(n,wtilde,z)
            xi = it%norm(n,z)
            gamma_prev = gamma
            theta_prev = theta
            theta = rho / (gamma_prev * abs(beta))
            gamma = one / sqrt(one + theta**2)

            if (abs(gamma) < gammatol) then
                info = -12
                return
            end if

            eta = eta * ( -(rho_prev / beta) * (gamma/gamma_prev)**2 )

            if (iter > 0) then
                d = d * (theta_prev * gamma)**2 + eta * p
                s = s * (theta_prev * gamma)**2 + eta * ptilde
            else
                d = eta*p
                s = eta*ptilde
            end if

            x = x + d
            r = r - s

            if (present(callback)) call callback(x)

        end do

        info = iter

        end associate
    end subroutine

end submodule
