module kiss_mkl_tests

    use kiss_mkl
    implicit none
    public

contains

    subroutine pcg_test

        integer, parameter :: n = 4
        real(dp) :: A(n,n), b(n), x(n)
        integer :: info

        type(rci_iss_cache) :: cache

        print '(A)', "PCG TEST"

        A = reshape([ real(dp) :: &
                     [4, 0, 1, 0], &
                     [0, 5, 0, 0], &
                     [1, 0, 3, 2], &
                     [0, 0, 2, 4]], shape=[n,n],order=[2,1])

        b = [ -1.0d0, -0.5d0, -1.0d0, 2.0d0 ]

        x = 0

        cache = new_cg_cache(n,pc=.false.)
        cache%matvec => matvec

        print *, allocated(cache%tmp)
        call pcg(cache,x,b,info,atol=1.0d-5)

        call print_residual_norm(x)
        write(*,'("info = ", I0)') info

    contains

        subroutine matvec(n,x,y)
            integer, intent(in) :: n
            real(dp), intent(in) :: x(n)
            real(dp), intent(out) :: y(n)
            y = matmul(A(1:n,1:n),x)
        end subroutine

        subroutine print_residual_norm(x)
            real(dp), intent(in) :: x(:)
            print *, "residual norm = ", norm2(b - matmul(A,x))
        end subroutine

    end subroutine


    subroutine fgmres_test

        integer, parameter :: n = 3
        real(dp) :: A(n,n), b(n), x(n)
        integer :: info

        type(rci_iss_cache) :: cache

        print '(A)', "FGMRES TEST"

        A = reshape([ real(dp) :: &
                     [3, 2, 0], &
                     [1,-1, 0], &
                     [0, 5, 1]], shape=[n,n],order=[2,1])

        b = [ real(dp) :: 2, 4, -1]

        x = 0

        cache = new_fgmres_cache(n,nrestart=min(20,n))
        cache%matvec => matvec

        print *, allocated(cache%tmp)
        call fgmres(cache,x,b,info,atol=1.0d-5)

        call print_residual_norm(x)
        write(*,'("info = ", I0)') info

    contains

        subroutine matvec(n,x,y)
            integer, intent(in) :: n
            real(dp), intent(in) :: x(n)
            real(dp), intent(out) :: y(n)
            y = matmul(A(1:n,1:n),x)
        end subroutine

        subroutine print_residual_norm(x)
            real(dp), intent(in) :: x(:)
            print *, "residual norm = ", norm2(b - matmul(A,x))
        end subroutine

    end subroutine
end module

program test_bicgstab

    use kiss_mkl_tests
    implicit none

    call pcg_test
    call fgmres_test

end program