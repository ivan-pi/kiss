module kiss_tests
    use kiss
    implicit none
    public
contains

    subroutine bicgstab_test

        integer, parameter :: n = 4
        real(wp) :: A(n,n), b(n), x(n)
        integer :: info

        print '(A)', "BICGSTAB TEST"

        A = reshape([ real(wp) :: &
                     [4,2,0,1], &
                     [3,0,0,2], &
                     [0,1,1,1], &
                     [0,2,1,0]], shape=[n,n],order=[2,1])

        b = [-1.0_wp, -0.5_wp, -1.0_wp, 2.0_wp]

        x = 0
        call bicgstab_dense(A,b,x,atol=1.0e-5_wp,&
            callback=print_residual_norm,info=info)

        call print_residual_norm(x)
        write(*,'("info = ", I0)') info

    contains

        subroutine print_residual_norm(x)
            real(wp), intent(in) :: x(:)
            print *, "residual norm = ", norm2(b - matmul(A,x))
        end subroutine

    end subroutine


    subroutine bicg_test

        integer, parameter :: n = 3
        real(wp) :: A(n,n), b(n), x(n)
        integer :: info

        print '(A)', "BICG TEST"

        A = reshape([ real(wp) :: &
                     [3, 2,0], &
                     [1,-1,0], &
                     [0, 5,1]], shape=[n,n],order=[2,1])

        b = [ real(wp) :: 2, 4, -1 ]

        x = 0
        call bicg_dense(A,b,x,atol=1.0e-5_wp,&
            callback=print_residual_norm,info=info)

        call print_residual_norm(x)
        write(*,'("info = ", I0)') info

    contains

        subroutine print_residual_norm(x)
            real(wp), intent(in) :: x(:)
            print *, "residual norm = ", norm2(b - matmul(A,x))
        end subroutine

    end subroutine


    subroutine cg_test

        integer, parameter :: n = 4
        real(wp) :: A(n,n), b(n), x(n)
        integer :: info

        print '(A)', "CG TEST"

        A = reshape([ real(wp) :: &
                     [4, 0, 1, 0], &
                     [0, 5, 0, 0], &
                     [1, 0, 3, 2], &
                     [0, 0, 2, 4]], shape=[n,n],order=[2,1])

        b = [ -1.0_wp, -0.5_wp, -1.0_wp, 2.0_wp ]

        x = 0
        call cg_dense(A,b,x,atol=1.0e-5_wp,&
            callback=print_residual_norm,info=info)

        call print_residual_norm(x)
        write(*,'("info = ", I0)') info

    contains

        subroutine print_residual_norm(x)
            real(wp), intent(in) :: x(:)
            print *, "residual norm = ", norm2(b - matmul(A,x))
        end subroutine

    end subroutine


    subroutine cgs_test

        integer, parameter :: n = 4
        real(wp) :: A(n,n), b(n), x(n)
        integer :: info

        print '(A)', "CGS TEST"

        A = reshape([ real(wp) :: &
                     [4, 2, 0, 1], &
                     [3, 0, 0, 2], &
                     [0, 1, 1, 1], &
                     [0, 2, 1, 0]], shape=[n,n],order=[2,1])

        b = [ -1.0_wp, -0.5_wp, -1.0_wp, 2.0_wp ]

        x = 0
        call cgs_dense(A,b,x,atol=1.0e-5_wp,&
            callback=print_residual_norm,info=info)

        call print_residual_norm(x)
        write(*,'("info = ", I0)') info

    contains

        subroutine print_residual_norm(x)
            real(wp), intent(in) :: x(:)
            print *, "residual norm = ", norm2(b - matmul(A,x))
        end subroutine

    end subroutine


    subroutine gmres_test

        integer, parameter :: n = 3
        real(wp) :: A(n,n), b(n), x(n)
        integer :: info

        print '(A)', "GMRES TEST"

        A = reshape([ real(wp) :: &
                     [ 3, 2, 0], &
                     [ 1,-1, 0], &
                     [ 0, 5, 1]], shape=[n,n],order=[2,1])

        b = [ real(wp) :: 2, 4, -1]

        x = 0
        call gmres_dense(A,b,x,atol=1.0e-5_wp,&
            callback=print_residual_norm,info=info)

        call print_residual_norm(x)
        write(*,'("info = ", I0)') info

    contains

        subroutine print_residual_norm(x)
            real(wp), intent(in) :: x(:)
            print *, "residual norm = ", norm2(b - matmul(A,x))
        end subroutine

    end subroutine


    subroutine qmr_test

        integer, parameter :: n = 3
        real(wp) :: A(n,n), b(n), x(n)
        integer :: info

        print '(A)', "QMR TEST"

        A = reshape([ real(wp) :: &
                     [3, 2, 0], &
                     [1,-1, 0], &
                     [0, 5, 1]], shape=[n,n],order=[2,1])

        b = [ real(wp) :: 2, 4, -1]

        x = 0
        call qmr_dense(A,b,x,atol=1.0e-5_wp,&
            callback=print_residual_norm,info=info)

        call print_residual_norm(x)
        write(*,'("info = ", I0)') info

    contains

        subroutine print_residual_norm(x)
            real(wp), intent(in) :: x(:)
            print *, "residual norm = ", norm2(b - matmul(A,x))
        end subroutine

    end subroutine


    subroutine tfqmr_test

        integer, parameter :: n = 3
        real(wp) :: A(n,n), b(n), x(n)
        integer :: info

        print '(A)', "TFQMR TEST"

        A = reshape([ real(wp) :: &
                     [3, 2, 0], &
                     [1,-1, 0], &
                     [0, 5, 1]], shape=[n,n],order=[2,1])

        b = [ real(wp) :: 2, 4, -1]

        x = 0
        call tfqmr_dense(A,b,x,atol=0.0_wp,&
            callback=print_residual_norm,info=info)

        call print_residual_norm(x)
        write(*,'("info = ", I0)') info

    contains

        subroutine print_residual_norm(x)
            real(wp), intent(in) :: x(:)
            print *, "residual norm = ", norm2(b - matmul(A,x))
        end subroutine

    end subroutine
end module

program test_bicgstab

    use kiss_tests
    implicit none

    call bicg_test
    call bicgstab_test
    call cg_test
    call cgs_test
    call gmres_test
    call qmr_test
    call tfqmr_test

end program
