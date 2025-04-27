module kiss_linop_tests
    use kiss
    implicit none
    public

    type, extends(kiss_linop) :: dense_op
        real(wp), allocatable :: A(:,:)
    contains
        procedure :: apply => apply_dense
    end type

    !type, extends(kiss_linop) :: stencil_op
    !contains
    !    procedure :: apply => apply_dense
    !end type

contains

    subroutine apply_dense(op,x,y)
        class(dense_op), intent(in) :: op
        real(wp), intent(in) :: x(op%n)
        real(wp), intent(out) :: y(op%n)
        y = matmul(op%A,x)
    end subroutine

    subroutine bicgstab_test

        integer, parameter :: n = 4
        type(dense_op) :: A
        real(wp) :: b(n), x(n)
        integer :: info

        print '(A)', "BICGSTAB TEST"

        A = dense_op(n,reshape([ real(wp) :: &
                     [4,2,0,1], &
                     [3,0,0,2], &
                     [0,1,1,1], &
                     [0,2,1,0]], shape=[4,4],order=[2,1]))

        b = [-1.0_wp, -0.5_wp, -1.0_wp, 2.0_wp]

        x = 0
        call bicgstab_linop(A,b,x,atol=1.0e-5_wp,&
            callback=print_residual_norm,info=info)

        call print_residual_norm(x)
        write(*,'("info = ", I0)') info

    contains

        subroutine print_residual_norm(x)
            real(wp), intent(in) :: x(:)
            real(wp) :: ax(size(x))
            call A%apply(x,ax)
            print *, "residual norm = ", norm2(b - ax)
        end subroutine

    end subroutine

end module

program test_linop

    use kiss_linop_tests
    implicit none

    call bicgstab_test

end program