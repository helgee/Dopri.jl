program testrunner

use iso_c_binding, only: c_ptr, c_loc, c_f_pointer
use dopri

implicit none

double precision :: x
double precision :: xend
double precision, dimension(1) :: rtol
double precision, dimension(1) :: atol
double precision, dimension(6) :: y
double precision, dimension(6) :: y1
double precision, dimension(6) :: f
double precision, dimension(200) :: work
integer, dimension(100) :: iwork
integer :: n
integer :: itol
integer :: iout
integer :: idid
integer :: lwork
integer :: liwork
integer :: stat

type params
    double precision :: mu
end type params

type(params), target :: p
p%mu = 398600.4415d0
y = [-1814.0d0, -3708.0d0, 5153.0d0, 6.512d0, -4.229d0, -0.744d0]
xend = 5402.582703094263d0
x = 0d0
n = 6
itol = 0
iout = 1
atol = 1.4901161193847656d-8
rtol = 1d-6
lwork = 200
liwork = 100
work = 0d0
iwork = 0

call dopri5(n,newton,x,y,xend,rtol,atol,itol,solout,iout,work,lwork,iwork,liwork,c_loc(p),idid)

write(*,*) "ind1"
work = 0d0
iwork = 0
x = 0d0
y = [-1814.0d0, -3708.0d0, 5153.0d0, 6.512d0, -4.229d0, -0.744d0]
call dop853(n,newton,x,y,xend,rtol,atol,itol,solout,iout,work,lwork,iwork,liwork,c_loc(p),idid)

write(*,*) "ind2"
iout = 2
work = 0d0
iwork = 0
iwork(5) = n
x = 0d0
y = [-1814.0d0, -3708.0d0, 5153.0d0, 6.512d0, -4.229d0, -0.744d0]
call dopri5(n,newton,x,y,xend,rtol,atol,itol,solout5spec,iout,work,lwork,iwork,liwork,c_loc(p),idid)

write(*,*) "ind3"
work = 0d0
iwork = 0
iwork(5) = n
x = 0d0
y = [-1814.0d0, -3708.0d0, 5153.0d0, 6.512d0, -4.229d0, -0.744d0]
call dop853(n,newton,x,y,xend,rtol,atol,itol,solout8spec,iout,work,lwork,iwork,liwork,c_loc(p),idid)

write(*,*) "ind4"
work = 0d0
iwork = 0
iwork(5) = n
x = 0d0
y = [-1814.0d0, -3708.0d0, 5153.0d0, 6.512d0, -4.229d0, -0.744d0]
call dopri5(n,newton,x,y,xend,rtol,atol,itol,solout5all,iout,work,lwork,iwork,liwork,c_loc(p),idid)

write(*,*) "ind5"
work = 0d0
iwork = 0
iwork(5) = n
x = 0d0
y = [-1814.0d0, -3708.0d0, 5153.0d0, 6.512d0, -4.229d0, -0.744d0]
call dop853(n,newton,x,y,xend,rtol,atol,itol,solout8all,iout,work,lwork,iwork,liwork,c_loc(p),idid)

contains
    subroutine newton(n,x,y,f,pc)
        integer, intent(in) :: n
        double precision, intent(in) :: x
        double precision, dimension(n), intent(in) :: y
        type(c_ptr), intent(in) :: pc
        double precision, dimension(n), intent(out) :: f

        type(params), pointer :: pf
        double precision :: mu
        double precision :: r
        double precision :: r1

        call c_f_pointer(pc, pf)
        mu = pf%mu
        r = sqrt(sum(y(1:3)**2))
        f(1:3) = y(4:6)
        f(4:6) = -mu*y(1:3)/r**3
    end subroutine

    subroutine solout(nr, xold, x, y, n, con, icomp,&
        nd, pc, irtrn, xout)
                integer, intent(in) :: n
                integer, intent(in) :: nr
                integer, intent(in) :: nd
                integer, intent(in) :: irtrn
                integer, dimension(nd), intent(in) :: icomp
                double precision, intent(in) :: xold
                double precision, intent(in) :: x
                double precision, dimension(n), intent(in) :: y
                double precision, dimension(8*nd), intent(in) :: con
                double precision, intent(inout) :: xout
                type(c_ptr), intent(in) :: pc

                write(*,"(7E24.16)") x, y
    end subroutine

    subroutine solout5spec(nr, xold, x, y, n, con, icomp,&
        nd, pc, irtrn, xout)
                integer, intent(in) :: n
                integer, intent(in) :: nr
                integer, intent(in) :: nd
                integer, intent(in) :: irtrn
                integer, dimension(nd), intent(in) :: icomp
                double precision, intent(in) :: xold
                double precision, intent(in) :: x
                double precision, dimension(n), intent(in) :: y
                double precision, dimension(8*nd), intent(in) :: con
                double precision, intent(inout) :: xout
                type(c_ptr), intent(in) :: pc

                integer :: t
                integer :: told
                integer :: i
                integer :: j
                double precision, dimension(7) :: arr
                if (x > 1d0) then
                    told = ceiling(xold)
                    t = floor(x)
                    do i = told, t
                        arr(1) = dble(i)
                        do j = 1, 6
                            arr(j+1) = contd5(j, dble(i), con, icomp, nd)
                        end do
                        write(*,"(7E24.16)") arr
                    end do
                else if (x == 0d0) then
                    write(*,"(7E24.16)") x, y
                end if
    end subroutine

    subroutine solout8spec(nr, xold, x, y, n, con, icomp,&
        nd, pc, irtrn, xout)
                integer, intent(in) :: n
                integer, intent(in) :: nr
                integer, intent(in) :: nd
                integer, intent(in) :: irtrn
                integer, dimension(nd), intent(in) :: icomp
                double precision, intent(in) :: xold
                double precision, intent(in) :: x
                double precision, dimension(n), intent(in) :: y
                double precision, dimension(8*nd), intent(in) :: con
                double precision, intent(inout) :: xout
                type(c_ptr), intent(in) :: pc

                integer :: t
                integer :: told
                integer :: i
                integer :: j
                double precision, dimension(7) :: arr
                if (x > 1d0) then
                    told = ceiling(xold)
                    t = floor(x)
                    do i = told, t
                        arr(1) = dble(i)
                        do j = 1, 6
                            arr(j+1) = contd8(j, dble(i), con, icomp, nd)
                        end do
                        write(*,"(7E24.16)") arr
                    end do
                else if (x == 0d0) then
                    write(*,"(7E24.16)") x, y
                end if
    end subroutine

    subroutine solout5all(nr, xold, x, y, n, con, icomp,&
        nd, pc, irtrn, xout)
                integer, intent(in) :: n
                integer, intent(in) :: nr
                integer, intent(in) :: nd
                integer, intent(in) :: irtrn
                integer, dimension(nd), intent(in) :: icomp
                double precision, intent(in) :: xold
                double precision, intent(in) :: x
                double precision, dimension(n), intent(in) :: y
                double precision, dimension(8*nd), intent(in) :: con
                double precision, intent(inout) :: xout
                type(c_ptr), intent(in) :: pc

                integer :: t
                integer :: told
                integer :: i
                integer :: j
                double precision, dimension(7) :: arr
                if (x > 1d0) then
                    told = ceiling(xold)
                    t = floor(x)
                    do i = told, t
                        arr(1) = dble(i)
                        do j = 1, 6
                            arr(j+1) = contd5(j, dble(i), con, icomp, nd)
                        end do
                        write(*,"(7E24.16)") arr
                    end do
                end if
                write(*,"(7E24.16)") x, y
    end subroutine

    subroutine solout8all(nr, xold, x, y, n, con, icomp,&
        nd, pc, irtrn, xout)
                integer, intent(in) :: n
                integer, intent(in) :: nr
                integer, intent(in) :: nd
                integer, intent(in) :: irtrn
                integer, dimension(nd), intent(in) :: icomp
                double precision, intent(in) :: xold
                double precision, intent(in) :: x
                double precision, dimension(n), intent(in) :: y
                double precision, dimension(8*nd), intent(in) :: con
                double precision, intent(inout) :: xout
                type(c_ptr), intent(in) :: pc

                integer :: t
                integer :: told
                integer :: i
                integer :: j
                double precision, dimension(7) :: arr
                if (x > 1d0) then
                    told = ceiling(xold)
                    t = floor(x)
                    do i = told, t
                        arr(1) = dble(i)
                        do j = 1, 6
                            arr(j+1) = contd8(j, dble(i), con, icomp, nd)
                        end do
                        write(*,"(7E24.16)") arr
                    end do
                end if
                write(*,"(7E24.16)") x, y
    end subroutine
end program
