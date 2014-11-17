module dopri

use iso_c_binding, only: c_double, c_int, c_funptr, c_f_procpointer

implicit none

interface
    subroutine dop853(n, fcn, x, y, xend, rtol, atol,&
            itol, solout, iout, work, lwork, iwork,&
            liwork, rpar, ipar, idid)
        interface
            subroutine fcn(n, x, y, f, rpar, ipar)
                integer, intent(in) :: n
                integer, dimension(*),intent(inout) :: ipar
                double precision, intent(in) :: x
                double precision, dimension(n), intent(in) :: y
                double precision, dimension(*), intent(inout) :: rpar
                double precision, dimension(n), intent(out) :: f
            end subroutine fcn
            subroutine solout(nr, xold, x, y, n, con, icomp,&
                    nd, rpar, ipar, irtrn, xout)
                integer, intent(in) :: n
                integer, intent(in) :: nr
                integer, intent(in) :: nd
                integer, intent(in) :: irtrn
                integer, dimension(*), intent(inout) :: ipar
                integer, dimension(nd), intent(in) :: icomp
                double precision, intent(in) :: xold
                double precision, intent(in) :: x
                double precision, dimension(n), intent(in) :: y
                double precision, dimension(8*nd), intent(in) :: con
                double precision, dimension(*), intent(inout) :: rpar
                double precision, intent(inout) :: xout
            end subroutine solout
        end interface
        integer, intent(in) :: n
        integer, intent(in) :: itol
        integer, intent(in) :: iout
        integer, intent(in) :: lwork
        integer, intent(in) :: liwork
        integer, dimension(*),intent(inout) :: ipar
        integer, intent(out) :: idid
        double precision, intent(in) :: xend
        double precision, dimension(n), intent(in) :: rtol
        double precision, dimension(n), intent(in) :: atol
        double precision, dimension(*), intent(inout) :: rpar
        double precision, dimension(lwork), intent(inout) :: work
        integer, dimension(liwork), intent(inout) :: iwork
        double precision, intent(inout) :: x
        double precision, dimension(n), intent(inout) :: y
    end subroutine dop853

    double precision function contd8(ii, x, con, icomp, nd)
        integer, intent(in) :: ii
        double precision, intent(in) :: x
        double precision, dimension(8*nd), intent(in) :: con
        integer, dimension(nd), intent(in) :: icomp
        integer, intent(in) :: nd
    end function contd8

    subroutine dopri5(n, fcn, x, y, xend, rtol, atol,&
            itol, solout, iout, work, lwork, iwork,&
            liwork, rpar, ipar, idid)
        interface
            subroutine fcn(n, x, y, f, rpar, ipar)
                integer, intent(in) :: n
                integer, dimension(*),intent(inout) :: ipar
                double precision, intent(in) :: x
                double precision, dimension(n), intent(in) :: y
                double precision, dimension(*), intent(inout) :: rpar
                double precision, dimension(n), intent(out) :: f
            end subroutine fcn
            subroutine solout(nr, xold, x, y, n, con, icomp,&
                    nd, rpar, ipar, irtrn, xout)
                integer, intent(in) :: n
                integer, intent(in) :: nr
                integer, intent(in) :: nd
                integer, intent(in) :: irtrn
                integer, dimension(*), intent(inout) :: ipar
                integer, dimension(nd), intent(in) :: icomp
                double precision, intent(in) :: xold
                double precision, intent(in) :: x
                double precision, dimension(n), intent(in) :: y
                double precision, dimension(8*nd), intent(in) :: con
                double precision, dimension(*), intent(inout) :: rpar
                double precision, intent(inout) :: xout
            end subroutine solout
        end interface
        integer, intent(in) :: n
        integer, intent(in) :: itol
        integer, intent(in) :: iout
        integer, intent(in) :: lwork
        integer, intent(in) :: liwork
        integer, dimension(*),intent(inout) :: ipar
        integer, intent(out) :: idid
        double precision, intent(in) :: xend
        double precision, dimension(n), intent(in) :: rtol
        double precision, dimension(n), intent(in) :: atol
        double precision, dimension(*), intent(inout) :: rpar
        double precision, dimension(lwork), intent(inout) :: work
        integer, dimension(liwork), intent(inout) :: iwork
        double precision, intent(inout) :: x
        double precision, dimension(n), intent(inout) :: y
    end subroutine dopri5

    double precision function contd5(ii, x, con, icomp, nd)
        integer, intent(in) :: ii
        double precision, intent(in) :: x
        double precision, dimension(5*nd), intent(in) :: con
        integer, dimension(nd), intent(in) :: icomp
        integer, intent(in) :: nd
    end function contd5
end interface

abstract interface
    subroutine c_fcn(n, x, y, f, rpar, ipar)
        import :: c_int
        import :: c_double
        integer(c_int), intent(in) :: n
        integer(c_int), dimension(*),intent(inout) :: ipar
        real(c_double), intent(in) :: x
        real(c_double), dimension(n), intent(in) :: y
        real(c_double), dimension(*), intent(inout) :: rpar
        real(c_double), dimension(n), intent(out) :: f
    end subroutine c_fcn
    subroutine c_solout(nr, xold, x, y, n, con, icomp,&
            nd, rpar, ipar, irtrn, xout)
        import :: c_int
        import :: c_double
        integer(c_int), intent(in) :: n
        integer(c_int), intent(in) :: nr
        integer(c_int), intent(in) :: nd
        integer(c_int), intent(in) :: irtrn
        integer(c_int), dimension(*), intent(inout) :: ipar
        integer(c_int), dimension(nd), intent(in) :: icomp
        real(c_double), intent(in) :: xold
        real(c_double), intent(in) :: x
        real(c_double), dimension(n), intent(in) :: y
        real(c_double), dimension(8*nd), intent(in) :: con
        real(c_double), dimension(*), intent(inout) :: rpar
        real(c_double), intent(inout) :: xout
    end subroutine c_solout
end interface

contains

subroutine c_dop853(n, cfcn, x, y, xend, rtol, atol,&
        itol, csolout, iout, work, lwork, iwork,&
        liwork, rpar, ipar, idid) bind(c)
    integer(c_int), intent(in) :: n
    type(c_funptr), intent(in), value :: cfcn
    real(c_double), intent(inout) :: x
    real(c_double), dimension(n), intent(inout) :: y
    real(c_double), intent(in) :: xend
    real(c_double), dimension(n), intent(in) :: rtol
    real(c_double), dimension(n), intent(in) :: atol
    integer(c_int), intent(in) :: itol
    type(c_funptr), intent(in), value :: csolout
    integer(c_int), intent(in) :: iout
    real(c_double), dimension(lwork), intent(inout) :: work
    integer(c_int), intent(in) :: lwork
    integer(c_int), dimension(liwork), intent(inout) :: iwork
    integer(c_int), intent(in) :: liwork
    real(c_double), dimension(*), intent(inout) :: rpar
    integer(c_int), dimension(*),intent(inout) :: ipar
    integer(c_int), intent(out) :: idid

    procedure(c_fcn), pointer :: fcn
    procedure(c_solout), pointer :: solout

    real(c_double), dimension(n) :: f

    call c_f_procpointer(cfcn, fcn)
    call c_f_procpointer(csolout, solout)

    call dop853(n, fcn, x, y, xend, rtol, atol,&
            itol, solout, iout, work, lwork, iwork,&
            liwork, rpar, ipar, idid)
end subroutine c_dop853

function c_contd8(ii, x, con, icomp, nd) result(ret) bind(c)
    real(c_double) :: ret
    integer(c_int), intent(in) :: ii
    real(c_double), intent(in) :: x
    real(c_double), dimension(8*nd), intent(in) :: con
    integer(c_int), dimension(nd), intent(in) :: icomp
    integer(c_int), intent(in) :: nd
    ret = contd8(ii, x, con, icomp, nd)
end function c_contd8

subroutine c_dopri5(n, cfcn, x, y, xend, rtol, atol,&
        itol, csolout, iout, work, lwork, iwork,&
        liwork, rpar, ipar, idid) bind(c)
    type(c_funptr), intent(in), value :: cfcn
    type(c_funptr), intent(in), value :: csolout
    integer(c_int), intent(in) :: n
    integer(c_int), intent(in) :: itol
    integer(c_int), intent(in) :: iout
    integer(c_int), intent(in) :: lwork
    integer(c_int), intent(in) :: liwork
    integer(c_int), dimension(*),intent(inout) :: ipar
    integer(c_int), intent(out) :: idid
    real(c_double), intent(in) :: xend
    real(c_double), dimension(n), intent(in) :: rtol
    real(c_double), dimension(n), intent(in) :: atol
    real(c_double), dimension(*), intent(inout) :: rpar
    real(c_double), dimension(lwork), intent(inout) :: work
    integer(c_int), dimension(liwork), intent(inout) :: iwork
    real(c_double), intent(inout) :: x
    real(c_double), dimension(n), intent(inout) :: y

    procedure(c_fcn), pointer :: fcn
    procedure(c_solout), pointer :: solout

    call c_f_procpointer(cfcn, fcn)
    call c_f_procpointer(csolout, solout)
    call dopri5(n, fcn, x, y, xend, rtol, atol,&
            itol, solout, iout, work, lwork, iwork,&
            liwork, rpar, ipar, idid)
end subroutine c_dopri5

function c_contd5(ii, x, con, icomp, nd) result(ret) bind(c)
    real(c_double) :: ret
    integer(c_int), intent(in) :: ii
    real(c_double), intent(in) :: x
    real(c_double), dimension(5*nd), intent(in) :: con
    integer(c_int), dimension(nd), intent(in) :: icomp
    integer(c_int), intent(in) :: nd
    ret = contd5(ii, x, con, icomp, nd)
end function c_contd5

end module
