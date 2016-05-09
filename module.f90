 !*********************************************************************************
  MODULE constants
    implicit none
    real(8), parameter :: a=2d0
    real(8), parameter :: gperp=10**8d0 !#gamma perpendicular, loss rate
    real(8), parameter :: tempscale=1*(10.**6)/gperp !#scale to micro seconds
    real(8), parameter :: wscale=1000*gperp/(10.**6) !#scale frequency to khz
    real(8), parameter :: mu=0.25d0/10**4, Dphi0=0.0d0
    real(8), parameter :: k=0.9*10**7d0/gperp, g=2.5*10.**4/gperp, D0=a*k/mu!, w_res=sqrt(k*g*((D0*mu/k)-1.))*wscale, w=sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale, atest=D0*mu/k,
    real(8) :: d=1.0d0
    !real(8) :: m=0.02d0
    real(8) :: wf=0.0045d0
    real(8) :: w_res, w!, atest
    integer, parameter :: savefile=1


! TOL : Tolerance for functions and Romberg integration
! NRAD: Number of points used in radial Gauss quadratures
! NANG: Number of points used in angular Gauss quadratures
!
      DOUBLE COMPLEX, PARAMETER   :: IM = (0.,1.)
      real(8), PARAMETER :: pi = 3.141592653589793d0
      DOUBLE PRECISION, PARAMETER :: TOL = 1.d-12
      DOUBLE PRECISION, PARAMETER :: TINY = 1.d-12
      INTEGER, PARAMETER          :: NRAD = 87
      INTEGER, PARAMETER          :: NANG = 87
      SAVE

    contains

    subroutine comparams()
    !'''parameters to compare with the results'''
        w_res=sqrt(k*g*((D0*mu/k)-1.))*wscale !#resonance frequency
        !atest=D0*mu/k
        w=sqrt(k*g*(a-1.)-(g*g*a*a)/4)*wscale !#Relaxation oscilations frequency
    end subroutine

    subroutine saveparams()
        if (savefile.eq.1) then
            open (1,file="scales.in",form='unformatted')
                write(1)  wf*wscale, Dphi0, w_res , k, mu, d, g, D0, a, wf, wscale, tempscale
            close (1)
        endif
    end subroutine

  END MODULE constants
 !************************************************************************************
  MODULE rungekutta
!
! ord: order of the Runge-Kutta time integration
!
      INTEGER, PARAMETER :: ord = 4
      SAVE

  END MODULE
  !**************************************************************************************

 module funcs
    implicit none
    public :: linspace

    contains

    subroutine linspace(x,x_start, x_end, x_len,dir)
    !******************************************************************************
    !linearly spaced array named x, from x_start value to x_end, with #x_len values.
    !dir=1 ---> from min to max.    dir=2 ---> from max to min.
    !******************************************************************************
        real(8), dimension(:), intent(inout) :: x
        real(8) :: x_start, x_end
        integer :: x_len, i
        real(8) :: dx
        integer :: dir

        dx=(x_end - x_start)/(x_len-1)
        if (dir.eq.1) then
            do i=1,x_len
                x(i)=x_start+(i-1)*dx
            end do
        end if
        if (dir.eq.2) then
            do i=1,x_len
                x(i)=x_end-(i-1)*dx
            end do
        end if
        !x(1:x_len)=[(x_start+(i-1)*dx)),i=1,x_len]
    end subroutine
end module

  MODULE resolution !q and l to be called and changed as seen fit.

! q: maximum number of zeros in the radial direction
! l: maximum number of l in the angular direction
    implicit none
      INTEGER, PARAMETER :: q = 5
      INTEGER, PARAMETER :: l = 5
      INTEGER, PARAMETER :: i = 2
      SAVE

      contains

      subroutine save_resolution()
      open (unit=2,file='resolution.txt',action="write",status='replace')
          write (2,*) "q=",q
          write (2,*) "l=",l
          write (2,*) "i=",i
      close(2)
      end subroutine



  END MODULE
