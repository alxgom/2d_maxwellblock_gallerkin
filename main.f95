program maxblock_gallerkin

    !18/03/2016

    !1-Program Scheme. Generate the tables needed for the integrations.
    !2-From the tables: the absissas for rho(from the generalized laguerre rule), and the absissas for phi(from the chebychev rule), and the weights.
    !3-Generate a table for the nonlinear term
    !4-Proyect my functions with the laguerre-chevy interpolation.
    !5-evolve all the coefficients,
    !6-for each time, keep the coeffs in a table. And keep the funcition for a fixed phi=0.(how good is this??)

    !1) -To do: Create the table, calling the programs. And find a way to read the values.
        !Problems: a) The binomial coefficient doesn-t work well if there is a factorial > 16
                  !solution: gamma(n+1)=n!
                  !b)
    !2) -To do: Make two array with the absissas and Rearange them.
    !3) -To do: Use the cuadrature to generate the nonlinear table.
    !4) -To do: Use the quadratures, and make a table for the coefficients(n,l,i).
    !5) -To do: Loop runge kutta for the coeffs.
                !almost finished
    !6) -To do: keep everything nice and clean.



    !write(ext, fmtext) data
    !character(len=4) :: ext
    !character(len=6), save :: fmtext = '(i4.4)'(integer 4 digitos.)
    !file=trim(dir) // '/' // fname// '.' // ext // '.bin'

    use resolution
    use constants
    implicit none


    real*8 ::  pump
    external pump
    real*8 :: pump_coefs(0:l,0:q,i)
    integer :: j,x,z,n,s,c,c2
    integer,parameter :: rho_dim=1000 , phi_dim=500
    real*8 :: rho(rho_dim), phi(phi_dim)
    real*8,dimension(rho_dim,phi_dim) :: pump_func, pump_proy
    real ( kind = 8 ), dimension(rho_dim ,0:l) :: cx !cx array calls from 0 in the second column.
    real ( kind = 8 ), dimension(rho_dim , phi_dim ,0:l) :: cy !cx array calls from 0 in the second column.

    call testing_rk_evol()! este codigo ya casi tiene lista la integracion con RK o Euler.
!    call quad_test2()
!    call rieman_integration()
!    call rieman_integration_BONtest()

    call rieman_proyection(pump,pump_coefs)

    call linspace(rho, 0.d0 , 8.d0 , rho_dim, 1)
    call linspace(phi, 0.d0 , 2*pi , phi_dim, 1)

    print*, pump_coefs
    pump_proy(:,:)=0
    do j=1,i
        do x=0,q
            call Alm_polynomial ( rho_dim, phi_dim, l, x, i, rho, phi, cx, cy )
            do z=1,l
                pump_proy(:,:)=pump_proy(:,:)+(pump_coefs(z,x,j)*cy(:,:,z-1))
            end do
        end do
    end do
    OPEN(2,file='pump_proy.in',form='unformatted',status='replace',access='stream')
        write(2), pump_proy
    close(2)

    pump_func(:,:)=0
    do c=1,rho_dim
        do c2=1,phi_dim
            pump_func(c,c2)=pump(rho(c),phi(c2))
        end do
    end do
    OPEN(2,file='pump_func.in',form='unformatted',status='replace',access='stream')
        write(2), pump_func
    close(2)
end program

subroutine rieman_proyection(func,coefs)
    use constants
    use funcs
    use resolution

    implicit none

    integer :: j,z,s,smax,it
    integer,parameter :: rho_dim=1000 , phi_dim=500
    real*8 :: rho(rho_dim), phi(phi_dim)
    real*8,intent(out) :: coefs(0:l,0:q,i)
    real ( kind = 8 ), dimension(rho_dim ,0:l) :: cx1  !cx array calls from 0 in the second column.
    real ( kind = 8 ), dimension(rho_dim , phi_dim ,0:l) :: cy1 !cx array calls from 0 in the second column.
    real*8 :: sumation
    real*8 :: func
    external func
    integer ( kind = 4 ) m,n,c,c2
    real ( kind = 8 ), allocatable :: x(:)
    real*8 integ
   !********************************************************************************
    call save_resolution()
   !********************************************************************************
    print*, l,q,i
    call linspace(rho, 0.d0 , 30.d0 , rho_dim, 1)
    call linspace(phi, 0.d0 , 2*pi , phi_dim, 1)
    do j=1,i
        do z=0,q

            call Alm_polynomial( rho_dim, phi_dim, l, q, i, rho, phi, cx1, cy1 )

            do s=0,l
                sumation=0
                do c2=1,phi_dim-1
                    do c=1,rho_dim-1
                        sumation=sumation + func(c,c2)*cy1(c,c2,s)*(rho(c+1)-rho(c))*(phi(c2+1)-phi(c2))*rho(c)
                    end do
                end do
                coefs(s,q,i)=sumation
            end do

        end do
    end do

    print*, 'rieman proyection with func: ', sumation

end subroutine


subroutine rieman_integration()
    use constants
    use funcs
    use resolution

    implicit none

    integer :: j,z,s,smax,tsteps,it

    integer,parameter :: rho_dim=1000 , phi_dim=120
    real*8 :: rho(rho_dim), phi(phi_dim)
    real ( kind = 8 ), dimension(rho_dim ,0:l) :: cx1, cx2  !cx array calls from 0 in the second column.
    real ( kind = 8 ), dimension(rho_dim , phi_dim ,0:l) :: cy1, cy2  !cx array calls from 0 in the second column.
    real*8 :: sumation
    real*8 :: pump, loss, testfunc, titas,erre
    external pump, loss, testfunc, titas,erre
    integer ( kind = 4 ) :: mm
    real ( kind = 8 ),allocatable :: cx3(:,:),cx4(:,:)  !cx array calls from 0 in the second column.
    integer ( kind = 4 ) m,n,c,c2
    real ( kind = 8 ), allocatable :: x(:)
    real*8 integ

   !************************************************
    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.
    call save_resolution()

  !  allocate(coefs_A(0:q,0:l,i),coefs_A0(0:q,0:l),coefs(0:q,0:l,0:i))
   !********************************************************************************
   !1-Program Scheme. Generate the tables needed for the integrations.

    !check Rlm:
    cx1(:,:)=0
    cx2(:,:)=0
    print*, l,q,i
    call linspace(rho, 0.d0 , 30.d0 , rho_dim, 1)
    call linspace(phi, 0.d0 , 2*pi , phi_dim, 1)
    OPEN(2,file='rho.in',form='unformatted',status='replace',access='stream')
        write(2), rho
    close(2)
    OPEN(2,file='phi.in',form='unformatted',status='replace',access='stream')
        write(2), phi
    close(2)


    call Alm_polynomial( rho_dim, phi_dim, l, 1, 1, rho, phi, cx1, cy1 )
    OPEN(2,file='Alm1.in',form='unformatted',status='replace',access='stream')
        write(2), cy1
    close(2)

    call Alm_polynomial( rho_dim, phi_dim, l, 1, 2, rho, phi, cx2, cy2 )
    OPEN(2,file='Alm2.in',form='unformatted',status='replace',access='stream')
        write(2), cy2
    close(2)

    sumation=0
    do c2=1,phi_dim-1
        do c=1,rho_dim-1
            sumation=sumation + cy1(c,c2,1)*cy2(c,c2,2)*(rho(c+1)-rho(c))*(phi(c2+1)-phi(c2))*rho(c)
        end do
    end do
    print*, 'integ con Alm m:', sumation

end subroutine




subroutine quad_test2()
    use constants
    use funcs
    use resolution

    implicit none

    !field Vars
 !   real(8), allocatable, dimension(:,:,:) :: coefs_A !coeficcients.  Table of  coefs(p,m,i)
 !   real(8), allocatable, dimension(:,:) :: coefs_A0 !coeficcients.  Table of  coefs(p,m,0)
    !functions
    !real(8), allocatable, dimension (:,:) :: loss, pump, Atest
    !space
 !   real*8 :: c0(l,q,i,5), c1(l,q,i,5)
  !  real*8, dimension(:), allocatable :: timetest
    integer :: j,z,s,smax,tsteps,it
  !  real*8 :: t0,tmin,tmax
  !  character(len=20) :: filename
 !   real*8,dimension(l,q,i) ::  phi_r,phi_i,pol_r,pol_i,pop
    integer,parameter :: rho_dim=1000 , phi_dim=120
    real*8 :: rho(rho_dim), phi(phi_dim)
    real ( kind = 8 ), dimension(rho_dim ,0:l) :: cx1, cx2  !cx array calls from 0 in the second column.
    real ( kind = 8 ), dimension(rho_dim , phi_dim ,0:l) :: cy1, cy2  !cx array calls from 0 in the second column.
    real*8 :: sumation
 !   real*8,dimension(rho_dim,phi_dim) :: intensity, er, ei, population
       !real ( kind = 8 ) :: dt=0.5
!    real(8) :: rn,phin
!    real(8) :: rho_0=0.,1
!    integer :: numkeep
    !laguerre vars
  !  real(8), allocatable, dimension(:,:,:) :: coefs !C(p,m,i) table
    real*8 :: pump, loss, testfunc, titas,erre
    external pump, loss, testfunc, titas,erre
!!    integer :: z,c,j
  !  real(8) :: tempvar
    integer ( kind = 4 ) :: mm
    real ( kind = 8 ),allocatable :: cx3(:,:),cx4(:,:)  !cx array calls from 0 in the second column.
    integer ( kind = 4 ) m,n,c,c2
    real ( kind = 8 ), allocatable :: x(:)
    !   real ( kind = 8 ), allocatable :: y(:)
  !  integer, parameter :: debug=1
    real*8 integ

   !************************************************
    call comparams()                             !parameters to compare with the expected solutions
    call saveparams()                            !saves the used parameters to a bin file, to be read by python.
    call save_resolution()

  !  allocate(coefs_A(0:q,0:l,i),coefs_A0(0:q,0:l),coefs(0:q,0:l,0:i))
   !********************************************************************************
   !1-Program Scheme. Generate the tables needed for the integrations.

    !check Rlm:
    cx1(:,:)=0
    cx2(:,:)=0
    print*, l,q,i
    call linspace(rho, 0.d0 , 30.d0 , rho_dim, 1)
    call linspace(phi, 0.d0 , 2*pi , phi_dim, 1)
    OPEN(2,file='rho.in',form='unformatted',status='replace',access='stream')
        write(2), rho
    close(2)
    OPEN(2,file='phi.in',form='unformatted',status='replace',access='stream')
        write(2), phi
    close(2)

!check with python
    call lm_polynomial( rho_dim, l, 2, rho, cx1 )
    OPEN(2,file='lm.in',form='unformatted',status='replace',access='stream')
        write(2), cx1
    close(2)!check with python
!check with python
    call Rlm_polynomial( rho_dim, l, 2, rho, cx1 )
    OPEN(2,file='Rlm.in',form='unformatted',status='replace',access='stream')
        write(2), cx1
    close(2)

    call Rlm_polynomial( rho_dim, l, 2, rho, cx2 )

    sumation=0
    do c=1,rho_dim-1
        sumation=sumation + cx1(c,1)*cx1(c,2)*(rho(c+1)-rho(c))*rho(c)
    end do
    print*, 'integ con m:', sumation


!check with python
    call lf_function( rho_dim, l, 1, rho, cx1 )
    OPEN(2,file='lm_alpha.in',form='unformatted',status='replace',access='stream')
        write(2), cx1
    close(2)
    call Rlm_alpha_polynomial( rho_dim, l, 1, rho, cx1 )
    OPEN(2,file='Rlm_alpha.in',form='unformatted',status='replace',access='stream')
        write(2), cx1
    close(2)

    sumation=0
    do c=1,rho_dim-1
        sumation=sumation + cx1(c,1)*cx1(c,3)*rho(c)*(rho(c+1)-rho(c))*rho(c)
    end do
    print*, 'integ con alpha', sumation


    sumation=0   !la integ de riemman funciona
    do c=1,phi_dim-1
        sumation=sumation + titas(phi(c))*(phi(c+1)-phi(c))
    end do
    print*, 'integ con titas', sumation


    call Alm_polynomial( rho_dim, phi_dim, l, 1, 1, rho, phi, cx1, cy1 )
    OPEN(2,file='Alm1.in',form='unformatted',status='replace',access='stream')
        write(2), cy1
    close(2)
!
    call Alm_polynomial( rho_dim, phi_dim, l, 1, 2, rho, phi, cx2, cy2 )
    OPEN(2,file='Alm2.in',form='unformatted',status='replace',access='stream')
        write(2), cy2
    close(2)
!
    sumation=0
    do c2=1,phi_dim-1
        do c=1,rho_dim-1
        sumation=sumation + cy1(c,c2,1)*cy2(c,c2,2)*(rho(c+1)-rho(c))*(phi(c2+1)-phi(c2))*rho(c)
        end do
    end do
    print*, 'integ con Alm m:', sumation


!    n=400
!    m=100
!    !call gen_laguerre_vars(n,m)
!   ! call gen_laguerre_quad(n,m)
! ! call laguerre_polynomial_test05 ( )
! ! call laguerre_polynomial_test08 ( 5, 0 )  !usar como base para mis integrales
!    mm=1
!    allocate(x(mm))
!    x(1)=1d0
!    allocate(cx3(mm,n+1))
!    call Rlm_polynomial ( mm, n, m, x, cx3 )
!!    print*, cx3
!    it=0
!    call gl_quadrature(erre,n,m,it,integ)
!    Print*,'deberia dar  sin jacobiano: 6.28321383214923'
!    print*, 'deberia dar con jacobiano: 3.2449440566614527'
!    !call gl_quadrature_cilindrical(loss,n,m,i,integ)
!
!    print*, ''
!    print*, 'pump test:', pump(10d0,0d0)
!    print*, 'loss test:', loss(10d0,0d0)
!    print*, 'testfunc test:', testfunc(20d0,10d0)


    !call quadrature_chev1(titas,100, integ)
   ! call quadrature_chev2(titas,10, integ) !bad abcisas
    !call quadrature_gauchv(titas,20000, integ)

end subroutine



double precision function erre(x)
    use resolution
    implicit none
    real*8, intent(in) :: x
    real*8, parameter :: m=1
    real*8, parameter :: n=1
    real*8, dimension(1 ,0:l-1) :: cx  !cx array calls from 0 in the second column.


    call Rlm_polynomial( 1, l, m, x, cx )
    erre=cx(1,n)
end function


!double precision function loss(x,y)
! implicit none
!    real*8, intent(in) :: x,y
!    real*8, parameter :: rho_0=1.d0
!
!    loss=5d0+4d0*tanh(5d0*(x-rho_0))
!end function

!double precision function testfunc(x,y)
! implicit none
!    real*8, intent(in) :: x,y
!    real*8, parameter :: rho_0=1.d0
!
!    testfunc=x
!end function

subroutine laguerre_polynomial_test05 ( )  !Lnm(x)

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_TEST05 tests LM_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) e
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ), parameter :: mm = 1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: v(:,:)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_TEST05:'
  write ( *, '(a)' ) '  LM_POLYNOMIAL_VALUES stores values of'
  write ( *, '(a)' ) '  the Laguerre polynomial Lm(n,m,x)'
  write ( *, '(a)' ) '  LM_POLYNOMIAL evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                                 Tabulated                 Computed'
  write ( *, '(a)' ) '     N     M        X            Lm(N,M,X)                 Lm(N,M,X)               Error'

  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lm_polynomial_values ( n_data, n, m, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( v(mm,n+1) )
    call lm_polynomial ( mm, n, m, x, v )
    fx2 = v(1,n+1)
    deallocate ( v )

    e = fx1 - fx2

    write ( *, '(2x,i4,2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, m, x, fx1, fx2, e

  end do

  return
end

subroutine laguerre_polynomial_test08 ( p, e )

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_TEST08 tests L_POWER_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polynomial
!    factors.
!
!    Input, integer ( kind = 4 ) E, the exponent of X.
!
  implicit none

  integer ( kind = 4 ) e
  integer ( kind = 4 ) p
  real ( kind = 8 ), allocatable :: table(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_TEST08'
  write ( *, '(a)' ) '  Compute a power product table for L(n,x):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( 0 <= x < +oo ) x^e L(i,x) L(j,x) exp(-x) dx'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponent of X, E = ', e

  allocate ( table(0:p,0:p) )

  call l_power_product ( p, e, table )

  call r8mat_print ( p + 1, p + 1, table, '  Power product table:' )

  deallocate ( table )

  return
end




!subroutine quad_test1()
!    use constants
!    use funcs
!    use resolution
!   ! use binom
!
!    implicit none
!
!    !integration vars
!    real(8) :: intime
!    real(8) :: timeinit=0.d0
!    real(8), allocatable, dimension(:) :: time
!    !real(8), dimension(9) :: yinit
!    integer :: numstep
!    integer :: indexkeep
!
!    !field Vars
!    real(8), allocatable, dimension(:) :: exr,exi,eyr,eyi,rxr,rxi,ryr,ryi,pop  !notation: exr-> real part of Ex electric field. rxr->real part of Rx polarization field. ->pop= population inversion.
!    real(8), allocatable, dimension(:) :: intensity_x2, intensity_y2, intensity
!    real(8), allocatable, dimension(:,:,:) :: coefs_A !coeficcients.  Table of  coefs(p,m,i)
!    real(8), allocatable, dimension(:,:) :: coefs_A0 !coeficcients.  Table of  coefs(p,m,0)
!
!    !functions
!  !  real(8), allocatable, dimension (:,:) :: loss, pump, Atest
!
!    !space
!    integer :: rho_dim=300
!    integer :: phi_dim=100
!    real(8), allocatable, dimension(:) :: rho, phi
!    real(8) :: rn,phin
!    real(8) :: rho_0=0.1
!    integer :: numkeep
!
!    !laguerre vars
!    real(8), allocatable, dimension(:,:,:) :: coefs !C(p,m,i) table
!!    integer :: p,m,i
!    real(8) :: ln, Alag
!    integer :: factorial,factorial_ros, choose
!    real(8) :: binomial
!    external factorial_ros, choose
!    real*8 :: pump, loss, testfunc, titas
!    external pump, loss, testfunc, titas
!    !real(8) :: ln
!
!!!    integer :: z,c,j
!    real(8) :: tempvar
!    integer, parameter :: debug=1
!
!    integer ( kind = 4 ) :: mm
!    integer ( kind = 4 ) :: n
!
!    real ( kind = 8 ),allocatable :: cx(:,:)  !cx array calls from 0 in the second column.
!    integer ( kind = 4 ) :: j,it
!    integer ( kind = 4 ) m
!
!    real ( kind = 8 ), allocatable :: x(:)
!    !   real ( kind = 8 ), allocatable :: y(:)
!
!    real*8 integ
!
!   !************************************************
!    call comparams()                             !parameters to compare with the expected solutions
!    call saveparams()                            !saves the used parameters to a bin file, to be read by python.
!
!    allocate(coefs_A(0:q,0:l,i),coefs_A0(0:q,0:l),coefs(0:q,0:l,0:i))
!   !********************************************************************************
!   !1-Program Scheme. Generate the tables needed for the integrations.
!
!    n=40
!    m=3
!
!    call gen_laguerre_vars(n,m)
!    call gen_laguerre_quad(n,m)
!
!  !  call laguerre_polynomial_test05 ( )
!  !  call laguerre_polynomial_test08 ( 5, 0 )  !usar como base para mis integrales
!
!    mm=1
!    allocate(x(mm))
!    x(1)=1d0
!
!    allocate(cx(mm,n+1))
!    call Rlm_polynomial ( mm, n, m, x, cx )
!    print*, cx
!
!    it=0
!    call gl_quadrature(pump,n,m,it,integ)
!    Print*,'deberia dar  sin jacobiano: 6.28321383214923'
!    print*, 'deberia dar con jacobiano: 3.2449440566614527'
!    !call gl_quadrature_cilindrical(loss,n,m,i,integ)
!
!    print*, ''
!    print*, 'pump test:', pump(10d0,0d0)
!    print*, 'loss test:', loss(10d0,0d0)
!    print*, 'testfunc test:', testfunc(20d0,10d0)
!
!
!    !call quadrature_chev1(titas,100, integ)
!   ! call quadrature_chev2(titas,10, integ) !bad abcisas
!    !call quadrature_gauchv(titas,20000, integ)
!
!end subroutine






!subroutine rungesfield(c0,c1)
!    use constants
!    use resolution
!    use rungekutta
!    implicit none
!
!    real*8, intent(inout) :: c0(l,q,i,9)
!    real*8, intent(out) :: c1(l,q,i,9)
!    real*8 :: rhs(l,q,i,9)
!    real*8 :: dt=0.005d0
!    integer :: j,x,z,order
!    real*8 :: rkcoef
!
!
!    do order=ord,1,-1
!        rkcoef=1./order
!        print*, order
!        call compute_RHS(c0,rhs)
!        !c0=c1+rhs*dt*rkcoef
!        c1=c0
!    end do
!end subroutine

!subroutine compute_RHS(c0,rhs)
!    use constants
!    use resolution
!    use rungekutta
!    implicit none
!
!    real*8, intent(inout) :: c0(l,q,i,9)
!    real*8, intent(out) :: rhs(l,q,i,9)
!    real*8 ::  dt=0.005d0
!    integer :: j,x,z
!
!    do j=1,i
!        do x=1,q
!            do z=1,l
!                rhs(z,x,j,1)=-k*c0(z,x,j,1)+(-a+d)*c0(z,x,j,2)-2*k*c0(z,x,j,5)
!                rhs(z,x,j,2)=-k*c0(z,x,j,2)+(a-d)*c0(z,x,j,1)-2*k*c0(z,x,j,6)
!                rhs(z,x,j,3)=-k*c0(z,x,j,3)+(-a+d)*c0(z,x,j,4)-2*k*c0(z,x,j,7)
!                rhs(z,x,j,4)=-k*c0(z,x,j,4)+(a-d)*c0(z,x,j,3)-2*k*c0(z,x,j,8)
!                rhs(z,x,j,5)=c0(z,x,j,5)-d*c0(z,x,j,6)
!                rhs(z,x,j,6)=c0(z,x,j,6)+d*c0(z,x,j,5)
!                rhs(z,x,j,7)=c0(z,x,j,7)-d*c0(z,x,j,8)
!                rhs(z,x,j,8)=c0(z,x,j,8)+d*c0(z,x,j,7)
!                rhs(z,x,j,9)=-g*(-c0(z,x,j,9))
!            end do
!        end do
!    end do
!
!end subroutine








!    allocate(rho(rho_dim))
!    call linspace(rho,0d0,20d0,rho_dim,1)
!    allocate(phi(phi_dim))
!    call linspace(phi,0d0,2*pi,phi_dim,1)
!
!    allocate(loss(phi_dim,rho_dim))!phi en las filas, rho en las colmnas
!    do c=1,rho_dim
!        do j=1,phi_dim
!             loss(j,c)=5d0+4d0*tanh(5d0*(rho(c)-rho_0))
!        end do
!    end do
!
!
!    allocate(pump(phi_dim,rho_dim))!phi en las filas, rho en las colmnas
!    do c=1,rho_dim
!        do j=1,phi_dim
!             pump(j,c)=5d0-tanh(5d0*(rho(c)-rho_0))/2d0
!        end do
!    end do
!
!!    call gen_laguerre_test() !generalized lagurre rule(quadrature)
!!    call cheby1()           !order n chebychev rule
!
!    print*,shape(loss)
!    print*, 'factorial', factorial(19)
!    print*, 'factorial_ros', factorial_ros(19)
!    print*, 'binomial', binomial(26,8)
!    print*, 'choose', choose(26,8)
!    print*, 'ln', ln(1,1,1d0)
!    print*, 'Alag', Alag(1,1,1,0d0,1d0)
!    !print*, 'binomial from stack', combo(8,15)
!    !test to plot some A
!    allocate(Atest(phi_dim,rho_dim))
!    do c=1,rho_dim
!        do j=1, phi_dim
!            Atest(j,c)=Alag(0,0,0,rho(c),phi(j))
!        end do
!    end do
!
!    print*, shape(Atest)
!    open(1,file='atest.in', form='unformatted')
!        write(1) Atest
!    close(1)
! 43727855 sistem.
! 43124191 rochela carlos
