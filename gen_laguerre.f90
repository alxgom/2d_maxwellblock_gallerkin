
!************************ My subroutines***********************
subroutine gen_laguerre_vars(p,m)

!*****************************************************************************80
!
!! MAIN is the main program for GEN_LAGUERRE_RULE.
!
!  Discussion:
!
!    This program computes a generalized Gauss-Laguerre quadrature rule
!    and writes it to a file.
!
!    The user specifies: AS INPUT
!    * the ORDER (number of points) in the rule;
!    * ALPHA, the exponent of |X|;
!    * A, the left endpoint;
!    * B, the scale factor;
!    * FILENAME, the root name of the output files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21/03/2016
!  Author:
!
!    John Burkardt, mod by alexis gomel
!
  implicit none

  integer,intent(in) :: p, m
  real ( kind = 8 ) :: a=0
  real ( kind = 8 ) alpha
 ! integer ( kind = 4 ) arg_num
  real ( kind = 8 ) :: b=1
  real ( kind = 8 ) beta
  character ( len = 255 )filename
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) last
  integer ( kind = 4 ) order
  real ( kind = 8 ) r(2)
  real ( kind = 8 ) r8_huge
  character ( len = 255 ) string
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  !call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEN_LAGUERRE_RULE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Compute a generalized Gauss-Laguerre rule for approximating'
  write ( *, '(a)' ) &
    '    Integral ( a <= x < +oo ) |x-a|^ALPHA exp(-B*(x-a)) f(x) dx'
  write ( *, '(a)' ) '  of order ORDER.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The user specifies ORDER, ALPHA, A, B, and FILENAME.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ORDER is the number of points.'
  write ( *, '(a)' ) '  ALPHA is the exponent of |X|.'
  write ( *, '(a)' ) '  A is the left endpoint, (typically 0).'
  write ( *, '(a)' ) '  B is the scale factor, (typically 1).'
  write ( *, '(a)' ) '  FILENAME is used to generate 3 files:'
  write ( *, '(a)' ) '  * filename_w.txt - the weight file'
  write ( *, '(a)' ) '  * filename_x.txt - the abscissa file.'
  write ( *, '(a)' ) '  * filename_r.txt - the region file.'

  alpha=m
  order=p
!
!  Initialize parameters.
!
  beta = 0.0D+00
  filename='laguerre_vars'
!
!  Input summary.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ORDER = ', order
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a)' ) '  FILENAME = "' // trim ( filename ) // '".'
!
!  Construct the rule.
!
  allocate ( w(order) )
  allocate ( x(order) )

  kind = 5
  call cgqf ( order, kind, alpha, beta, a, b, x, w )


!
!  Write the rule.
!
  r(1) = a
  r(2) = r8_huge ( )

  call rule_write ( order, x, w, r, filename )
!
!  Free memory.
!
  deallocate ( w )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEN_LAGUERRE_RULE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  !call timestamp ( )

  return
end subroutine
subroutine gen_laguerre_quad(p,m)

!*****************************************************************************80
!
!! MAIN is the main program for GEN_LAGUERRE_RULE.
!
!  Discussion:
!
!    This program computes a generalized Gauss-Laguerre quadrature rule
!    and writes it to a file.
!
!    The user specifies: AS INPUT
!    * the ORDER (number of points) in the rule;
!    * ALPHA, the exponent of |X|;
!    * A, the left endpoint;
!    * B, the scale factor;
!    * FILENAME, the root name of the output files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21/03/2016
!  Author:
!
!    John Burkardt, mod by alexis gomel
!
  implicit none

  integer,intent(in) :: p, m
  real ( kind = 8 ) :: a=0
  real ( kind = 8 ) alpha
 ! integer ( kind = 4 ) arg_num
  real ( kind = 8 ) :: b=1
  real ( kind = 8 ) beta
  character ( len = 255 )filename
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) last
  integer ( kind = 4 ) order
  real ( kind = 8 ) r(2)
  real ( kind = 8 ) r8_huge
  character ( len = 255 ) string
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  !call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEN_LAGUERRE_quad'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Compute a generalized Gauss-Laguerre rule for approximating'
  write ( *, '(a)' ) &
    '    Integral ( a <= x < +oo ) |x-a|^ALPHA exp(-B*(x-a)) f(x) dx'
  write ( *, '(a)' ) '  of order ORDER.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The user specifies ORDER, ALPHA, A, B, and FILENAME.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ORDER is the number of points.'
  write ( *, '(a)' ) '  ALPHA is the exponent of |X|.'
  write ( *, '(a)' ) '  A is the left endpoint, (typically 0).'
  write ( *, '(a)' ) '  B is the scale factor, (typically 1).'
  write ( *, '(a)' ) '  FILENAME is used to generate 3 files:'
  write ( *, '(a)' ) '  * filename_w.txt - the weight file'
  write ( *, '(a)' ) '  * filename_x.txt - the abscissa file.'
  write ( *, '(a)' ) '  * filename_r.txt - the region file.'

  alpha=m
  order=p
!
!  Initialize parameters.
!
  beta = 0.0D+00
  filename='laguerre_quad'
!
!  Input summary.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ORDER = ', order
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a)' ) '  FILENAME = "' // trim ( filename ) // '".'
!
!  Construct the rule.
!
  allocate ( w(order) )
  allocate ( x(order) )

  kind = 5
  call cdgqf ( order, kind, alpha, beta, a, b, x, w )


!
!  Write the rule.
!
  r(1) = a
  r(2) = r8_huge ( )

  call rule_write ( order, x, w, r, filename )
!
!  Free memory.
!
  deallocate ( w )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEN_LAGUERRE_quad:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  !call timestamp ( )

  return
end subroutine
subroutine gen_laguerre_integ(func,p,m)

!*****************************************************************************80
!
!! MAIN is the main program for GEN_LAGUERRE_RULE.
!
!  Discussion:
!
!    This program computes a generalized Gauss-Laguerre quadrature rule
!    and writes it to a file.
!
!    The user specifies: AS INPUT
!    * the ORDER (number of points) in the rule;
!    * ALPHA, the exponent of |X|;
!    * A, the left endpoint;
!    * B, the scale factor;
!    * FILENAME, the root name of the output files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21/03/2016
!  Author:
!
!    John Burkardt, mod by alexis gomel
!
  implicit none

  integer,intent(in) :: p, m
  real ( kind = 8 ) :: a=0
  real ( kind = 8 ) alpha
 ! integer ( kind = 4 ) arg_num
  real ( kind = 8 ) :: b=1
  real ( kind = 8 ) beta
  character ( len = 255 )filename
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) last
  integer ( kind = 4 ) order
  real ( kind = 8 ) r(2)
  real ( kind = 8 ) r8_huge
  character ( len = 255 ) string
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  external func
  !call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEN_LAGUERRE_RULE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Compute a generalized Gauss-Laguerre rule for approximating'
  write ( *, '(a)' ) &
    '    Integral ( a <= x < +oo ) |x-a|^ALPHA exp(-B*(x-a)) f(x) dx'
  write ( *, '(a)' ) '  of order ORDER.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The user specifies ORDER, ALPHA, A, B, and FILENAME.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ORDER is the number of points.'
  write ( *, '(a)' ) '  ALPHA is the exponent of |X|.'
  write ( *, '(a)' ) '  A is the left endpoint, (typically 0).'
  write ( *, '(a)' ) '  B is the scale factor, (typically 1).'
  write ( *, '(a)' ) '  FILENAME is used to generate 3 files:'
  write ( *, '(a)' ) '  * filename_w.txt - the weight file'
  write ( *, '(a)' ) '  * filename_x.txt - the abscissa file.'
  write ( *, '(a)' ) '  * filename_r.txt - the region file.'

  alpha=m
  order=p
!
!  Initialize parameters.
!
  beta = 0.0D+00
  filename='laguerre_vars'
!
!  Input summary.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ORDER = ', order
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a)' ) '  FILENAME = "' // trim ( filename ) // '".'
!
!  Construct the rule.
!
  allocate ( w(order) )
  allocate ( x(order) )

  kind = 5
  call cgqf ( order, kind, alpha, beta, a, b, x, w )

!
!  Write the rule.
!
  r(1) = a
  r(2) = r8_huge ( )

  call rule_write ( order, x, w, r, filename )
!
!  Free memory.
!
  deallocate ( w )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEN_LAGUERRE_RULE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  !call timestamp ( )

  return
end subroutine
subroutine Rlm_alpha_polynomial ( mm, n, m, x, cx ) !Rnm=2*((2x**2)**(m/2))*((facoriales)**(1/2))*exp(-x**2)*Lnm(2x**2) gives 0 to n Rjm(x)

!*****************************************************************************
!  Special values:
!
!    For M = 0, the associated Laguerre polynomials Lm(N,M,X) are equal
!    to the Laguerre polynomials L(N,X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22/3/2016
!
!  Author:
!
!    John Burkardt mod Alexis Gomel
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, integer ( kind = 4 ) M, the parameter.  M must be nonnegative.
!
!    Input, real ( kind = 8 ) X(MM), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the associated Laguerre polynomials
!    of degrees 0 through N evaluated at the evaluation points.
!
    implicit none

    integer ( kind = 4 ) mm
    integer ( kind = 4 ),intent(in) :: n,m
    real ( kind = 8 ) cx(mm,0:n)  !cx array calls from 0 in the second column.
    integer ( kind = 4 ) i,j
    real ( kind = 8 ) x(mm)
    real ( kind = 8 ) y(mm)
    real*8 div

    y=x**2
    call  lf_function ( mm, n, m, 2*y,cx )
    do j=0,n !take a look at how i implemet the do.
        div=1./gamma(j+m+1.)
        do i=1,mm
            cx(i,j)= 2*((2*y(i))**(m/2))*sqrt(gamma(j+1.)*div)*cx(i,j)*exp(-y(i))
        end do
    end do
    !  return
end subroutine
subroutine gl_quadrature(func,n,m,i, integral)
    use constants
    implicit none
    real*8 func
    integer, intent(in) :: n,m,i
    real*8,allocatable, dimension(:) :: x, xx, y, yy , test, wtest
    real*8,allocatable :: wx(:), wy(:)
    real*8, intent(out) :: integral
    integer c,j, mtest
    integer ( kind = 4 ) kind
    real ( kind = 8 ),allocatable :: cx(:,:)  !cx array calls from 0 in the second column.

    allocate (x(n),wx(n),xx(n+1), y(m), wy(m), test(n), wtest(n))
    xx(1)=0

    !zeros Lnm:
    call lm_quadrature_rule ( n, m, x, w )
    kind=2
    call cgqf (m, kind, 0.d0, 0.d0, 0d0, 2*pi, y, wy )

 !   print*,'laguerre abcisas:', x
!    print*,'chebychev abcisas:', y

!    kind=5
!    call cgqf (n, kind, 5d0 , 0.d0, 0do, 1d0, test, wtest )
!    print*,'laguerre abcisas test:', test

    integral=0.d0
    allocate(cx(n,n+1))
    call Rlm_polynomial ( n, n, m, x, cx )

    do c=1,n
        do j=1,m
            integral= integral + func(x(c),y(j))*wy(j)*wx(c)!*cx(c,m)  !le flta los de laguere.
        end do
    end do

    print*, 'integral: ', integral

end subroutine
subroutine gl_quadrature_cilindrical(func,n,m,i, integral) !for runctions. !doesn't work
    use constants
    implicit none
    real*8 func
    external func
    integer, intent(in) :: n,m,i
    real*8,allocatable, dimension(:) :: x, xx, y, yy , test, wtest
    real*8,allocatable :: wx(:), wy(:)
    real*8, intent(out) :: integral
    integer c,j, mtest
    integer ( kind = 4 ) kind
    real ( kind = 8 ),allocatable :: cx(:,:)  !cx array calls from 0 in the second column.

    allocate (x(n),wx(n),xx(n+1), y(m), wy(m), test(n), wtest(n))
    xx(1)=0

    !zeros Lnm:
    call lm_quadrature_rule ( n, m, x, w )
    kind=2
    call cgqf (m, kind, 0.d0, 0.d0, 0d0, 2*pi, y, wy )

    print*,'laguerre abcisas:', x
    print*,'chebychev abcisas:', y

!    kind=5
!    call cgqf (n, kind, 5d0 , 0.d0, 0do, 1d0, test, wtest )
!    print*,'laguerre abcisas test:', test

    integral=0.d0
!    allocate(cx(n,n+1))
!    call Rlm_polynomial ( n, n, m, x, cx )

    do c=1,n
        do j=1,m
            integral= integral + x(c)*func(x(c),y(j))*wy(j)*wx(c)!*cx(c,m)  !le flta los de laguere.
          !  print*, 'x','', 'y','', 'func(x,y)','', 'wy'
           ! print*, x(c), y(j), func(x(c),y(j)), wy(j)
        end do
    end do

print*, 'integral: ', integral
end subroutine
subroutine quadrature_chev1(func,m, integral) !for runctions. !doesn't work
    use constants
    implicit none

    real*8 func
    external func
    integer, intent(in) :: m
    real*8,allocatable, dimension(:) :: y
    real*8,allocatable :: wy(:)
    real*8, intent(out) :: integral
    integer :: j
    integer ( kind = 4 ) :: kind

    allocate (y(m), wy(m))

    kind=2
    call cgqf (m, kind, 0.d0, 0.d0, 0d0, 2*pi, y, wy )
!    print*, ''
!    print*,'chebychev abcisas:', y

    integral=0.d0
        do j=1,m
            integral= integral + func(y(j))*wy(j)
        end do

    print*, ''
    print*, 'integral by chebychev type 1 in [0,2pi]: ', integral
end subroutine
subroutine quadrature_gauchv(func,m, integral) !for runctions. !doesn't work
    use constants
    implicit none

    real*8 func
    external func
    integer, intent(in) :: m
    real*8,allocatable, dimension(:) :: y
    real*8,allocatable :: wy(:)
    real*8, intent(out) :: integral
    integer :: j
    integer ( kind = 4 ) :: kind

    allocate (y(m), wy(m))

    call gauchv(y,wy,m)
!    print*, ''
!    print*,'chebychev abcisas:', y

    integral=0.d0
        do j=1,m
            integral= integral + func(y(j))*wy(j)
        end do

    print*, ''
    print*, 'integral by chebychev type 1 in [-1,1]: ', integral
end subroutine
subroutine quadrature_chev2(func,m, integral) !for runctions. !doesn't work
    use constants
    implicit none

    real*8 func
    external func
    integer, intent(in) :: m
    real*8,allocatable, dimension(:) :: y
    real*8,allocatable :: wy(:)
    real*8, intent(out) :: integral
    integer :: j
    integer ( kind = 4 ) :: kind

    allocate (y(m), wy(m))

    kind=9
    call cgqf (m, kind, 0.d0, 0.d0, 0d0, 2*pi, y, wy )
    !print*, ''
   ! print*,'chebychev abcisas:', y

    integral=0.d0
        do j=1,m
            integral= integral + func(y(j))*wy(j)
        end do

    print*, ''
    print*, 'integral by chebychev type 2 in [0,2pi]: ', integral
end subroutine
double precision function titas(x)
 implicit none
    real*8, intent(in) :: x

    titas=sin(20.*x)*sin(20.*x)
end function
subroutine testing_newtonevol()

    use constants
    use resolution
    implicit none

    real*8 :: c0(l,q,i,9), c1(l,q,i,9)
    real*8, dimension(:), allocatable :: timetest
    integer :: j,x,z,n,s,smax,tsteps
    real*8 :: t0,tmin,tmax
    character(len=20) :: filename

    integer,parameter :: rho_dim=200
    integer,parameter :: phi_dim=30
    real*8 :: rho(rho_dim)
    real*8 :: phi(phi_dim)
    real*8 :: intensity(rho_dim,phi_dim)
    real*8 :: cx(rho_dim ,0:l-1)  !cx array calls from 0 in the second column.
    real*8 :: cy(rho_dim , phi_dim ,0:l-1)  !cx array calls from 0 in the second column.
    real*8 :: dt=0.5

    print*, ''
    print*, 'Time Integration  Test'
    print*, ''

    do j=1,i
        do x=1,q
            do z=1,l
                do n=1,9
                    c0(z,x,j,n)=1
                end do
            end do
        end do
    end do

    OPEN(2,file='initial_coeffs.in',form='unformatted',status='replace',access='stream')
    write(2), c0
    close(2)

    OPEN(2,file='initial_coeffs.in',form='unformatted',status='unknown',access='stream')
        read(2), c0
    close(2)

    tsteps=20
    allocate(timetest(tsteps)) !imprimo 20 timepos
    call linspace(timetest, 0.d0 , 50.d0 , tsteps, 1)

    !*********** Here is  the newton integration***********
    smax=tsteps
    do s=1,smax  !aca va shape(time)
        !set integration times
        tmin=timetest(s)
        tmax=timetest(s+1)

        t0=tmin
        do while (t0<=tmax)  ! Stop if we've exceeded TMAX.

            print*, 't=', t0
            !set initial coefs
            OPEN(2,file='initial_coeffs.in',form='unformatted',status='unknown',access='stream')
                read(2), c0
            close(2)

            call newton(c0,c1)
            c0=c1  !set new c0
            t0=t0 + dt  !dt !set next step of integration

        enddo

        OPEN(2,file='initial_coeffs.in',form='unformatted',status='replace',access='stream')
            write(2), c1
        close(2)

        !Print coefs values at time(s+1)
        write (filename, '("coeff_time_",i4.4,".in")') , s
        OPEN(2,file='filename',form='unformatted',status='replace',access='stream')
            write(2), c1
        close(2)
        print*, 'saved ', filename

    end do

    call linspace(rho, 0.d0 , 50.d0 , rho_dim, 1)
    call linspace(phi, 0.d0 , 2*pi , phi_dim, 1)

    OPEN(2,file='rho.in',form='unformatted',status='replace',access='stream')
    write(2), rho
    close(2)
    OPEN(2,file='phi.in',form='unformatted',status='replace',access='stream')
    write(2), phi
    close(2)

    call lm_polynomial ( rho_dim, l, 1, rho, cx )
    !call Rlm_polynomial ( rho_dim, l, 0, rho, cx )
    !print*, cx
    OPEN(2,file='lm.in',form='unformatted',status='replace',access='stream')
    write(2), cx
    close(2)
    print*,cx(:,0)
    print*, 'hasta aca llegamo-'

    !call Alm_polynomial ( rho_dim, phi_dim, l, q, i, rho, phi, cx, cy )

end subroutine

subroutine testing_rk_evol()

    use constants
    use resolution
    implicit none

    real*8,dimension(l,q,i,5) :: c0, c1
    real*8, dimension(:), allocatable :: timetest
    integer :: j,x,z,n,s,smax,tsteps
    real*8 :: t0,tmin,tmax
    character(len=20) :: filename
    real*8,dimension(l,q,i) ::  phi_r,phi_i,pol_r,pol_i,pop
    real*8,dimension(l,q,i) ::  new_phi_r,new_phi_i,new_pol_r,new_pol_i,new_pop

    integer,parameter :: rho_dim=1000
    integer,parameter :: phi_dim=500
    real*8 :: rho(rho_dim)
    real*8 :: phi(phi_dim)
    real*8,dimension(rho_dim,phi_dim) :: intensity, er, ei, population
    real ( kind = 8 ) cx(rho_dim ,0:l-1)  !cx array calls from 0 in the second column.
    real ( kind = 8 ) cy(rho_dim , phi_dim ,0:l-1)  !cx array calls from 0 in the second column.
    real ( kind = 8 ) :: dt=0.5

    print*, ''
    print*, 'Time Integration  Test'
    print*, ''
    call save_resolution() !save resolucion values to txt, to pass to python plot.
    !******************set initial field values.***************************
    phi_r(:,:,:)=1
    phi_i(:,:,:)=1
    pol_r(:,:,:)=1
    pol_i(:,:,:)=1
    pop(:,:,:)=1
    !***************** Print initial values to .bin************************
    call save_fields()
    !*****************************Read initial values to .bin************************
    call read_fields()

    tsteps=10
    allocate(timetest(tsteps)) !imprimo 20 timepos
    call linspace(timetest, 0.d0 , 100.d0 , tsteps, 1)

    !*********** Here is  the newton integration***********
    smax=tsteps
    do s=1,smax  !aca va shape(time)
        !set integration times
        tmin=timetest(s)
        tmax=timetest(s+1)

        t0=tmin
        do while (t0<=tmax)  ! Stop if we've exceeded TMAX.
            print*, 't=', t0
            !set initial coefs
            call read_fields()

            c0(:,:,:,1)=phi_r
            c0(:,:,:,2)=phi_i
            c0(:,:,:,3)=pol_r
            c0(:,:,:,4)=pol_i
            c0(:,:,:,5)=pop

            call newton(c0,c1)

            new_phi_r=c1(:,:,:,1)
            new_phi_i=c1(:,:,:,2)
            new_pol_r=c1(:,:,:,3)
            new_pol_i=c1(:,:,:,4)
            new_pop=c1(:,:,:,5)
            phi_r=new_phi_r
            phi_i=new_phi_i
            pol_r=new_pol_r
            pol_i=new_pol_i
            pop=new_pop

            !c0=c1  !set new c0
            t0=t0 + dt  !dt !set next step of integration
        enddo

        call save_new_fields

        OPEN(2,file='initial_coeffs.in',form='unformatted',status='replace',access='stream')
            write(2), c1
        close(2)
       !Print coefs values at time(s+1)
        write (filename, '("coeff_time_",i4.4,".in")') , s
        OPEN(2,file='filename',form='unformatted',status='replace',access='stream')
            write(2), c1
        close(2)
        print*, 'saved ', filename

    end do

    call linspace(rho, 0.d0 , 8.d0 , rho_dim, 1)
    call linspace(phi, 0.d0 , 2*pi , phi_dim, 1)
    OPEN(2,file='rho.in',form='unformatted',status='replace',access='stream')
        write(2), rho
    close(2)
    OPEN(2,file='phi.in',form='unformatted',status='replace',access='stream')
        write(2), phi
    close(2)

    intensity(:,:)=0

    do x=0,q
        do j=1,i
            call Alm_polynomial ( rho_dim, phi_dim, l, x, i, rho, phi, cx, cy )
            OPEN(2,file='Alm.in',form='unformatted',status='replace',access='stream')
                write(2), cy
            close(2)
            do z=1,l
!                print*, phi_r(z,x,j)
                er(:,:)=er(:,:)+(phi_r(z,x,j)*cy(:,:,z-1))
                ei(:,:)=ei(:,:)+(phi_r(z,x,j)*cy(:,:,z-1))
                population(:,:)=population(:,:)+(pop(z,x,j)*cy(:,:,z-1))
                intensity=sqrt(er**2+ei**2)
            end do
        end do
    end do

    OPEN(2,file='intensity.in',form='unformatted',status='replace',access='stream')
    write(2), intensity
    close(2)
    OPEN(2,file='population.in',form='unformatted',status='replace',access='stream')
    write(2), population
    close(2)

    print*, 'hasta aca llegamo-'


end subroutine
subroutine rieman_integration_BONtest()
    use constants
    use funcs
    use resolution

    implicit none

    integer,parameter :: rho_dim=2000 , phi_dim=2000
    real*8 :: rho(rho_dim), phi(phi_dim)
    real ( kind = 8 ), dimension(rho_dim ,0:l) :: cx1, cx2  !cx array calls from 0 in the second column.
    real ( kind = 8 ), dimension(rho_dim , phi_dim ,0:l) :: cy1, cy2  !cx array calls from 0 in the second column.
    real*8 :: sumation
    integer ( kind = 4 ) c,c2
    real*8 integ
    real*8 :: coefs(0:l,0:q,i)
    integer :: j,x,z,n,s,j2,x2,z2
   !********************************************************************************
    call linspace(rho, 0.d0 , 40.d0 , rho_dim, 1)
    call linspace(phi, 0.d0 , 2*pi , phi_dim, 1)
    open (unit=2,file='BON_check.txt',action="write",status='replace')
    write (2,'("integration with")')
do z=1,i
    do x=0,q
        do z2=1,i
            do x2=0,q
                call Alm_polynomial( rho_dim, phi_dim, l, x, z, rho, phi, cx1, cy1 )
                call Alm_polynomial( rho_dim, phi_dim, l, x2, z2, rho, phi, cx2, cy2 )
                do j=0,l
                    do j2=0,l
                        sumation=0
                        do c2=1,phi_dim-1
                            do c=1,rho_dim-1
                                sumation=sumation + cy1(c,c2,j)*cy2(c,c2,j2)*(rho(c+1)-rho(c))*(phi(c2+1)-phi(c2))*rho(c)
                            end do
                        end do
                        write (2,'("n= ", i4," m= ", i4," i= ", i4,  /, "n*=", i4," m*=",i4," i*=", i4)') , j,x,z,j2,x2,z2
                        write (2,*) sumation
                    end do
                end do
            end do
        end do
    end do
end do
close(2)
end subroutine

subroutine newton(c0,c1)
    use constants
    use resolution
    implicit none

    real*8, intent(in) :: c0(l,q,i,5)
    real*8, intent(out) ::  c1(l,q,i,5)
    real*8 ::  dt=0.005d0
    integer :: j,x,z

    do j=1,i
        do x=1,q
            do z=1,l
                c1(z,x,j,1)=c0(z,x,j,1)+dt*(-k*c0(z,x,j,1)+(-a+d)*c0(z,x,j,2)-2*k*c0(z,x,j,3))
                c1(z,x,j,2)=c0(z,x,j,1)+dt*(-k*c0(z,x,j,2)+(a-d)*c0(z,x,j,1)-2*k*c0(z,x,j,4))
                c1(z,x,j,3)=c0(z,x,j,1)+dt*(c0(z,x,j,3)-d*c0(z,x,j,4))
                c1(z,x,j,4)=c0(z,x,j,1)+dt*(c0(z,x,j,4)+d*c0(z,x,j,3))
                c1(z,x,j,5)=c0(z,x,j,1)+dt*(-g*(-c0(z,x,j,5)))
            end do
        end do
    end do
end subroutine

subroutine runges(c0,c1)
    use constants
    use resolution
    use rungekutta
    implicit none

    real*8, intent(in) :: c0(l,q,i,5)
    real*8, intent(out) :: c1(l,q,i,5)
    real*8 ::  dt=0.005d0
    integer :: j,x,z,order
    real*8 :: rkcoef


    do order=ord,1,-1
        rkcoef=1./order
        do j=1,i
            do x=1,q
                do z=1,l
                    c1(z,x,j,1)=c0(z,x,j,1)+dt*(-k*c0(z,x,j,1)+(-a+d)*c0(z,x,j,2)-2*k*c0(z,x,j,3))*rkcoef
                    c1(z,x,j,2)=c0(z,x,j,1)+dt*(-k*c0(z,x,j,2)+(a-d)*c0(z,x,j,1)-2*k*c0(z,x,j,4))*rkcoef
                    c1(z,x,j,4)=c0(z,x,j,1)+dt*(c0(z,x,j,3)-d*c0(z,x,j,4))*rkcoef
                    c1(z,x,j,4)=c0(z,x,j,1)+dt*(c0(z,x,j,4)+d*c0(z,x,j,3))*rkcoef
                    c1(z,x,j,5)=c0(z,x,j,1)+dt*(-g*(-c0(z,x,j,5)))*rkcoef
                end do
            end do
        end do
    end do
end subroutine

subroutine linspace(x,x_start, x_end, x_len,dir)
!******************************************************************************
!linearly spaced array named x, from x_start value to x_end, with #x_len values.
!dir=1 ---> from min to max.    dir=2 ---> from max to min.
!******************************************************************************
    implicit none
    real(8), dimension(x_len), intent(inout) :: x
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
subroutine Alm_polynomial ( mm, ii, n, m, i, x, y, cx, cy ) !Anm=2*((2x**2)**(m/2))*((facoriales)**(1/2))*exp(-x**2)*Lnm(2x**2)*exp(j*m*phi) gives 0 to n Ajm(rho,phi)

!*****************************************************************************
!  Special values:
!
!    For M = 0, the associated Laguerre polynomials Lm(N,M,X) are equal
!    to the Laguerre polynomials L(N,X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22/3/2016
!
!  Author:
!
!    John Burkardt mod Alexis Gomel
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, integer ( kind = 4 ) M, the parameter.  M must be nonnegative.
!
!    Input, real ( kind = 8 ) X(MM), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the associated Laguerre polynomials
!    of degrees 0 through N evaluated at the evaluation points.
!

    use constants

    implicit none

    integer ( kind = 4 ),intent(in) :: mm,ii   !steps in rho,steps in phi dir
    integer ( kind = 4 ),intent(in) :: n,m,i
    real ( kind = 8 ),intent(in) :: cx(mm ,0:n)  !cx array calls from 0 in the second column.
    real ( kind = 8 ),intent(out) ::  cy(mm , ii ,0:n)  !cx array calls from 0 in the second column.
    integer ( kind = 4 ) q,j,s
    real ( kind = 8 ),intent(in) ::  x(mm)
    real ( kind = 8 ),intent(in) ::  y(ii)
    real*8 :: norm

    call  Rlm_polynomial ( mm, n, m, x, cx )

    if (m.eq.0) then
        norm=sqrt(1/(2*pi))
        do j=0,n !take a look at how i implemet the do.
            do s=1,mm
                do q=1,ii
                    cy(s,q,j)= norm*cx(s,j)
                end do
            end do
        end do
    elseif (i.eq.1) then
        norm=sqrt(1/(pi))
        do j=0,n !take a look at how i implemet the do.
            do s=1,mm
                do q=1,ii
                    cy(s,q,j)= norm*cx(s,j)*sin(m*y(q))
                end do
            end do
        end do
    elseif (i.eq.2) then
        norm=sqrt(1/(pi))
        do j=0,n !take a look at how i implemet the do.
            do s=1,mm
                do q=1,ii
                    cy(s,q,j)= norm*cx(s,j)*cos(m*y(q))
                end do
            end do
        end do
    elseif (i>3) then
        print*, 'Error: i > 2'
        stop
    elseif (i<0) then
        print*, 'Error: i < 2'
        stop
    end if
    return
end
subroutine Rlm_polynomial ( mm, n, m, x, cx ) !Rnm=2*((2x**2)**(m/2))*((facoriales)**(1/2))*exp(-x**2)*Lnm(2x**2) gives 0 to n Rjm(x)

!*****************************************************************************
!  Special values:
!
!    For M = 0, the associated Laguerre polynomials Lm(N,M,X) are equal
!    to the Laguerre polynomials L(N,X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22/3/2016
!
!  Author:
!
!    John Burkardt mod Alexis Gomel
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, integer ( kind = 4 ) M, the parameter.  M must be nonnegative.
!
!    Input, real ( kind = 8 ) X(MM), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the associated Laguerre polynomials
!    of degrees 0 through N evaluated at the evaluation points.
!
    implicit none

    integer ( kind = 4 ) mm
    integer ( kind = 4 ),intent(in) :: n,m
    real ( kind = 8 ) cx(mm,0:n)  !cx array calls from 0 in the second column.
    integer ( kind = 4 ) i,j
    real ( kind = 8 ) x(mm)
    real ( kind = 8 ) y(mm)
    real*8 div

    y=x**2
    call  lm_polynomial ( mm, n, m, 2*y,cx )
    do j=0,n !take a look at how i implemet the do.
        div=1./gamma(j+m+1.)
        do i=1,mm
            cx(i,j)= 2*((2*y(i))**(m/2))*sqrt(gamma(j+1.)*div)*cx(i,j)*exp(-y(i))
        end do
    end do
    !  return
end subroutine
subroutine save_fields()

    use resolution
    implicit none
    real*8 ::   phi_r(l,q,i),phi_i(l,q,i),pol_r(l,q,i),pol_i(l,q,i),pop(l,q,i)

        OPEN(2,file='initial_phir.in',form='unformatted',status='replace',access='stream')
            write(2), phi_r
        close(2)
        OPEN(2,file='initial_phii.in',form='unformatted',status='replace',access='stream')
            write(2), phi_i
        close(2)
        OPEN(2,file='initial_polr.in',form='unformatted',status='replace',access='stream')
            write(2), pol_r
        close(2)
        OPEN(2,file='initial_poli.in',form='unformatted',status='replace',access='stream')
            write(2), pol_i
        close(2)
        OPEN(2,file='initial_pop.in',form='unformatted',status='replace',access='stream')
            write(2), pop
        close(2)
    end subroutine
subroutine save_new_fields()

    use resolution
    implicit none
    real*8 ::   new_phi_r(l,q,i),new_phi_i(l,q,i),new_pol_r(l,q,i),new_pol_i(l,q,i),new_pop(l,q,i)
        OPEN(2,file='initial_phir.in',form='unformatted',status='replace',access='stream')
            write(2), new_phi_r
        close(2)
        OPEN(2,file='initial_phii.in',form='unformatted',status='replace',access='stream')
            write(2), new_phi_i
        close(2)
        OPEN(2,file='initial_polr.in',form='unformatted',status='replace',access='stream')
            write(2), new_pol_r
        close(2)
        OPEN(2,file='initial_poli.in',form='unformatted',status='replace',access='stream')
            write(2), new_pol_i
        close(2)
        OPEN(2,file='initial_pop.in',form='unformatted',status='replace',access='stream')
            write(2), new_pop
        close(2)
    end subroutine
subroutine read_fields()

      use resolution
      implicit none
      real*8 ::   phi_r(l,q,i),phi_i(l,q,i),pol_r(l,q,i),pol_i(l,q,i),pop(l,q,i)
        OPEN(2,file='initial_phir.in',form='unformatted',status='unknown',access='stream')
            read(2), phi_r
        close(2)
        OPEN(2,file='initial_phii.in',form='unformatted',status='unknown',access='stream')
            read(2), phi_i
        close(2)
        OPEN(2,file='initial_polr.in',form='unformatted',status='unknown',access='stream')
            read(2), pol_r
        close(2)
        OPEN(2,file='initial_poli.in',form='unformatted',status='unknown',access='stream')
            read(2), pol_i
        close(2)
        OPEN(2,file='initial_pop.in',form='unformatted',status='unknown',access='stream')
            read(2), pop
        close(2)
    end subroutine
double precision function pump(x,y)
    implicit none
    real*8, intent(in) :: x,y
    real*8, parameter :: rho_0=1.d0

    pump=.5d0-tanh(5d0*(x-rho_0))/2d0
end function

!*********************** Rountines i need to do*****************
subroutine proyection(func,n,m,i, integral)  !function proyection in Almi

    !call Almi_policnimial(cx)
    !call quadrature(function,x,cx)
     use constants
    implicit none
    real*8 func
    external func
    integer, intent(in) :: n,m,i
    real*8,allocatable, dimension(:) :: x, xx, y, yy , test, wtest
    real*8,allocatable :: wx(:), wy(:)
    real*8, intent(out) :: integral
    integer c,j, mtest
    integer ( kind = 4 ) kind
    real ( kind = 8 ),allocatable :: cx(:,:)  !cx array calls from 0 in the second column.

    character(len=4) :: ext
    character(len=6), save :: fmtext = '(i4.4)'!(integer 4 digitos.)

    allocate (x(n),wx(n),xx(n+1), y(m), wy(m), test(n), wtest(n))
    xx(1)=0

    !zeros Lnm:
    call lm_quadrature_rule ( n, m, x, w )
    kind=2
    call cgqf (m, kind, 0.d0, 0.d0, 0d0, 2*pi, y, wy )

    print*,'laguerre abcisas:', x
    print*,'chebychev abcisas:', y

!    kind=5
!    call cgqf (n, kind, 5d0 , 0.d0, 0do, 1d0, test, wtest )
!    print*,'laguerre abcisas test:', test

    integral=0.d0
!    allocate(cx(n,n+1))
!    call Rlm_polynomial ( n, n, m, x, cx )

    do c=1,n
        do j=1,m
            integral= integral + x(c)*func(x(c),y(j))*wy(j)*wx(c)*cx(c,m)  !le flta los de laguere.
!              IF (i.eq.0) then
!                  do j=0,n !take a look at how i implemet the do.
!                    do i=1,mm
!                        do q=1,ii
!                            cy(i,q,j)= 1/(2*pi)*cx(i,j)
!                        end do
!                    end do
!                  end do
!              elseif (i.eq.1) then
!                  do j=0,n !take a look at how i implemet the do.
!                    do i=1,mm
!                        do q=1,ii
!                            cy(i,q,j)= (1/(pi))*cx(i,j)*sin(m*y(q))
!                        end do
!                    end do
!                  end do
!              elseif (i.eq.2) then
!                  do j=0,n !take a look at how i implemet the do.
!                      do i=1,mm
!                          do q=1,ii
!                              cy(i,q,j)= (1/(pi))*cx(i,j)*cos(m*y(q))
!                          end do
!                      end do
!                  end do
!              elseif (i>3) then
!                print*, 'Error: i > 2'
!                stop
!            elseif (i<0) then
!                print*, 'Error: i < 2'
!                stop
!            END if
           !  print*, 'x','', 'y','', 'func(x,y)','', 'wy'
           ! print*, x(c), y(j), func(x(c),y(j)), wy(j)
        end do
    end do

    write(ext, fmtext) n, m, i
  !  file=trim(dir) // '/' // 'coefs'// '.' // ext // '.bin'

end subroutine
real function z(p,m,i)
integer :: p,m,i
!call nonlinear_coupling(p,m,i,z)
z(p,m,i)=1
end function
subroutine nonlinear_coupling() !integration of the coupling between the coefficients.
    !input <-- pm
    !output --> modemode(p,m,i) Table.
    !Modemode

    !call quadrature(ro*Apmi*Appmmii*Apppmmiii,modemode)

    !modemode= sum(pp,mm,ii,ppp,mmm,iii)

    !save
end subroutine

!****************************** gen_laguerre(kind 5) & chebychev(kind 2) subruotines*******************************************
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: itn = 30
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) prec
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end
subroutine gen_laguerre_test()
!*****************************************************************************80
!
!! MAIN is the main program for GEN_LAGUERRE_RULE.
!
!  Discussion:
!
!    This program computes a generalized Gauss-Laguerre quadrature rule
!    and writes it to a file.
!
!    The user specifies:
!    * the ORDER (number of points) in the rule;
!    * ALPHA, the exponent of |X|;
!    * A, the left endpoint;
!    * B, the scale factor;
!    * FILENAME, the root name of the output files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  integer ( kind = 4 ) arg_num
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  character ( len = 255 )filename
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) last
  integer ( kind = 4 ) order
  real ( kind = 8 ) r(2)
  real ( kind = 8 ) r8_huge
  character ( len = 255 ) string
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  !call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEN_LAGUERRE_RULE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Compute a generalized Gauss-Laguerre rule for approximating'
  write ( *, '(a)' ) &
    '    Integral ( a <= x < +oo ) |x-a|^ALPHA exp(-B*(x-a)) f(x) dx'
  write ( *, '(a)' ) '  of order ORDER.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The user specifies ORDER, ALPHA, A, B, and FILENAME.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ORDER is the number of points.'
  write ( *, '(a)' ) '  ALPHA is the exponent of |X|.'
  write ( *, '(a)' ) '  A is the left endpoint, (typically 0).'
  write ( *, '(a)' ) '  B is the scale factor, (typically 1).'
  write ( *, '(a)' ) '  FILENAME is used to generate 3 files:'
  write ( *, '(a)' ) '  * filename_w.txt - the weight file'
  write ( *, '(a)' ) '  * filename_x.txt - the abscissa file.'
  write ( *, '(a)' ) '  * filename_r.txt - the region file.'
!
!  Initialize parameters.
!
  beta = 0.0D+00
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get ORDER.
!
  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, order, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the rule order ORDER:'
    read ( *, * ) order
  end if
!
!  Get ALPHA.
!
  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, string )
    call s_to_r8 ( string, alpha, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the value of the parameter ALPHA:'
    write ( *, '(a)' ) '  (requirement: -1.0 < ALPHA )'
    read ( *, * ) alpha
  end if
!
!  Get A.
!
  if ( 3 <= arg_num ) then
    iarg = 3
    call getarg ( iarg, string )
    call s_to_r8 ( string, a, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the left endpoint, A:'
    read ( *, * ) a
  end if
!
!  Get B.
!
  if ( 4 <= arg_num ) then
    iarg = 4
    call getarg ( iarg, string )
    call s_to_r8 ( string, b, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the scale factor B:'
    read ( *, * ) b
  end if
!
!  Get FILENAME.
!
  if ( 5 <= arg_num ) then
    iarg = 5
    call getarg ( iarg, filename )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the "root name" of the quadrature files).'
    read ( *, '(a)' ) filename
  end if
!
!  Input summary.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ORDER = ', order
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a)' ) '  FILENAME = "' // trim ( filename ) // '".'
!
!  Construct the rule.
!
  allocate ( w(order) )
  allocate ( x(order) )

  kind = 5
  call cgqf ( order, kind, alpha, beta, a, b, x, w )
!
!  Write the rule.
!
  r(1) = a
  r(2) = r8_huge ( )

  call rule_write ( order, x, w, r, filename )
!
!  Free memory.
!
  deallocate ( w )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEN_LAGUERRE_RULE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  !call timestamp ( )



  stop
end subroutine
subroutine cdgqf ( nt, kind, alpha, beta, t, wts )

!*****************************************************************************80
!
!! CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with a classical weight function with default values for A and B,
!    and only simple knots.
!
!    There are no moments checks and no printing is done.
!
!    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) aj(nt)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) kind
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) zemu

  call parchk ( kind, 2 * nt, alpha, beta )
!
!  Get the Jacobi matrix and zero-th moment.
!
  call class_matrix ( kind, nt, alpha, beta, aj, bj, zemu )
!
!  Compute the knots and weights.
!
  call sgqf ( nt, aj, bj, zemu, t, wts )

  return
end
subroutine cgqf ( nt, kind, alpha, beta, a, b, t, wts )

!*****************************************************************************80
!
!! CGQF computes knots and weights of a Gauss quadrature formula.
!
!  Discussion:
!
!    The user may specify the interval (A,B).
!
!    Only simple knots are produced.
!
!    Use routine EIQFS to evaluate this quadrature formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints, or
!    other parameters.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ), allocatable :: ndx(:)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
!
!  Compute the Gauss quadrature formula for default values of A and B.
!
  call cdgqf ( nt, kind, alpha, beta, t, wts )
!
!  Prepare to scale the quadrature formula to other weight function with
!  valid A and B.
!
  allocate ( mlt(1:nt) )

  mlt(1:nt) = 1

  allocate ( ndx(1:nt) )

  do i = 1, nt
    ndx(i) = i
  end do

  call scqf ( nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b )

  deallocate ( mlt )
  deallocate ( ndx )

  return
end
subroutine scqf ( nt, t, mlt, wts, nwts, ndx, swts, st, kind, alpha, beta, a, &
  b )

!*****************************************************************************80
!
!! SCQF scales a quadrature formula to a nonstandard interval.
!
!  Discussion:
!
!    The arrays WTS and SWTS may coincide.
!
!    The arrays T and ST may coincide.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the original knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, real ( kind = 8 ) WTS(NWTS), the weights.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input, integer ( kind = 4 ) NDX(NT), used to index the array WTS.
!    For more details see the comments in CAWIQ.
!
!    Output, real ( kind = 8 ) SWTS(NWTS), the scaled weights.
!
!    Output, real ( kind = 8 ) ST(NT), the scaled knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( kind = 8 ) a
  real ( kind = 8 ) al
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) be
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) ndx(nt)
  real ( kind = 8 ) p
  real ( kind = 8 ) shft
  real ( kind = 8 ) slp
  real ( kind = 8 ) st(nt)
  real ( kind = 8 ) swts(nwts)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) temp
  real ( kind = 8 ) tmp
  real ( kind = 8 ) wts(nwts)

  temp = epsilon ( temp )

  call parchk ( kind, 1, alpha, beta )

  if ( kind == 1 ) then

    al = 0.0D+00
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 2 ) then

    al = -0.5D+00
    be = -0.5D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 3 ) then

    al = alpha
    be = alpha

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 4 ) then

    al = alpha
    be = beta

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 5 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  B <= 0'
      stop
    end if

    shft = a
    slp = 1.0D+00 / b
    al = alpha
    be = 0.0D+00

  else if ( kind == 6 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  B <= 0.'
      stop
    end if

    shft = a
    slp = 1.0D+00 / sqrt ( b )
    al = alpha
    be = 0.0D+00

  else if ( kind == 7 ) then

    al = alpha
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 8 ) then

    if ( a + b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  A + B <= 0.'
      stop
    end if

    shft = a
    slp = a + b
    al = alpha
    be = beta

  else if ( kind == 9 ) then

    al = 0.5D+00
    be = 0.5D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  end if

  p = slp**( al + be + 1.0D+00 )

  do k = 1, nt

    st(k) = shft + slp * t(k)
    l = abs ( ndx(k) )

    if ( l /= 0 ) then
      tmp = p
      do i = l, l + mlt(k) - 1
        swts(i) = wts(i) * tmp
        tmp = tmp * slp
      end do
    end if

  end do

  return
end
subroutine sgqf ( nt, aj, bj, zemu, t, wts )

!*****************************************************************************80
!
!! SGQF computes knots and weights of a Gauss Quadrature formula.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with simple knots from the Jacobi matrix and the zero-th
!    moment of the weight function, using the Golub-Welsch technique.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) AJ(NT), the diagonal of the Jacobi matrix.
!
!    Input/output, real ( kind = 8 ) BJ(NT), the subdiagonal of the Jacobi
!    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
!
!    Input, real ( kind = 8 ) ZEMU, the zero-th moment of the weight function.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) aj(nt)
  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) zemu
!
!  Exit if the zero-th moment is not positive.
!
  if ( zemu <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGQF - Fatal error!'
    write ( *, '(a)' ) '  ZEMU <= 0.'
    stop
  end if
!
!  Set up vectors for IMTQLX.
!
  t(1:nt) = aj(1:nt)

  wts(1) = sqrt ( zemu )
  wts(2:nt) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt)**2

  return
end
subroutine parchk ( kind, m, alpha, beta )

!*****************************************************************************80
!
!! PARCHK checks parameters ALPHA and BETA for classical weight functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!
!    Input, integer ( kind = 4 ) M, the order of the highest moment to
!    be calculated.  This value is only needed when KIND = 8.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters, if required
!    by the value of KIND.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) m
  real ( kind = 8 ) tmp

  if ( kind <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  KIND <= 0.'
    stop
  end if
!
!  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
!
  if ( 3 <= kind .and. alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  3 <= KIND and ALPHA <= -1.'
    stop
  end if
!
!  Check BETA for Jacobi.
!
  if ( kind == 4 .and. beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  KIND == 4 and BETA <= -1.0.'
    stop
  end if
!
!  Check ALPHA and BETA for rational.
!
  if ( kind == 8 ) then
    tmp = alpha + beta + m + 1.0D+00
    if ( 0.0D+00 <= tmp .or. tmp <= beta ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PARCHK - Fatal error!'
      write ( *, '(a)' ) '  KIND == 8 but condition on ALPHA and BETA fails.'
      stop
    end if
  end if

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine class_matrix ( kind, m, alpha, beta, aj, bj, zemu )

!*****************************************************************************80
!
!! CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
!
!  Discussion:
!
!    This routine computes the diagonal AJ and sub-diagonal BJ
!    elements of the order M tridiagonal symmetric Jacobi matrix
!    associated with the polynomials orthogonal with respect to
!    the weight function specified by KIND.
!
!    For weight functions 1-7, M elements are defined in BJ even
!    though only M-1 are needed.  For weight function 8, BJ(M) is
!    set to zero.
!
!    The zero-th moment of the weight function is returned in ZEMU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!
!    Input, integer ( kind = 4 ) M, the order of the Jacobi matrix.
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) AJ(M), BJ(M), the diagonal and subdiagonal
!    of the Jacobi matrix.
!
!    Output, real ( kind = 8 ) ZEMU, the zero-th moment.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a2b2
  real ( kind = 8 ) ab
  real ( kind = 8 ) aba
  real ( kind = 8 ) abi
  real ( kind = 8 ) abj
  real ( kind = 8 ) abti
  real ( kind = 8 ) aj(m)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) apone
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp2
  real ( kind = 8 ) zemu

  temp = epsilon ( temp )

  call parchk ( kind, 2 * m - 1, alpha, beta )

  temp2 = 0.5D+00

  if ( 500.0D+00 * temp < abs ( ( r8_gamma ( temp2 ) )**2 - pi ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLASS_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  Gamma function does not match machine parameters.'
    stop
  end if

  if ( kind == 1 ) then

    ab = 0.0D+00

    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 2 ) then

    zemu = pi

    aj(1:m) = 0.0D+00

    bj(1) =  sqrt ( 0.5D+00 )
    bj(2:m) = 0.5D+00

  else if ( kind == 3 ) then

    ab = alpha * 2.0D+00
    zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 )**2 &
      / r8_gamma ( ab + 2.0D+00 )

    aj(1:m) = 0.0D+00
    bj(1) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
    do i = 2, m
      bj(i) = i * ( i + ab ) / ( 4.0D+00 * ( i + alpha )**2 - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 4 ) then

    ab = alpha + beta
    abi = 2.0D+00 + ab
    zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 ) &
      * r8_gamma ( beta + 1.0D+00 ) / r8_gamma ( abi )
    aj(1) = ( beta - alpha ) / abi
    bj(1) = 4.0D+00 * ( 1.0 + alpha ) * ( 1.0D+00 + beta ) &
      / ( ( abi + 1.0D+00 ) * abi * abi )
    a2b2 = beta * beta - alpha * alpha

    do i = 2, m
      abi = 2.0D+00 * i + ab
      aj(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
      abi = abi**2
      bj(i) = 4.0D+00 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) &
        / ( ( abi - 1.0D+00 ) * abi )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 5 ) then

    zemu = r8_gamma ( alpha + 1.0D+00 )

    do i = 1, m
      aj(i) = 2.0D+00 * i - 1.0D+00 + alpha
      bj(i) = i * ( i + alpha )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 6 ) then

    zemu = r8_gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      bj(i) = ( i + alpha * mod ( i, 2 ) ) / 2.0D+00
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 7 ) then

    ab = alpha
    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 8 ) then

    ab = alpha + beta
    zemu = r8_gamma ( alpha + 1.0D+00 ) * r8_gamma ( - ( ab + 1.0D+00 ) ) &
      / r8_gamma ( - beta )
    apone = alpha + 1.0D+00
    aba = ab * apone
    aj(1) = - apone / ( ab + 2.0D+00 )
    bj(1) = - aj(1) * ( beta + 1.0D+00 ) / ( ab + 2.0D+00 ) / ( ab + 3.0D+00 )
    do i = 2, m
      abti = ab + 2.0D+00 * i
      aj(i) = aba + 2.0D+00 * ( ab + i ) * ( i - 1 )
      aj(i) = - aj(i) / abti / ( abti - 2.0D+00 )
    end do

    do i = 2, m - 1
      abti = ab + 2.0D+00 * i
      bj(i) = i * ( alpha + i ) / ( abti - 1.0D+00 ) * ( beta + i ) &
        / ( abti**2 ) * ( ab + i ) / ( abti + 1.0D+00 )
    end do

    bj(m) = 0.0D+00
    bj(1:m) =  sqrt ( bj(1:m) )

  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
function r8_epsilon ( )

!*****************************************************************************80
!
!! R8_EPSILON returns the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_EPSILON, the round-off unit.
!
  implicit none

  real ( kind = 8 ) r8_epsilon

  r8_epsilon = 2.220446049250313D-016

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none

  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = 8 ), parameter :: eps = 2.22D-16
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 171.624D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.0D+30
  real ( kind = 8 ), parameter :: xminin = 2.23D-308
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real ( kind = 8 ) r8_huge

  r8_huge = 1.0D+30

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine rule_write ( order, x, w, r, filename )

!*****************************************************************************80
!
!! RULE_WRITE writes a quadrature rule to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) X(ORDER), the abscissas.
!
!    Input, real ( kind = 8 ) W(ORDER), the weights.
!
!    Input, real ( kind = 8 ) R(2), defines the region.
!
!    Input, character ( len = * ) FILENAME, specifies the output.
!    'filename_w.txt', 'filename_x.txt', 'filename_r.txt' defining weights,
!    abscissas, and region.
!
  implicit none

  integer ( kind = 4 ) order

  character ( len = * ) filename
  character ( len = 255 ) filename_r
  character ( len = 255 ) filename_w
  character ( len = 255 ) filename_x
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  real ( kind = 8 ) r(2)
  real ( kind = 8 ) w(order)
  real ( kind = 8 ) x(order)

  filename_w = trim ( filename ) // '_w.txt'
  filename_x = trim ( filename ) // '_x.txt'
  filename_r = trim ( filename ) // '_r.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Creating quadrature files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "Root" file name is   "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weight file will be   "' // trim ( filename_w ) // '".'
  write ( *, '(a)' ) '  Abscissa file will be "' // trim ( filename_x ) // '".'
  write ( *, '(a)' ) '  Region file will be   "' // trim ( filename_r ) // '".'

  call r8mat_write ( filename_w, 1, order, w )
  call r8mat_write ( filename_x, 1, order, x )
  call r8mat_write ( filename_r, 1, 2,     r )

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S
!    used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) length
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length+1 == nchar ) then
    length = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

  return
end

!*************************************** laguerre_pols subroutines**********************************
subroutine l_exponential_product ( p, b, table )

!*****************************************************************************80
!
!! L_EXPONENTIAL_PRODUCT: exponential product table for L(n,x).
!
!  Discussion:
!
!    Let L(n,x) represent the Laguerre polynomial of degree n.
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of exp(B*X) with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( 0 <= X < +oo ) exp(b*x) * L(i,x) * L(j,x) * exp (-x) dx
!
!    Because of the exponential factor, the quadrature will not be exact.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial
!    factors.  0 <= P.
!
!    Input, real ( kind = 8 ) B, the coefficient of X in the exponential factor.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.
!    TABLE(I,J) represents the weighted integral of exp(B*X) * L(i,x) * L(j,x).
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), allocatable :: l_table(:)
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = ( 3 * p + 4 ) / 2

  allocate ( x_table(1:order) )
  allocate ( w_table(1:order) )

  call l_quadrature_rule ( order, x_table, w_table )

  allocate ( l_table(0:p) )

  do k = 1, order

    x = x_table(k)
    call l_polynomial ( 1, p, x, l_table )

    do j = 0, p
      do i = 0, p
        table(i,j) = table(i,j) &
          + w_table(k) * exp ( b * x ) * l_table(i) * l_table(j)
      end do
    end do

  end do

  deallocate ( l_table )
  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine l_integral ( n, exact )

!*****************************************************************************80
!
!! L_INTEGRAL evaluates a monomial integral associated with L(n,x).
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the exponent.
!    0 <= N.
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial

  exact = r8_factorial ( n )

  return
end
subroutine l_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! L_POLYNOMIAL evaluates the Laguerre polynomial L(n,x).
!
!  First terms:
!
!      1
!     -X     +  1
!   (  X^2 -  4 X      +  2 ) / 2
!   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
!   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
!   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -   600 X    +  120 ) / 120
!   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 +  5400 X^2 -  4320 X     + 720 )
!     / 720
!   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3 + 52920 X^2 - 35280 X
!     + 5040 ) / 5040
!
!  Recursion:
!
!    L(0,X) = 1
!    L(1,X) = 1 - X
!    L(N,X) = (2*N-1-X)/N * L(N-1,X) - (N-1)/N * L(N-2,X)
!
!  Orthogonality:
!
!    Integral ( 0 <= X < oo ) exp ( - X ) * L(N,X) * L(M,X) dX = delta ( M, N )
!
!  Relations:
!
!    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
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
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  v(1:m,1) = 1.0D+00 - x(1:m)

  do j = 2, n

    v(1:m,j) = ( ( real ( 2 * j - 1, kind = 8 ) - x(1:m) ) * v(1:m,j-1)   &
                 - real (     j - 1, kind = 8 )            * v(1:m,j-2) ) &
                 / real (     j,     kind = 8 )

  end do

  return
end
subroutine l_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! L_POLYNOMIAL_COEFFICIENTS: coefficients of the Laguerre polynomial L(n,x).
!
!  First terms:
!
!    0: 1
!    1: 1  -1
!    2: 1  -2  1/2
!    3: 1  -3  3/2  1/6
!    4: 1  -4  4   -2/3  1/24
!    5: 1  -5  5   -5/3  5/24  -1/120
!
!  Recursion:
!
!    L(0) = ( 1,  0, 0, ..., 0 )
!    L(1) = ( 1, -1, 0, ..., 0 )
!    L(N) = (2*N-1-X) * L(N-1) - (N-1) * L(N-2) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the
!    Laguerre polynomials of degree 0 through N.  Each polynomial
!    is stored as a row.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0:n,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = -1.0D+00

  do i = 2, n

    c(i,1:n) = ( &
        real ( 2 * i - 1, kind = 8 ) * c(i-1,1:n)     &
      + real (   - i + 1, kind = 8 ) * c(i-2,1:n)     &
      -                                c(i-1,0:n-1) ) &
      / real (     i,     kind = 8 )

  end do

  return
end
subroutine l_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! L_POLYNOMIAL_VALUES: some values of the Laguerre polynomial L(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LaguerreL[n,x]
!
!  Differential equation:
!
!    X * Y'' + (1-X) * Y' + N * Y = 0
!
!  First terms:
!
!      1
!     -X    +  1
!   (  X^2 -  4 X     +  2 ) / 2
!   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
!   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
!   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
!   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
!   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3 + 52920 X^2 - 35280 X
!     + 5040 ) / 5040
!
!  Recursion:
!
!    L(0,X) = 1,
!    L(1,X) = 1-X,
!    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
!
!  Orthogonality:
!
!    Integral ( 0 <= X < oo ) exp ( - X ) * L(N,X) * L(M,X) dX
!    = 0 if N /= M
!    = 1 if N == M
!
!  Special values:
!
!    L(N,0) = 1.
!
!  Relations:
!
!    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.0000000000000000D+00, &
    -0.5000000000000000D+00, &
    -0.6666666666666667D+00, &
    -0.6250000000000000D+00, &
    -0.4666666666666667D+00, &
    -0.2569444444444444D+00, &
    -0.4047619047619048D-01, &
     0.1539930555555556D+00, &
     0.3097442680776014D+00, &
     0.4189459325396825D+00, &
     0.4801341790925124D+00, &
     0.4962122235082305D+00, &
    -0.4455729166666667D+00, &
     0.8500000000000000D+00, &
    -0.3166666666666667D+01, &
     0.3433333333333333D+02 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    0.5D+00, &
    3.0D+00, &
    5.0D+00, &
    1.0D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine l_polynomial_zeros ( n, x )

!*****************************************************************************80
!
!! L_POLYNOMIAL_ZEROS: zeros of the Laguerre polynomial L(n,x).
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
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X(N), the zeros.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = 1.0D+00
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    x(i) = real ( 2 * i - 1, kind = 8 )
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  return
end
subroutine l_power_product ( p, e, table )

!*****************************************************************************80
!
!! L_POWER_PRODUCT: power product table for L(n,x).
!
!  Discussion:
!
!    Let L(n,x) represent the Laguerre polynomial of degree n.
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of X^E with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( 0 <= X < +oo ) x^e * L(i,x) * L(j,x) * exp (-x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial
!    factors.  0 <= P.
!
!    Input, integer ( kind = 4 ) E, the exponent of X.
!    0 <= E.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.
!    TABLE(I,J) represents the weighted integral of x^E * L(i,x) * L(j,x).
!
  implicit none

  integer ( kind = 4 ) p

  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), allocatable :: l_table(:)
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = p + 1 + ( e + 1 ) / 2

  allocate ( x_table(1:order) )
  allocate ( w_table(1:order) )

  call l_quadrature_rule ( order, x_table, w_table )

  allocate ( l_table(0:p) )

  do k = 1, order

    x = x_table(k)
    call l_polynomial ( 1, p, x, l_table )

    if ( e == 0 ) then

      do j = 0, p
        do i = 0, p
          table(i,j) = table(i,j) + w_table(k) * l_table(i) * l_table(j)
        end do
      end do

    else

      do j = 0, p
        do i = 0, p
          table(i,j) = table(i,j) &
            + w_table(k) * ( x ** e ) * l_table(i) * l_table(j)
        end do
      end do

    end if

  end do

  deallocate ( l_table )
  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine l_quadrature_rule ( n, x, w )

!*****************************************************************************80
!
!! L_QUADRATURE_RULE: Gauss-Laguerre quadrature based on L(n,x).
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
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = 1.0D+00
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    x(i) = real ( 2 * i - 1, kind = 8 )
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine lf_function ( m, n, alpha, x, cx )

!*****************************************************************************80
!
!! LF_FUNCTION evaluates the Laguerre function Lf(n,alpha,x).
!
!  Recursion:
!
!    Lf(0,ALPHA,X) = 1
!    Lf(1,ALPHA,X) = 1+ALPHA-X
!
!    Lf(N,ALPHA,X) = (2*N-1+ALPHA-X)/N * Lf(N-1,ALPHA,X)
!                      - (N-1+ALPHA)/N * Lf(N-2,ALPHA,X)
!
!  Restrictions:
!
!    -1 < ALPHA
!
!  Special values:
!
!    Lf(N,0,X) = L(N,X).
!    Lf(N,ALPHA,X) = LM(N,ALPHA,X) for ALPHA integral.
!
!  Norm:
!
!    Integral ( 0 <= X < +oo ) exp ( - X ) * Lf(N,ALPHA,X)^2 dX
!    = Gamma ( N + ALPHA + 1 ) / N!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order function to compute.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.  -1 < ALPHA is required.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(1:M,0:N), the functions of
!    degrees 0 through N evaluated at the points X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) cx(1:m,0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(1:m)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LF_FUNCTION - Fatal error!'
    write ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
    write ( *, '(a)' ) '  but ALPHA must be greater than -1.'
    stop
  end if

  if ( n < 0 ) then
    return
  end if

  cx(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1:m,1) = 1.0D+00 + alpha - x(1:m)

  do i = 2, n
    cx(1:m,i) = ( &
      ( real ( 2 * i - 1, kind = 8 ) + alpha - x(1:m) ) * cx(1:m,i-1)   &
    + ( real (   - i + 1, kind = 8 ) - alpha          ) * cx(1:m,i-2) ) &
      / real (     i,     kind = 8 )
  end do

  return
end
subroutine lf_function_values ( n_data, n, a, x, fx )

!*****************************************************************************80
!
!! LF_FUNCTION_VALUES returns values of the Laguerre function Lf(n,alpha,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LaguerreL[n,a,x]
!
!    The functions satisfy the following differential equation:
!
!      X * Y'' + (ALPHA+1-X) * Y' + N * Y = 0
!
!    Function values can be generated by the recursion:
!
!      Lf(0,ALPHA,X) = 1
!      Lf(1,ALPHA,X) = 1+ALPHA-X
!
!      Lf(N,ALPHA,X) = ( (2*N-1+ALPHA-X) * Lf(N-1,ALPHA,X)
!                          - (N-1+ALPHA) * Lf(N-2,ALPHA,X) ) / N
!
!    The parameter ALPHA is required to be greater than -1.
!
!    For ALPHA = 0, the generalized Laguerre function Lf(N,ALPHA,X)
!    is equal to the Laguerre polynomial L(N,X).
!
!    For ALPHA integral, the generalized Laguerre function
!    Lf(N,ALPHA,X) equals the associated Laguerre polynomial Lm(N,ALPHA,X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) A, the parameter.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.00D+00, &
    0.25D+00, &
    0.50D+00, &
    0.75D+00, &
    1.50D+00, &
    2.50D+00, &
    5.00D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.3726399739583333D-01, &
     0.3494791666666667D+00, &
     0.8710042317708333D+00, &
     0.1672395833333333D+01, &
     0.6657625325520833D+01, &
     0.2395726725260417D+02, &
     0.2031344319661458D+03, &
     0.1284193996800000D+02, &
     0.5359924801587302D+01, &
     0.9204589064126984D+00, &
    -0.1341585114857143D+01, &
    -0.2119726307555556D+01, &
    -0.1959193658349206D+01, &
     0.1000000000000000D+01, &
     0.5450000000000000D+01, &
     0.1720125000000000D+02, &
     0.4110393750000000D+02, &
     0.8239745859375000D+02, &
     0.1460179186171875D+03, &
     0.2359204608298828D+03 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     8, &
     8, &
     8, &
     8, &
     8, &
     8, &
     0, &
     1, &
     2, &
     3, &
     4, &
     5, &
     6 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.00D+00, &
    0.20D+00, &
    0.40D+00, &
    0.60D+00, &
    0.80D+00, &
    1.00D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine lf_function_zeros ( n, alpha, x )

!*****************************************************************************80
!
!! LF_FUNCTION_ZEROS returns the zeros of Lf(n,alpha,x).
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
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    ALPHA must be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the zeros.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) i_r8
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_gamma ( alpha + 1.0D+00 )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    bj(i) = i_r8 * ( i_r8 + alpha )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    x(i) = 2.0D+00 * i_r8 - 1.0D+00 + alpha
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  return
end
subroutine lf_integral ( n, alpha, exact )

!*****************************************************************************80
!
!! LF_INTEGRAL evaluates a monomial integral associated with Lf(n,alpha,x).
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) x^n * x^alpha * exp ( -x ) dx
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the exponent.
!    0 <= N.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) exact
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_gamma

  arg = alpha + real ( n + 1, kind = 8 )

  exact = r8_gamma ( arg )

  return
end
subroutine lf_quadrature_rule ( n, alpha, x, w )

!*****************************************************************************80
!
!! LF_QUADRATURE_RULE: Gauss-Laguerre quadrature rule for Lf(n,alpha,x);
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) exp ( - x ) * x^alpha * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    ALPHA must be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) i_r8
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_gamma ( alpha + 1.0D+00 )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    bj(i) = i_r8 * ( i_r8 + alpha )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    x(i) = 2.0D+00 * i_r8 - 1.0D+00 + alpha
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine lm_integral ( n, m, exact ) !???

!*****************************************************************************80
!
!! LM_INTEGRAL evaluates a monomial integral associated with Lm(n,m,x).
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) x^n * x^m * exp ( -x ) dx
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the exponent.
!    0 <= N.
!
!    Input, integer ( kind = 4 ) M, the parameter.
!    0 <= M.
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial

  exact = r8_factorial ( n + m )

  return
end
subroutine lm_polynomial_coefficients ( n, m, c )

!*****************************************************************************80
!
!! LM_POLYNOMIAL_COEFFICIENTS: coefficients of Laguerre polynomial Lm(n,m,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, integer ( kind = 4 ) M, the parameter.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the
!    Laguerre polynomials of degree 0 through N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,0) = real ( m + 1, kind = 8 )
  c(1,1) = -1.0D+00

  do i = 2, n

    c(i,0:i) = ( real (   m + 2 * i - 1, kind = 8 ) * c(i-1,0:i)   &
               + real ( - m     - i + 1, kind = 8 ) * c(i-2,0:i) ) &
               / real (           i,     kind = 8 )

    c(i,1:i) = c(i,1:i) - c(i-1,0:i-1) / real ( i, kind = 8 )

  end do

  return
end
subroutine lm_polynomial_values ( n_data, n, m, x, fx ) !table values

!*****************************************************************************80
!
!! LM_POLYNOMIAL_VALUES returns values of Laguerre polynomials Lm(n,m,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LaguerreL[n,m,x]
!
!    The associated Laguerre polynomials may be generalized so that the
!    parameter M is allowed to take on arbitrary noninteger values.
!    The resulting function is known as the generalized Laguerre function.
!
!    The polynomials satisfy the differential equation:
!
!      X * Y'' + (M+1-X) * Y' + (N-M) * Y = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, integer ( kind = 4 ) M, the parameter.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1500000000000000D+01, &
    0.1625000000000000D+01, &
    0.1479166666666667D+01, &
    0.1148437500000000D+01, &
    0.4586666666666667D+00, &
    0.2878666666666667D+01, &
    0.8098666666666667D+01, &
    0.1711866666666667D+02, &
    0.1045328776041667D+02, &
    0.1329019368489583D+02, &
    0.5622453647189670D+02, &
    0.7484729341779436D+02, &
    0.3238912982762806D+03, &
    0.4426100000097533D+03, &
    0.1936876572288250D+04 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 0, 0, &
    0, 1, 1, 1, &
    1, 0, 1, 2, &
    3, 2, 2, 3, &
    3, 4, 4, 5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    1,  2,  3,  4, &
    5,  1,  2,  3, &
    4,  3,  3,  3, &
    3,  4,  5,  6, &
    7,  8,  9, 10 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine lm_polynomial_zeros ( n, m, x )

!*****************************************************************************80
!
!! LM_POLYNOMIAL_ZEROS returns the zeros for Lm(n,m,x).
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
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, integer ( kind = 4 ) M, the parameter.
!    0 <= M.
!
!    Output, real ( kind = 8 ) X(N), the zeros.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_factorial ( m )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i * ( i + m ), kind = 8 )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  do i = 1, n
    x(i) = real ( 2 * i - 1 + m, kind = 8 )
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  return
end
subroutine lm_quadrature_rule ( n, m, x, w )

!*****************************************************************************80
!
!! LM_QUADRATURE_RULE: Gauss-Laguerre quadrature rule for Lm(n,m,x);
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) exp ( - x ) * x^m * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
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
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, integer ( kind = 4 ) M, the parameter.
!    0 <= M.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_factorial ( m )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i * ( i + m ), kind = 8 )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  do i = 1, n
    x(i) = real ( 2 * i - 1 + m, kind = 8 )
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( n ) = product ( 1 <= i <= n ) i
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = 8 )
  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
subroutine lm_polynomial ( mm, n, m, x, cx )

!*****************************************************************************80
!
!! LM_POLYNOMIAL evaluates Laguerre polynomials Lm(n,m,x).
!
!  First terms:
!
!    M = 0
!
!    Lm(0,0,X) =   1
!    Lm(1,0,X) =  -X   +  1
!    Lm(2,0,X) =   X^2 -  4 X   +  2
!    Lm(3,0,X) =  -X^3 +  9 X^2 -  18 X   +    6
!    Lm(4,0,X) =   X^4 - 16 X^3 +  72 X^2 -   96 X +     24
!    Lm(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x   +  120
!    Lm(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720
!
!    M = 1
!
!    Lm(0,1,X) =    0
!    Lm(1,1,X) =   -1,
!    Lm(2,1,X) =    2 X - 4,
!    Lm(3,1,X) =   -3 X^2 + 18 X - 18,
!    Lm(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
!
!    M = 2
!
!    Lm(0,2,X) =    0
!    Lm(1,2,X) =    0,
!    Lm(2,2,X) =    2,
!    Lm(3,2,X) =   -6 X + 18,
!    Lm(4,2,X) =   12 X^2 - 96 X + 144
!
!    M = 3
!
!    Lm(0,3,X) =    0
!    Lm(1,3,X) =    0,
!    Lm(2,3,X) =    0,
!    Lm(3,3,X) =   -6,
!    Lm(4,3,X) =   24 X - 96
!
!    M = 4
!
!    Lm(0,4,X) =    0
!    Lm(1,4,X) =    0
!    Lm(2,4,X) =    0
!    Lm(3,4,X) =    0
!    Lm(4,4,X) =   24
!
!  Recursion:
!
!    Lm(0,M,X)   = 1
!    Lm(1,M,X)   = (M+1-X)
!
!    if 2 <= N:
!
!      Lm(N,M,X)   = ( (M+2*N-1-X) * Lm(N-1,M,X)
!                   +   (1-M-N)    * Lm(N-2,M,X) ) / N
!
!  Special values:
!
!    For M = 0, the associated Laguerre polynomials Lm(N,M,X) are equal
!    to the Laguerre polynomials L(N,X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, integer ( kind = 4 ) M, the parameter.  M must be nonnegative.
!
!    Input, real ( kind = 8 ) X(MM), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the associated Laguerre polynomials
!    of degrees 0 through N evaluated at the evaluation points.
!
  implicit none

  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(mm,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(mm)

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LM_POLYNOMIAL - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M = ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop
  end if

  if ( n < 0 ) then
    return
  end if

  cx(1:mm,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1:mm,1) = real ( m + 1, kind = 8 ) - x(1:mm)

  do i = 2, n
    cx(1:mm,i) = &
      ( ( real (   m + 2 * i - 1, kind = 8 ) - x(1:mm) ) * cx(1:mm,i-1)   &
        + real ( - m     - i + 1, kind = 8 )             * cx(1:mm,i-2) ) &
        / real (           i,     kind = 8 )
  end do

  return
end
!****************************** Subroutines de Pablo *********************************

  SUBROUTINE gauchv(x,wn,n)
!-----------------------------------------------------------------
!
! Computes the abscissas and weights for Gauss-Chebyshev n-point
! quadrature. The weigh function for Gauss-Chebyshev quadrature
! is W(x) = 1/sqrt(1-x^2), and the limits of integration are -1
! and 1.
!
! Parameters
!    x  : at the output contains the abscissas
!    w  : at the output contains the weights
!    n  : number of abscissas and weights to compute
!
      USE constants
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: n
      DOUBLE PRECISION, INTENT(OUT) :: x(n),wn(n)
      INTEGER                       :: i

      DO i = 1,n
         x(i) = cos(PI*(i-.5d0)/n)
         wn(i) = PI/n
      END DO

      RETURN
      END SUBROUTINE gauchv
  SUBROUTINE genext(ind,ext)
!*****************************************************************
!
! Converts the index ind into a character.
!
! Parameters
!    ind: index
!    ext: character
!
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: ind
      INTEGER                       :: ic,id,iu
      CHARACTER(len=3), INTENT(OUT) :: ext

      ic = 48+int(ind/100)
      id = 48+int(ind/10)-int(ind/100)*10
      iu = 48+int(ind)-int(ind/10)*10
      ext = char(ic) // char(id) // char(iu)

      RETURN
      END SUBROUTINE genext


