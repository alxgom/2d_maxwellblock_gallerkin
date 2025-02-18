module rk
    implicit none
    contains


    subroutine rk4 ( t0, u0, dt, f, u )

    !*****************************************************************************80
    !
    !! RK4 takes one Runge-Kutta step for a scalar ODE.
    !
    !  Discussion:
    !
    !    It is assumed that an initial value problem, of the form
    !
    !      du/dt = f ( t, u )
    !      u(t0) = u0
    !
    !    is being solved.
    !
    !    If the user can supply current values of t, u, a stepsize dt, and a
    !    function to evaluate the derivative, this function can compute the
    !    fourth-order Runge Kutta estimate to the solution at time t+dt.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    31 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T0, the current time.
    !
    !    Input, real ( kind = 8 ) U0, the solution estimate at the current time.
    !
    !    Input, real ( kind = 8 ) DT, the time step.
    !
    !    Input, external F, a subroutine of the form
    !      subroutine f ( t, u, uprime )
    !    which evaluates the derivative uprime given the time T and
    !    solution vector U.
    !
    !    Output, real ( kind = 8 ) U, the fourth-order Runge-Kutta solution
    !    estimate at time T0+DT.
    !
      implicit none

      real ( kind = 8 ) dt
      external f
      real ( kind = 8 ) f0
      real ( kind = 8 ) f1
      real ( kind = 8 ) f2
      real ( kind = 8 ) f3
      real ( kind = 8 ) t0
      real ( kind = 8 ) t1
      real ( kind = 8 ) t2
      real ( kind = 8 ) t3
      real ( kind = 8 ) u
      real ( kind = 8 ) u0
      real ( kind = 8 ) u1
      real ( kind = 8 ) u2
      real ( kind = 8 ) u3
    !
    !  Get four sample values of the derivative.
    !
      call f ( t0, u0, f0 )

      t1 = t0 + dt / 2.0D+00
      u1 = u0 + dt * f0 / 2.0D+00
      call f ( t1, u1, f1 )

      t2 = t0 + dt / 2.0D+00
      u2 = u0 + dt * f1 / 2.0D+00
      call f ( t2, u2, f2 )

      t3 = t0 + dt
      u3 = u0 + dt * f2
      call f ( t3, u3, f3 )
    !
    !  Combine them to estimate the solution U at time T1.
    !
      u = u0 + dt * ( f0 + 2.0D+00 * f1 + 2.0D+00 * f2 + f3 ) / 6.0D+00

      return
    end subroutine

    subroutine rk4vec ( t0, m, u0, dt, f, u )

    !*****************************************************************************80
    !
    !! RK4VEC takes one Runge-Kutta step for a vector ODE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    08 October 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T0, the current time.
    !
    !    Input, integer ( kind = 4 ) M, the dimension of the system.
    !
    !    Input, real ( kind = 8 ) U0(M), the solution estimate at the current time.
    !
    !    Input, real ( kind = 8 ) DT, the time step.
    !
    !    Input, external F, a subroutine of the form
    !      subroutine f ( t, m, u, uprime )
    !    which evaluates the derivative UPRIME(1:M) given the time T and
    !    solution vector U(1:M).
    !
    !    Output, real ( kind = 8 ) U(M), the fourth-order Runge-Kutta solution
    !    estimate at time T0+DT.
    !
      implicit none

      integer ( kind = 4 ) m

      real ( kind = 8 ) dt
      external f
      real ( kind = 8 ) f0(m)
      real ( kind = 8 ) f1(m)
      real ( kind = 8 ) f2(m)
      real ( kind = 8 ) f3(m)
      real ( kind = 8 ) t0
      real ( kind = 8 ) t1
      real ( kind = 8 ) t2
      real ( kind = 8 ) t3
      real ( kind = 8 ) u(m)
      real ( kind = 8 ) u0(m)
      real ( kind = 8 ) u1(m)
      real ( kind = 8 ) u2(m)
      real ( kind = 8 ) u3(m)
    !
    !  Get four sample values of the derivative.
    !
      call f ( t0, m, u0, f0 )

      t1 = t0 + dt / 2.0D+00
      u1(1:m) = u0(1:m) + dt * f0(1:m) / 2.0D+00
      call f ( t1, m, u1, f1 )

      t2 = t0 + dt / 2.0D+00
      u2(1:m) = u0(1:m) + dt * f1(1:m) / 2.0D+00
      call f ( t2, m, u2, f2 )

      t3 = t0 + dt
      u3(1:m) = u0(1:m) + dt * f2(1:m)
      call f ( t1, m, u1, f3 )
    !
    !  Combine them to estimate the solution U at time T1.
    !
      u(1:m) = u0(1:m) + dt * ( f0(1:m) + 2.0D+00 * f1(1:m) + 2.0D+00 * f2(1:m) &
        + f3(1:m) ) / 6.0D+00

      return
    end subroutine

end module
