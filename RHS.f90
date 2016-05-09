Module rhs
    !to do: complex version
    implicit none
    contains

    subroutine df_paper(coefs,p,m,j,t,dxdy)
        use constants
        use resolution

        implicit none

        integer*4, intent(in) :: p,m,j
        real*8 :: t
        real*8, dimension(q,l,i,5), intent(in) :: coefs
        real*8, dimension(q,l,i,5), intent(out) :: dxdy
        real*8, dimension(q,l,i,1) :: chi
        real*8, dimension(q,l,i,1) :: nonlincou

        dxdy(p,m,j,1)=-k*coefs(p,m,j,1)+(-a+d)*coefs(p,m,j,2)-2*k*coefs(p,m,j,3)
        dxdy(p,m,j,2)=-k*coefs(p,m,j,2)+(a-d)*coefs(p,m,j,1)-2*k*coefs(p,m,j,4)
        dxdy(p,m,j,3)=-nonlincou(p,m,j,1)+coefs(p,m,j,3)-d*coefs(p,m,j,4)
        dxdy(p,m,j,4)=-nonlincou(p,m,j,1)+coefs(p,m,j,4)+d*coefs(p,m,j,3)
        dxdy(p,m,j,5)=-g*(-1/2*(nonlincou(p,m,j,1)+nonlincou(p,m,j,1))+coefs(p,m,j,5)-chi(p,m,j,1))
    end subroutine

    subroutine df_test(c0,p,m,j,t,c1)
        use constants
        use resolution

        implicit none

        integer*4, intent(in) :: p,m,j
        real*8 :: t
        real*8, dimension(q,l,i,5), intent(in) :: c0
        real*8, dimension(q,l,i,5), intent(out) :: c1
        real*8, dimension(q,l,i,1) :: chi
        real*8, dimension(q,l,i,1) :: nonlincou

        c1(p,m,j,1)=-k*c0(p,m,j,1)+(-a+d)*c0(p,m,j,2)-2*k*c0(p,m,j,3)
        c1(p,m,j,2)=-k*c0(p,m,j,2)+(a-d)*c0(p,m,j,1)-2*k*c0(p,m,j,4)
        c1(p,m,j,3)=c0(p,m,j,3)-d*c0(p,m,j,4)
        c1(p,m,j,4)=c0(p,m,j,4)+d*c0(p,m,j,3)
        c1(p,m,j,5)=-g*(-1/2+c0(p,m,j,5))
    end subroutine

!subroutine df_fresnel(coefs,p,m,j,t,dxdy)
!        use constants
!        use resolution
!
!        implicit none
!        integer*4, intent(in) :: p,m,j
!        real*8 :: t
!        real*8, dimension(q,l,i,9), intent(in) :: coefs
!        real*8, dimension(q,l,i,1) :: dxdy1
!        real*8, dimension(q,l,i,1) :: chi
!        real*8, dimension(q,l,i,1) :: nonlincou
!
!
!        dxdy1(p,m,j,1)=-k*coefs(p,m,j,1)+(-a+d)*coefs(p,m,j,2)-2*k*coefs(p,m,j,5)
!        dxdy2(p,m,j,2)=-k*coefs(p,m,j,2)+(a-d)*coefs(p,m,j,1)-2*k*coefs(p,m,j,6)
!        dxdy3(p,m,j,3)=-k*coefs(p,m,j,3)+(-a+d)*coefs(p,m,j,4)-2*k*coefs(p,m,j,7)
!        dxdy4(p,m,j,4)=-k*coefs(p,m,j,4)+(a-d)*coefs(p,m,j,3)-2*k*coefs(p,m,j,8)
!        dxdy5(p,m,j,5)=-nonlincou(p,m,j)+coefs(p,m,j,5)-d*coefs(p,m,j,6)
!        dxdy6(p,m,j,6)=-nonlincou(p,m,j)+coefs(p,m,j,6)+d*coefs(p,m,j,5)
!        dxdy7(p,m,j,7)=-nonlincou(p,m,j)+coefs(p,m,j,7)-d*coefs(p,m,j,8)
!        dxdy8(p,m,j,8)=-nonlincou(p,m,j)+coefs(p,m,j,8)+d*coefs(p,m,j,7)
!        dxdy9(p,m,j,9)=-g*(-1/2*(nonlincou(p,m,j)+nonlincou(p,m,j))+coefs(p,m,j,9)-chi(p,m,j,1))
!    end subroutine


!subroutine df1(coef,p,m,j,t,dxdyout,n)
!        use constants
!
!        integer*4, intent(in) :: p,m,j,n
!        real*8 :: t
!        real*8, dimension(q,l,i,9), intent(in) :: coefs,dxdy
!        real*8, dimension(q,l,i,1) :: dxdyout
!        real*8, dimension(q,l,i,1) :: chi
!        real*8, dimension(q,l,i,1) :: nonlincou
!
!
!        dxdy(p,m,j,1)=-k*coefs(p,m,j,1)+(-a+d)*coefs(p,m,j,2)-2*k*coefs(p,m,j,5)
!        dxdy(p,m,j,2)=-k*coefs(p,m,j,2)+(a-d)*coefs(p,m,j,1)-2*k*coefs(p,m,j,6)
!        dxdy(p,m,j,3)=-k*coefs(p,m,j,3)+(-a+d)*coefs(p,m,j,4)-2*k*coefs(p,m,j,7)
!        dxdy(p,m,j,4)=-k*coefs(p,m,j,4)+(a-d)*coefs(p,m,j,3)-2*k*coefs(p,m,j,8)
!        dxdy(p,m,j,5)=-nonlincou(p,m,j)+coefs(p,m,j,5)-d*coefs(p,m,j,6)
!        dxdy(p,m,j,6)=-nonlincou(p,m,j)+coefs(p,m,j,6)+d*coefs(p,m,j,5)
!        dxdy(p,m,j,7)=-nonlincou(p,m,j)+coefs(p,m,j,7)-d*coefs(p,m,j,8)
!        dxdy(p,m,j,8)=-nonlincou(p,m,j)+coefs(p,m,j,8)+d*coefs(p,m,j,7)
!        dxdy(p,m,j,9)=-g*(-1/2*(nonlincou(p,m,j)+nonlincou(p,m,j))+coefs(p,m,j,9)-chi(p,m,j))
!
!        if (n.eq.1) then
!            dxdyout(p,m,j)= dxdy(p,m,j,n)
!        end if
!end subroutine

End Module
