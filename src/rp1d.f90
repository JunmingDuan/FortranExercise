!**
! @file rp1d.f90
! @brief A simple Fortran exercise for 1D Euler
! @author Duan Junming, duanjm@pku.edu.cn
! @version 1.0
! @date 2019-11-03
!/

Subroutine F0(x, pri)
  Use GlobalData
  Use Euler1d
  Implicit None
  Double Precision :: x, pri(NCOMP)
  pri(RHO) = 1 + 0.5*Sin(2*pi*x)
  pri(VX1) = 0.5
  pri(PRS) = 1.0
  ! If(x < 0.5) then
  !   pri(RHO) = 1
  !   pri(VX1) = 0
  !   pri(PRS) = 1
  ! Else
  !   pri(RHO) = 0.125
  !   pri(VX1) = 0
  !   pri(PRS) = 0.1
  ! End If
End Subroutine F0

Subroutine GetInit(left,right,res)
  Use GlobalData
  Implicit None
  Double precision, intent(in)::left, right
  Double precision, intent(out)::res(NCOMP)
  Double precision::gp(3),gw(3),tmp(NCOMP)
  Real:: h,center
  Integer :: i
  gp=(/-0.7745966692,          0.0,  0.7745966692/)
  gw=(/ 0.5555555556, 0.8888888889,  0.5555555556/)
  h=right-left
  center=(right+left)/2.0
  Do i=1,3
    Call F0(center+gp(i)*h,tmp)
    res=res+gw(i)*tmp  ! Variable initial
  End Do
End subroutine GetInit

Subroutine Init
  Use GlobalData
  Implicit None
  Integer :: i
  Do i = nx_beg, nx_end
    ! Call F0(0.5*(grid(i-1)+grid(i)), prim(i,:))
    Call GetInit(grid(i-1),grid(i),prim(i,:))
    Call Pri2Con(prim(i,:), adiabatic_index, cons(i,:))
  End Do
End Subroutine Init

Subroutine Boundary
  Use GlobalData
  Implicit None
  Integer :: i

  Do i = nx_glbbeg, nx_beg-1
    If(xbeg_bc_type == 0) Then
      prim(i,:) = prim(i+nx,:)
      cons(i,:) = cons(i+nx,:)
    Else If(xbeg_bc_type == 1) Then
      prim(i,:) = prim(nx_beg,:)
      cons(i,:) = cons(nx_beg,:)
    End If
  End Do

  Do i = nx_end+1, nx_glbend
    If(xend_bc_type == 0) Then
      prim(i,:) = prim(i-nx,:)
      cons(i,:) = cons(i-nx,:)
    Else If(xend_bc_type == 1) Then
      prim(i,:) = prim(nx_end,:)
      cons(i,:) = cons(nx_end,:)
    End If
  End Do
End Subroutine Boundary

Subroutine Reconstruct
  Use GlobalData
  Use Reconstruction
  Implicit None
  Integer :: i

  Do i = nx_beg-1, nx_end
    Call Weno5(cons(i-nx_ghost+1:i+nx_ghost-1,:), rec_cons_m(i,:), NCOMP)
    Call Weno5(cons(i+nx_ghost:i-nx_ghost+2:-1,:), rec_cons_p(i,:), NCOMP)
  End Do
  !Pause
End Subroutine Reconstruct

Subroutine CompFlux(stage)
  Use GlobalData
  Implicit None
  Integer :: i, stage
  Double Precision :: alpha

  If(stage == 1) g_dt = 1E10
  Do i = nx_beg-1, nx_end
    Call LLf(rec_cons_m(i,:), rec_cons_p(i,:), adiabatic_index, alpha, flux(i,:))
    !Call LLf(cons(i,:), cons(i+1,:), adiabatic_index, alpha, flux(i,:))
    If(stage == 1) Then
      g_dt = Min(g_dt, cfl * (grid(i) - grid(i-1)) / alpha)
    End If
  End Do
  If(stage == 1 .And. g_t + g_dt > t_end) Then
    g_dt = t_end - g_t
  End If
End Subroutine CompFlux

!Subroutine ComputError
  !Use GlobalData
  !Implicit None
  !Integer :: i
  !Double Precision :: error(3)

  !error = 0;
  !Open(1, File="error.dat")
  !Do i = nx_beg, nx_end
    !0.5*(grid(i-1)+grid(i)), prim(i,:)
  !End Do
  !Write(1,*) nx, error
  !Close(1)
!End Subroutine PrintTecplot

Subroutine PrintTecplot
  Use GlobalData
  Implicit None
  Integer :: i

  Open(1, File="sol.dat")
  Write(1,*) "Title="//'"mesh and solution"'
  Write(1,*) "variables="//'"x",'//'"rho",'//'"vx",'//'"prs"'
  Write(1,*) 'zone t="', g_t, '", I=', nx
  Do i = nx_beg, nx_end
    Write(1,*) 0.5*(grid(i-1)+grid(i)), prim(i,:)
  End Do
  Close(1)
End Subroutine PrintTecplot

Program rp1d
  Use GlobalData
  Use TimeInfo
  Implicit None
  Integer :: i

  Print *, "It is Riemann 1d."
  Call PrintTime
  Call ReadParameter
  Call InitMesh
  Call Init

  ! Do While(g_t < t_end - 1E-14)
  !   Call RK3
  !   !Call ForwardEuler
  !   g_t = g_t + g_dt
  !   g_step = g_step + 1
  !   Print *, "step = ", g_step, ", t = ", g_t, ", dt = ", g_dt
  ! End Do

  !Call ComputError
  Call PrintTecplot
  Call DestroyMesh
End Program rp1d

