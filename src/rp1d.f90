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

Subroutine ExactFun(x, t, pri)
  Use GlobalData
  Use Euler1d
  Implicit None
  Double Precision :: x, t, pri(NCOMP)
  ! Write(*,*) pi ! no error???????

  pri(RHO) = 1.0 + 0.5*Sin(2*pi*(x-0.5*t))
  pri(VX1) = 0.5
  pri(PRS) = 1.0
End Subroutine

Subroutine GuassIntegration(left, right, t, res)
  Use GlobalData
  Implicit None
  Double precision, intent(in)  :: left, right, t
  Double precision, intent(out) :: res(NCOMP)
  Double precision :: gp(3), gw(3), tmp_pri(NCOMP), tmp_con(NCOMP)
  Double precision :: half_h, center
  Integer :: i
  ! gp=(/-0.7745966692,          0.0,  0.7745966692/)
  ! gw=(/ 0.5555555556, 0.8888888889,  0.5555555556/)
  gp = [-Sqrt(3.0/5.0), 0.0, Sqrt(3.0/5.0) ]
  gw = [5.d0/9.0, 8.d0/9.0, 5.d0/9.0 ]
  half_h = (right - left)/2.0
  center = (right + left)/2.0
  res = 0
  Do i=1, 3
    Call ExactFun(center + gp(i)*half_h, t, tmp_pri)
    Call Pri2Con(tmp_pri, adiabatic_index, tmp_con)
    res = res + gw(i)*tmp_con
  End Do
  res = res/2.0
End subroutine GuassIntegration

Subroutine Init
  Use GlobalData
  Implicit None
  Integer :: i
  Double precision::t
  t = 0
  Do i = nx_beg, nx_end
    ! Call F0(0.5*(grid(i-1)+grid(i)), prim(i,:))
    Call GuassIntegration(grid(i-1), grid(i), t, cons(i,:))
    Call Con2Pri(cons(i,:), adiabatic_index, prim(i,:))
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

Subroutine EigMatrix(cell_prim, cell_cons, matrixl, matrixr)
  Use GlobalData
  Implicit None
  Double Precision :: h, c, cvalue
  Double Precision :: cell_prim(NCOMP), cell_cons(NCOMP), matrixl(NCOMP,NCOMP), matrixr(NCOMP, NCOMP), ans(NCOMP, NCOMP)

  h = (cell_cons(PRS)+cell_prim(PRS))/cell_prim(RHO)
  c = sqrt(adiabatic_index*cell_prim(PRS)/cell_prim(RHO))
  cvalue = -cell_prim(VX1)**2 + 2*h; 

  matrixr(1,:) = 1.0
  matrixr(2,:) = [cell_prim(VX1)-c,   cell_prim(VX1),        cell_prim(VX1)+c]
  matrixr(3,:) = [h-cell_prim(VX1)*c, 0.5*cell_prim(VX1)**2, h+cell_prim(VX1)*c]

  matrixl(1,:) = [(-cell_prim(VX1)**3 + c*cell_prim(VX1)**2 + 2*h*cell_prim(VX1))/(2*c*cvalue),&
                 -(-cell_prim(VX1)**2 + 2*c*cell_prim(VX1) + 2*h)/(2*c*cvalue), 1.0/cvalue]
  matrixl(2,:) = [2*(-cell_prim(VX1)**2 + h)/cvalue, 2*cell_prim(VX1)/cvalue, -2.0/cvalue]
  matrixl(3,:) = [(cell_prim(VX1)**3 + c*cell_prim(VX1)**2 - 2*h*cell_prim(VX1))/(2*c*cvalue),&
                  -(cell_prim(VX1)**2 + 2*c*cell_prim(VX1) - 2*h)/(2*c*cvalue), 1.0/cvalue]

!ans = matmul(matrixl, matrixr)
!write(*,'(3(3f6.2,/))') ans 

End Subroutine EigMatrix

Subroutine Reconstruct
  Use GlobalData
  Use Reconstruction
  Implicit None
  Integer :: i
  Double Precision :: chara(2*nx_ghost-1, NCOMP), rec_chara(NCOMP), matl(NCOMP,NCOMP), matr(NCOMP,NCOMP)

  Do i = nx_beg-1, nx_end
    ! Call EigMatrix(prim(i,:), cons(i,:), matl, matr)
    ! chara = matmul(cons(i-nx_ghost+1:i+nx_ghost-1,:), transpose(matl))
    ! Call Weno5(chara, rec_chara, NCOMP)
    ! rec_cons_m(i,:) = matmul(matr, rec_chara)
    
    ! Call EigMatrix(prim(i+1,:), cons(i+1,:), matl, matr)
    ! chara = matmul(cons(i+nx_ghost:i-nx_ghost+2:-1,:), transpose(matl))
    ! Call Weno5(chara, rec_chara, NCOMP)
    ! rec_cons_p(i,:) = matmul(matr, rec_chara)

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
    ! Call LLf(cons(i,:), cons(i+1,:), adiabatic_index, alpha, flux(i,:))
    If(stage == 1) Then
      g_dt = Min(g_dt, cfl * (grid(i) - grid(i-1)) / alpha)
    End If
  End Do
  If(stage == 1 .And. g_t + g_dt > t_end) Then
    g_dt = t_end - g_t
  End If
End Subroutine CompFlux

Subroutine ComputError
  Use GlobalData
  Use Euler1d
  Implicit None
  Integer :: i
  Double Precision :: error(3)
  Double precision :: exact(nx_beg:nx_end,NCOMP)
  Double precision :: t, tmp
  
  t=g_t
  exact=0
  Do i = nx_beg, nx_end
    Call GuassIntegration(grid(i-1), grid(i), t, exact(i,:))
    ! Write(*,*) "i=",i,"exact=",exact(i,Rho)
  End Do
  error = 0
  ! Open(1, File="error.dat")
  Do i = nx_beg, nx_end
    tmp = Abs(exact(i,RHO) - prim(i,RHO))
    ! Write(*,*) "i=",i,"tmp",tmp,"h",(grid(i)-grid(i-1))
    error(1) = error(1) + (grid(i)-grid(i-1))*tmp
    error(2) = error(2) + (grid(i)-grid(i-1))*tmp*tmp
    If(error(3) < tmp) Then
      error(3) = tmp
    End If
  End Do
  error(2) = Sqrt(error(2))
  Write(*,*) nx, error
  ! Close(1)
End Subroutine ComputError

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

Subroutine rp1d(N)
  Use GlobalData
  Use TimeInfo
  Implicit None
  Integer :: i
  Integer :: N

  ! Call PrintTime
  Call ReadParameter
  nx=N
  Call InitMesh
  Call Init

  g_t=0
  Do While(g_t < t_end - 1E-14)
    Call RK3
    ! Call ForwardEuler
    g_t = g_t + g_dt
    g_step = g_step + 1
    ! Print *, "step = ", g_step, ", t = ", g_t, ", dt = ", g_dt
  End Do

  Call ComputError
  Call PrintTecplot
  Call DestroyMesh
End Subroutine

program main
  implicit none
  real*8 :: x
  Integer :: N
  N = 20
  Do While(N <= 640)
    Call rp1d(N)
    N = N*2
  End Do
end program main

