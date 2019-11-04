!**
! @file rp1d.f90
! @brief A simple Fortran exercise for 1D Euler
! @author Duan Junming, duanjm@pku.edu.cn
! @version 1.0
! @date 2019-11-03
!/

Module TimeInfo
  Implicit None
  Integer time_info_spe(8)
  Character*10 time_info(3)
Contains
  Subroutine PrintTime
    Call date_and_time(time_info(1), time_info(2), time_info(3), time_info_spe)
    Print *, "Time: ", time_info(1), time_info(2)
  End Subroutine PrintTime
End Module TimeInfo

Module Euler1d
  Implicit None
  Integer :: RHO = 1, VX1 = 2, PRS = 3, NCOMP = 3
  Contains

  Subroutine Pri2Con(prim, adiabatic_index, cons)
    Double Precision :: cons(NCOMP), prim(NCOMP), adiabatic_index
    cons(RHO) = prim(RHO)
    cons(VX1) = prim(RHO) * prim(VX1)
    cons(PRS) = 0.5 * prim(RHO) * prim(VX1)**2 + prim(PRS)/(adiabatic_index - 1.0)
  End Subroutine Pri2Con

  Subroutine Con2Pri(cons, adiabatic_index, prim)
    Double Precision :: cons(NCOMP), prim(NCOMP), adiabatic_index
    prim(RHO) = cons(RHO)
    prim(VX1) = cons(VX1) / cons(RHO)
    prim(PRS) = (cons(PRS) - 0.5 * cons(VX1)**2 / cons(RHO)) * (adiabatic_index - 1.0)
  End Subroutine Con2Pri

  Subroutine MaxLambda(prim, cons, adiabatic_index, lambda)
    Double Precision :: prim(NCOMP), cons(NCOMP), lambda, adiabatic_index, c
    c = Sqrt(adiabatic_index*prim(PRS)/prim(RHO))
    lambda = Abs(prim(VX1)) + c
  End Subroutine MaxLambda

  Subroutine NormalFlux(prim, cons, f)
    Double Precision :: prim(NCOMP), cons(NCOMP), f(NCOMP)
    f(RHO) = cons(VX1)
    f(VX1) = cons(VX1)**2 / cons(RHO) + prim(PRS)
    f(PRS) = (cons(PRS) + prim(PRS)) * prim(VX1)
  End Subroutine NormalFlux

  Subroutine LLf(ul, ur, adiabatic_index, alpha, f)
    Double Precision :: ul(NCOMP), ur(NCOMP), adiabatic_index, f(NCOMP)
    Double Precision :: vl(NCOMP), vr(NCOMP), laml, lamr, alpha, fl(NCOMP), fr(NCOMP)
    Call Con2Pri(ul, adiabatic_index, vl)
    Call Con2Pri(ur, adiabatic_index, vr)
    Call MaxLambda(vl, ul, adiabatic_index, laml)
    Call MaxLambda(vr, ur, adiabatic_index, lamr)
    Call NormalFlux(vl, ul, fl)
    Call NormalFlux(vr, ur, fr)
    alpha = Max(laml, lamr)
    f = 0.5*(fl+fr) - 0.5*alpha*(ur-ul)
  End Subroutine LLf

End Module Euler1d

Module Reconstruction
  Implicit None
  Double Precision :: thirdteen_12 = 13.0/12
  Double Precision :: eps = 1E-6
Contains
  Subroutine Weno5(cell_av, limit, ncomp)
    Integer :: ncomp, i
    Double Precision :: cell_av(5,ncomp), limit(ncomp)
    Double Precision :: a0, a1, a2, b0, b1, b2, sum_a
    Do i = 1, ncomp
      a0 = cell_av(1,i) - 2.0*cell_av(2,i) +     cell_av(3,i)
      a1 = cell_av(1,i) - 4.0*cell_av(2,i) + 3.0*cell_av(3,i)
      b0 = thirdteen_12*a0**2 + 0.25*a1**2

      a0 = cell_av(2,i) - 2.0*cell_av(3,i) +     cell_av(4,i)
      a1 = cell_av(2,i)                    -     cell_av(4,i)
      b1 = thirdteen_12*a0**2 + 0.25*a1**2

      a0 =     cell_av(3,i) - 2.0*cell_av(4,i) + cell_av(5,i)
      a1 = 3.0*cell_av(3,i) - 4.0*cell_av(4,i) + cell_av(5,i)
      b2 = thirdteen_12*a0**2 + 0.25*a1**2

      a0 = 0.1/(b0+eps)**2
      a1 = 0.6/(b1+eps)**2
      a2 = 0.3/(b2+eps)**2

      sum_a = 1.0/(a0 + a1 + a2)
      a0 = a0*sum_a
      a1 = a1*sum_a
      a2 = a2*sum_a

      b0 = 2.0*cell_av(1,i) - 7.0*cell_av(2,i) + 11.0*cell_av(3,i)
      b1 =    -cell_av(2,i) + 5.0*cell_av(3,i) +  2.0*cell_av(4,i)
      b2 = 2.0*cell_av(3,i) + 5.0*cell_av(4,i) -      cell_av(5,i)

      limit(i) = (a0*b0 + a1*b1 + a2*b2)/6.0
    End Do
  End Subroutine Weno5
End Module

Module GlobalData
  Use Euler1d
  Implicit None
  Integer :: nx, nx_ghost, nx_beg, nx_end, nx_glbbeg, nx_glbend, xbeg_bc_type, xend_bc_type
  Double Precision :: x_beg, x_end, x_glbbeg, x_glbend, g_t, g_dt, g_step, t_end
  Double Precision :: cfl, adiabatic_index, pi
  Double Precision, Pointer :: grid(:), cons(:,:), prim(:,:), rec_cons_m(:,:), rec_cons_p(:,:), flux(:,:), cons_1(:,:), cons_2(:,:)
Contains

  Subroutine ReadParameter
    pi = 4*Atan(1.d0)
    Print *, "Read parameters ..."
    Open(1, File="para.ini")
    Read(1, *) nx;              Write(*, *) "nx = ", nx
    Read(1, *) x_beg;           Write(*, *) "x_beg = ", x_beg
    Read(1, *) x_end;           Write(*, *) "x_end = ", x_end
    Read(1, *) nx_ghost;        Write(*, *) "nx_ghost = ", nx_ghost
    Read(1, *) adiabatic_index; Write(*, *) "adiabatic_index = ", adiabatic_index
    Read(1, *) cfl;             Write(*, *) "cfl = ", cfl
    Read(1, *) t_end;           Write(*, *) "t_end = ", t_end
    Read(1, *) xbeg_bc_type;    Write(*, *) "xbeg_bc_type = ", xbeg_bc_type
    Read(1, *) xend_bc_type;    Write(*, *) "xend_bc_type = ", xend_bc_type
!100 format(a20, 10i8)
!101 format(a20, 10f8)
!102 format(a20, 10a8)
    Close(1)
  End Subroutine ReadParameter

  Subroutine InitMesh
    Integer :: i
    Double Precision :: hx
    nx_beg = nx_ghost+1
    nx_end = nx + nx_ghost
    nx_glbbeg = 1
    nx_glbend = nx + 2*nx_ghost
    Allocate(grid(nx_glbbeg-1 : nx_glbend))
    Allocate(cons(nx_glbbeg : nx_glbend, NCOMP))
    Allocate(prim(nx_glbbeg : nx_glbend, NCOMP))
    Allocate(rec_cons_m(nx_glbbeg-1 : nx_glbend, NCOMP))
    Allocate(rec_cons_p(nx_glbbeg-1 : nx_glbend, NCOMP))
    Allocate(flux(nx_glbbeg-1 : nx_glbend, NCOMP))
    Allocate(cons_1(nx_glbbeg : nx_glbend, NCOMP))
    Allocate(cons_2(nx_glbbeg : nx_glbend, NCOMP))
    hx = (x_end - x_beg)/nx
    Do i = nx_glbbeg-1, nx_glbend
      grid(i) = x_beg + (i-nx_beg+1)*hx
    EndDo
  End Subroutine InitMesh

  Subroutine DestroyMesh
    Deallocate(grid)
    Deallocate(cons)
    Deallocate(prim)
    Deallocate(rec_cons_m)
    Deallocate(rec_cons_p)
    Deallocate(flux)
    Deallocate(cons_1)
    Deallocate(cons_2)
  End Subroutine DestroyMesh

End Module GlobalData

Subroutine F0(x, pri)
  Use GlobalData
  Use Euler1d
  Implicit None
  Double Precision :: x, pri(NCOMP)
  !pri(RHO) = 1 + 0.5*Sin(2*pi*x)
  !pri(VX1) = 0.5
  !pri(PRS) = 1.0
  If(x < 0.5) then
    pri(RHO) = 1
    pri(VX1) = 0
    pri(PRS) = 1
  Else
    pri(RHO) = 0.125
    pri(VX1) = 0
    pri(PRS) = 0.1
  End If
End Subroutine F0

Subroutine Init
  Use GlobalData
  Implicit None
  Integer :: i
  Do i = nx_beg, nx_end
    Call F0(0.5*(grid(i-1)+grid(i)), prim(i,:))
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

Subroutine ForwardEuler
  Use GlobalData
  Implicit None
  Integer :: i
  !1st stage
  Call Boundary
  Call Reconstruct
  Call CompFlux(1)
  Do i = nx_beg, nx_end
    cons(i,:) = cons(i,:) - g_dt*(flux(i,:) - flux(i-1,:))/(grid(i) - grid(i-1))
    Call Con2Pri(cons(i,:), adiabatic_index, prim(i,:))
  End Do
End Subroutine ForwardEuler

Subroutine RK3
  Use GlobalData
  Implicit None
  Integer :: i
  cons_1 = cons
  !1st stage
  Call Boundary
  Call Reconstruct
  Call CompFlux(1)
  Do i = nx_beg, nx_end
    cons(i,:) = cons_1(i,:) - g_dt*(flux(i,:) - flux(i-1,:))/(grid(i) - grid(i-1))
    Call Con2Pri(cons(i,:), adiabatic_index, prim(i,:))
  End Do
  !2nd stage
  Call Boundary
  Call Reconstruct
  Call CompFlux(2)
  Do i = nx_beg, nx_end
    cons(i,:) = 0.75*cons_1(i,:) + 0.25*(cons(i,:) - g_dt*(flux(i,:) - flux(i-1,:))/(grid(i) - grid(i-1)))
    Call Con2Pri(cons(i,:), adiabatic_index, prim(i,:))
  End Do
  !3rd stage
  Call Boundary
  Call Reconstruct
  Call CompFlux(3)
  Do i = nx_beg, nx_end
    cons(i,:) = 1.0/3*cons_1(i,:) + 2.0/3*(cons(i,:) - g_dt*(flux(i,:) - flux(i-1,:))/(grid(i) - grid(i-1)))
    Call Con2Pri(cons(i,:), adiabatic_index, prim(i,:))
  End Do
End Subroutine RK3

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

  Do While(g_t < t_end - 1E-14)
    Call RK3
    !Call ForwardEuler
    g_t = g_t + g_dt
    g_step = g_step + 1
    Print *, "step = ", g_step, ", t = ", g_t, ", dt = ", g_dt
  End Do

  !Call ComputError
  Call PrintTecplot
  Call DestroyMesh
End Program rp1d

