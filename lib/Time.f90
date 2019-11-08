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
    cons(i,:) = 1.0_8/3*cons_1(i,:) + 2.0_8/3*(cons(i,:) - g_dt*(flux(i,:) - flux(i-1,:))/(grid(i) - grid(i-1)))
    Call Con2Pri(cons(i,:), adiabatic_index, prim(i,:))
  End Do
End Subroutine RK3