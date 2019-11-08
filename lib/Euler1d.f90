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