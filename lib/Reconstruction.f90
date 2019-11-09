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