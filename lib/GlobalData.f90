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
    ! Print *, "Read parameters ..."
    Open(1, File="../src/para.ini")
    Read(1, *) nx;              !Write(*, "('nx =',1XI6)") nx
    Read(1, *) x_beg;           !Write(*, "('x_beg =',1XF10.5)") x_beg
    Read(1, *) x_end;           !Write(*, "('x_end =',1XF10.5)") x_end
    Read(1, *) nx_ghost;        !Write(*, "('nx_ghost =',1XI4)") nx_ghost
    Read(1, *) adiabatic_index; !Write(*, "('adiabatic_index =',1XF10.5)") adiabatic_index
    Read(1, *) cfl;             !Write(*, "('cfl =',1XF10.5)")  cfl
    Read(1, *) t_end;           !Write(*, "('t_end =',1XF10.5)") t_end
    Read(1, *) xbeg_bc_type;    !Write(*, "('xbeg_bc_type =',1XI4)") xbeg_bc_type
    Read(1, *) xend_bc_type;    !Write(*, "('xend_bc_type =',1XI4)") xend_bc_type
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
      ! Write(*,*) "grid(i)=",grid(i),i
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