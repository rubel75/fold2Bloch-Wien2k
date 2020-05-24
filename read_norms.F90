SUBROUTINE read_norms(unitnorm, unitnorm2, ne, &! <-- in
    sonorm, sonorm2) ! --> out

implicit none
INTEGER, INTENT(in) :: &
    unitnorm, unitnorm2, ne
DOUBLE PRECISION, INTENT(out) :: &
    sonorm(ne), sonorm2(ne)
! Internal
INTEGER :: &
    i

!! Read norms from files for all bands
!! (format from $WIENROOT/SRC_lapwso/kptout.F)

read(unitnorm,'(4e20.12)') (sonorm(i), i=1,ne) ! loop over all eigenstates
read(unitnorm2,'(4e20.12)') (sonorm2(i), i=1,ne)

RETURN
END SUBROUTINE read_norms
