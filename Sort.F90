!!!#################################################################################
!!! Sort
!!!     Sorts energy coefficients (complex number form) into appropriate groups 
!!!	according to relative location to Gamma point. Then calculates the weight 
!!!	of each group.
!!!#################################################################################
SUBROUTINE Sort(ks, kp, vscale, Dp2s, toldk, G, &! <-- args in
    pwcoeffr, pwcoeffz, &! <-- opt. args in
    NV, Orb, &! <-- args in
    w) ! --> args out
use util, only: inverse3x3
implicit none
! external vars
integer, intent(in) :: &
    vscale, &! primit. -> superc. volume scale (real space)
    NV, &! length of PW coefficient vector
    Orb, &! number of local orbitals in the vector (at the end)
    G(3,NV) ! matrix of PW lattice vectors
double precision, intent(in) :: &
    ks(3), &! k point in supercell
    kp(3,vscale), &! ks point unfolded to primitive BZ
    Dp2s(3,3), &! primit. -> superc. transform matrix (real space)
    toldk ! tolerance for finding unique k points
double precision, intent(in), optional :: &
    pwcoeffr(NV) ! plane wave coefficients r/z (real/complex)
double complex, intent(in), optional :: &
    pwcoeffz(NV) ! plane wave coefficients r/z (real/complex)
double precision, intent(out) :: &
    w(vscale) ! weights after unfolding (0-1)
! internal vars
double precision :: &
    Ds2p(3,3), &! matrix to transf. direct supercell to primitive vectors
    ksG(3), &! ks + G
    kptmpi(3) ! intermediate (non unique new k point)
integer :: &
    i, j, &! counter
    countw(vscale) ! count entries into weight bins
logical :: &
    matchfound ! used to identify when a k point match found in a group 'kp'

! check if one of optional variables is defined
if (present(pwcoeffr) .and. present(pwcoeffz)) then
    ! both arguments defined
    write(*,*) 'ERROR in subroutine Sort: both pwcoeffr and pwcoeffz ',&
        'arguments are defined.'
    ERROR STOP 1
else if ( .not.(present(pwcoeffr)) .and. .not.(present(pwcoeffz))) then
    ! none of two arguments defined
    write(*,*) 'ERROR in subroutine Sort: none of pwcoeffr and pwcoeffz ',&
        'arguments are defined.'
    ERROR STOP 1
endif

! construct Ds2p matrix
call inverse3x3(Dp2s, Ds2p)

! initialize weights and its counter
w = DBLE(0)
countw = 0

! sorts the PW coefficients into bins of weights
do i=1, (NV-Orb)
    ksG = ks + (/G(1,i), G(2,i), G(3,i)/)
    kptmpi = MATMUL(ksG, Ds2p) ! convert supercell -> primitive basis
    ! bring k points into the range [0,1)
    do j=1,3
        kptmpi(j) = MODULO( kptmpi(j), DBLE(1))
        if (1-kptmpi(j) .lt. toldk) then ! 0.99999 -> 1 -> 0
            kptmpi(j) = DBLE(0)
        endif
    enddo ! j
    ! compare the kptmpi point with the list 'kp' and find which group 
    ! it belongs
    matchfound = .false.
    do j=1, vscale
        if ( ALL(ABS(kptmpi- kp(:,j)) < toldk) ) then
            ! point matched to one on the 'kp' list
            matchfound = .true.
            ! store absolute value squared of the PW coefficient in
            ! appropriate bin 'w'
            if (present(pwcoeffr)) then
                ! PW coeff. are real
                w(j) = w(j) + pwcoeffr(i)*pwcoeffr(i)
            else if (present(pwcoeffz)) then
                ! PW coeff. are complex
                w(j) = w(j) + DBLE(pwcoeffz(i)*CONJG(pwcoeffz(i)))
            endif
            countw(j) = countw(j) + 1 ! update counter
            exit ! loop
        endif
    enddo ! j (kp)
    if (.not.matchfound) then
        write(*,*) 'ERROR in subroutine Sort: unable to match k point (', &
            kptmpi, ') to a point in the group kp listed below'
        do j=1, vscale
            write(*,*) 'kp(', j, ') = ', kp(:,j)
        enddo
        ERROR STOP 1
    endif
enddo ! i (PW)

! check if all bins got weights
do i=1, vscale
    if (countw(j) .eq. 0) then
        write(*,*) 'ERROR in subroutine Sort: unable to populate k point (', &
            kp(:,i), ') from the group kp listed below. This is unlikely.'
        do j=1, vscale
            write(*,*) 'kp(', j, ') = ', kp(:,j)
        enddo
        ERROR STOP 1
    endif
enddo

! normalize sum weights to 1
w = w/SUM(w)

END SUBROUTINE Sort
