!!!#############################################################################
!!! NewK
!!!   Transform one k point [ks1,ks2,ks3] (fractional coordinates) in supercell
!!!   BZ into 'vscale' new k points [kp1,kp2,kp3], ... in primitive cell BZ
!!!   (fractional coordinates)
!!!Input:    ks1, ks2, ks3, vscale, Dp2s, toldk
!!!Output:   kp(3,vscale)
!!!#############################################################################

SUBROUTINE NewK (ks, vscale, Dp2s, toldk, &! <-- vars in
      kp) ! --> vars out
use util, only: inverse3x3
implicit none
! external vars
double precision, intent(in) :: &
   ks(3), &! k point in supercell
   Dp2s(3,3), &! primit. -> superc. transform matrix (real space)
   toldk ! tolerance for finding unique k points
integer, intent(in) :: &
   vscale ! primit. -> superc. volume scale (real space)
double precision, intent(out) :: &
   kp(3,vscale) ! unique set of new (unfolded) k points in the primitive BZ
! internal vars
double precision :: &
   kptmp(3,(vscale+1)**3), &! intermediate (non unique new k points) in 
                             ! the primitive BZ
   Ds2p(3,3), &! matrix to tranf. direct supercell to primitive vectors
   ksG(3), &! ks + G
   kptmpi(3) ! intermediate (non unique new k point)
integer :: &
   i, j, l, nunique, &! counter
   G1, G2, G3 ! G vector, e.g. [0 2 3] in primitive BZ in units of its 
              ! reciprocal lattice vectors
logical :: &
   unique ! used to identify unique k points
logical, save :: &
   writeverbose1 = .true., &
   writeverbose2 = .true.

! construct Ds2p matrix
call inverse3x3(Dp2s, Ds2p)
if (writeverbose1) then
   write (*,'(A)') 'Ds2p = inv(Dp2s) matrix is successfully generated. '//&
      'It will be used to convert k points from supercell BZ into '//&
      'primitive BZ (k_s -> k_p) '
   write (*,'(A,3(F8.4,X),A)') '       | ', Ds2p(1,:), '|'
   write (*,'(A,3(F8.4,X),A)') 'Ds2p = | ', Ds2p(2,:), '|'
   write (*,'(A,3(F8.4,X),A)') '       | ', Ds2p(3,:), '|'
   write (*,'(A,3(F8.4,X),A)') 'k_s = { ', ks, '} in units '//&
   'of the supercell reciprocal lattice vectors'
endif

! loop over G vectors of the supercell BZ to create (vscale+1)**3 possible 
! unfolded k points and transform them into the primitive BZ
i = 0 ! counter
do G1=0, vscale
   do G2=0, vscale
      do G3=0, vscale
         i = i + 1
         ksG = ks + (/G1, G2, G3/)
         kptmpi = MATMUL(ksG, Ds2p) ! convert supercell -> primitive basis
         if (writeverbose1) then
            write (*,'(A,3(I0,X),A)') 'G = { ', G1, G2, G3, '}'
            write (*,'(A,3(F8.4,X),A)') 'k_s + G = { ', ksG, '}'
            write (*,'(A)') 'k_p = (k_s + G) * [Ds2p]'
            write (*,'(A,3(F8.4,X),A)') 'k_p = { ', &
               kptmpi, '} in units of the '//&
               'primitive reciprocal lattice vectors'
         endif
         ! bring k points into the range [0,1)
         do j=1,3
            ! take mod(x,1)
            kptmpi(j) = MODULO( kptmpi(j), DBLE(1))
            if (1-kptmpi(j) .lt. toldk) then ! 0.99999 -> 1 -> 0
               kptmpi(j) = DBLE(0)
            endif
         enddo ! j
         if (writeverbose1) then
            write (*,'(A,3(F8.4,X),A)') 'k_p = { ', &
               kptmpi, '} after MOD(kp,1) operator '//&
               '(needed for large G to bring k values in the range [0,1))'
            write (*,'(A,3(I0,X),A)') 'The presess is repeated '//&
               '(not shown here) with various G vectors up to { ', &
               vscale, vscale, vscale,'}. Only unique k_p vectors are collected.'
            write (*,'(A,I0)') 'The number of unique unfolded k points '//&
               'k_p should be strictly equal to vscale = ', vscale
            writeverbose1 = .false. ! stop writing details
         endif
         kptmp(:,i) = kptmpi ! store the new (unfolded) k point
      enddo ! G3
   enddo ! G2
enddo ! G1

! take (vscale+1)**3 possible unfolded k points and find only unique ones
kp(:,1) = kptmp(:,1) ! ferst is surely unique
nunique = 1 ! counter for unique points
do i=2, (vscale+1)**3
   unique = .true.
   do j=1, nunique
      if ( ALL(ABS(kptmp(:,i)- kp(:,j)) < toldk) ) then ! point is not unique
         unique = .false.
         exit ! loop
      endif
   enddo ! j
   ! if it makes to this point without detecting similar entries, then 
   ! the point is unique
   if (unique) then
      nunique = nunique + 1
      if (nunique .gt. vscale) then
         write(*,*) 'Here is a list of non-unique unfolded k points'
         do l=1,(vscale+1)**3
            write(*,'(A,I0,A,3(F8.4,X))') 'k_p(',l,')=', kptmp(:,l)
         enddo
         write(*,*) 'Here is a list of unique unfolded k points'
         do l=1,nunique-1
            write(*,'(A,I0,A,3(F8.4,X))') 'k_p(',l,')=', kp(:,l)
         enddo
         write(*,*) 'Tolerance used toldk = ', toldk
         write(*,'(A,3(F8.4,X),A,I0,A,I0)') &
            'ERROR in NewK: The code found an additional unique '//&
            'k point = (', kptmp(:,i), ') which will bring their '//&
            'total number to ', nunique, &
            ', which is greater than the volume scale when constructed the '//&
            'supercell vscale = ', vscale
         ERROR STOP 1
      endif
      kp(:,nunique) = kptmp(:,i)
   endif
enddo ! i

! checkpoint: the number of new unique unfolded k points should be _exactly_
! equal to the volume scale from primitive to supercell (real space)
if ( nunique .ne. vscale) then
   write(*,*) 'Here is a list of non-unique unfolded k points'
   do l=1,(vscale+1)**3
      write(*,'(A,I0,A,3(F8.4,X))') 'k_p(',l,')=', kptmp(:,l)
   enddo
   write(*,*) 'Here is a list of unique unfolded k points'
   do l=1,nunique
      write(*,'(A,I0,A,3(F8.4,X))') 'k_p(',l,')=', kp(:,l)
   enddo
   write(*,*) 'Tolerance used toldk = ', toldk
   write(*,*) 'ERROR in NewK: The number of unique k points = ', nunique, &
      ' is not equal to the volume scale when constructed the supercell ', &
      'vscale =', vscale
endif

if (writeverbose2) then
   write(*,'(A)') 'Here is a list of unique unfolded k points'
   do l=1,nunique
      write(*,'(A,I0,A,3(F8.4,X))') 'k_p(',l,')=', kp(:,l)
   enddo
   write (*,'(A)') 'Further detailed output will be supressed.'
   writeverbose2 = .false. ! stop writing details
endif

END SUBROUTINE NewK
