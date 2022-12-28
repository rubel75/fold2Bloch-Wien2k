! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! fold2Bloch.f90		
! 
!    Unfold the data in case.vector[_1] file and compute Bloch spectral weight
!
!    Usage:
!       # get help
!       fold2Bloch -h
!
!       # real calculation (inversion symm.)
!       fold2Bloch -r case.vector[_1] ...
!           "'P11 P12 P13:P21 P22 P23:P31 P32 P33'"
!
!       # complex calc. (no inv. symm.)
!       fold2Bloch -c case.vector[_1] ...
!           "'P11 P12 P13:P21 P22 P23:P31 P32 P33'"
!
!       # complex calc. implied
!       fold2Bloch case.vector[_1] ...
!           "'P11 P12 P13:P21 P22 P23:P31 P32 P33'"
!
!       # spin-orbit coupling without -sp
!       fold2Bloch -so case.vectorso[_1] case.vectorsodn[_1] ...
!           case.normsoup[_1] case.normsodn[_1] ...
!           "'P11 P12 P13:P21 P22 P23:P31 P32 P33'" 
!
!       # spin-orbit coupling with -sp
!       fold2Bloch -so case.vectorsoup[_1] case.vectorsodn[_1] ...
!           case.normsoup[_1] case.normsodn[_1] ...
!           "'P11 P12 P13:P21 P22 P23:P31 P32 P33'" 
!
!   Compilation:
!       see https://github.com/rubel75/fold2Bloch-Wien2k/wiki/Installation
!
! Copyright 2022 Oleg Rubel, Anton Bokhanchuk, Elias Assmann
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                               MAIN
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PROGRAM fold2Bloch
use structmod, only: struct, struct_read

implicit none

INTERFACE
    SUBROUTINE Sort(ks, kp, vscale, Dp2s, toldk, G, &! <-- args in
                    pwcoeffr, pwcoeffz, &! <-- opt. args in
                    NV, Orb, &! <-- args in
                    w) ! --> args out
        integer, intent(in) :: &
            vscale, &! primit. -> superc. volume scale (real space)
            Dp2s(3,3), &! primit. -> superc. transform matrix (real space)
            NV, &! length of PW coefficient vector
            Orb, &! number of local orbitals in the vector (at the end)
            G(3,NV) ! matrix of PW lattice vectors
        double precision, intent(in) :: &
            ks(3), &! k point in supercell
            kp(3,vscale), &! ks point unfolded to primitive BZ
            toldk ! tolerance for finding unique k points
        double precision, intent(in), optional :: &
            pwcoeffr(NV) ! plane wave coefficients r/z (real/complex)
        double complex, intent(in), optional :: &
            pwcoeffz(NV) ! plane wave coefficients r/z (real/complex)
        double precision, intent(out) :: &
            w(vscale) ! weights after unfolding (0-1)
    END SUBROUTINE Sort
END INTERFACE

INTEGER :: &
    nargs, & !command line input argument counter
    jatom,i,k,jl,jj,l,q, &
    i2, & ! spinor dn component
    lmax, lomax, nloat, &
    unitklist, unitstruct, unitvector, unitout, & ! I/O unit numbers
    unitvector2, & ! spinor dn component
    unitnorm, unitnorm2, & ! fine I/O unit for spinor norm components
    line_count, &
    NE, NV, ORB, &
    NE2, NV2, ORB2, & ! spinor dn component
    ios, iocplx, nkcount, nkpoints, &
    ios2, &! spinor dn component
    vscale, &! volume scale primitive -> supercell
    Dp2s(3,3) ! matrix to convert primitive cell to a supercell 
              ! in direct (real) space
INTEGER, ALLOCATABLE :: &
    Vector(:,:), &
    Vector2(:,:) ! spinor dn component
DOUBLE PRECISION :: &
    KX, KY, KZ, WEIGHT, & ! k-point coordinates and weight
    KX2, KY2, KZ2, WEIGHT2 ! spinor dn component
DOUBLE PRECISION, ALLOCATABLE :: &
    EIGVAL(:), &
    EIGVAL2(:), & ! spinor dn component
    E(:,:),ELO(:,:,:), &
    E2(:,:),ELO2(:,:,:), & ! spinor dn component
    Coef(:), & ! PW coefficients
    sonorm(:), sonorm2(:), &! norm of spin up/dn components for a spinor function
    Weights(:), &
    Weights2(:), &! spinor dn component
    NKVal(:,:) ! k point unfolded from supercell into the primitive cell
DOUBLE COMPLEX, ALLOCATABLE :: &
    CoefC(:), & ! PW coefficients (complex)
    CoefC2(:) ! spinor dn component
CHARACTER :: &
    casename*100, vectorname*100, normname*100, &
    vectorname2*100, normname2*100, & ! spinor dn component
    folds*400, &
    KNAME*10, & ! k-point name
    KNAME2*10, & ! spinor dn component
    endout*5, &
    bar
CHARACTER(400), ALLOCATABLE :: &
    args(:)
LOGICAL :: &
    usecomplex, dir, &
    lso ! spin-orbit coupling

TYPE(struct) stru ! structure specific info (see unil.F file)

!! Set I/O unit numbers (SOC files will be defined later)

unitvector = 11 ! vector file
unitout = 12 ! output file
unitklist = 13 ! case.klist
unitstruct = 14 ! case struct

!parameters
LMAX = 13
LOMAX = 3
nloat = 3

write(*,*) '**************************'
write(*,*) '**     fold2Bloch       **'
write(*,*) '** version Jul 12, 2022 **'
write(*,*) '**************************'

!! command line arguments read-in
  
nargs=iargc() ! number of input arguments
allocate (args(nargs))
write(*,'(A,I0,A)') 'Detected ', nargs,' input arguments'
if (nargs==1) then ! get help
    CALL GETARG(1,args(1))
    if ( (args(1)=='-h') .or. (args(1)=='--help') ) then
        GOTO 910 ! print help and exit
    else ! unknown
        GOTO 912 ! print error, usage options, and exit
    endif
elseif (nargs==2) then ! 2 input arguments (vector file and folds)
    CALL GETARG(1,args(1))
    CALL GETARG(2,args(2))
    read(args(1),*) vectorname ! 1st argument is the vector file name
    read(args(2),*) folds ! 2nd argument is Dp2s matrix
    usecomplex = .true. ! complex calculation implied
    lso = .false. ! no spin-orbit coupling
elseif (nargs==3) then ! 3 arguments (-r/-c, vector file, folds)
    CALL GETARG(1,args(1))
    CALL GETARG(2,args(2))
    CALL GETARG(3,args(3))
    if (args(1)=='-r') then
        usecomplex = .false.
    elseif (args(1)=='-c') then
        usecomplex = .true.
    else ! imposible
        GOTO 912 ! print error, usage options, and exit
    endif
    read(args(2),*) vectorname ! 2nd argument is the vector file name
    read(args(3),*) folds ! 3rd argument is Dp2s matrix
    lso = .false. ! no spin-orbit coupling
elseif (nargs==6) then ! 6 arguments (spin-orbit coupling)
    do i=1,nargs
        CALL GETARG(i,args(i))
    end do
    if (args(1)=='-so') then
        lso = .true. ! enable SOC
        write(*,'(A)') 'Spin-orbit calculation mode'
    else ! imposible
        GOTO 912 ! print error, usage options, and exit
    endif
    read(args(2),*) vectorname ! should be case.vectorsoup[_X]
    read(args(3),*) vectorname2 ! should be case.vectorsodn[_X]
    read(args(4),*) normname ! should be case.normsoup[_X]
    read(args(5),*) normname2 ! should be case.normsodn[_X]
    read(args(6),*) folds ! 6th argument is Dp2s matrix
    usecomplex = .true. ! complex calcilation implied
else ! impossible
    write(*,'(A,I0,A)') 'Detected ', nargs,' input arguments'
    GOTO 912 ! print error, usage options, and exit
endif

!! Check if vector file(s) and norm files (SOC only) are present

write(*,'(A)') 'Checking prerequisite files...' 
inquire(file=vectorname, exist=dir)
if (.not.(dir)) then
    write(*,'(A)') trim(vectorname)//&
        '  vector file does not exist or invalid file name entered'
    GOTO 912 ! print error, usage options, and exit
else
    write(*,'(A)') trim(vectorname)//&
        '  vector file found'
endif
if (lso) then ! SOC only
    inquire(file=vectorname2, exist=dir) ! case.vectorsodn[_X]
    if (.not.(dir)) then
        write(*,'(A)') trim(vectorname2)//&
            '  vector file does not exist or invalid file name entered'
        GOTO 912 ! print error, usage options, and exit
    else
        write(*,'(A)') trim(vectorname2)//&
            '  vector file found'
    endif
    inquire(file=normname, exist=dir) ! case.normsoup[_X]
    if (.not.(dir)) then
        write(*,'(A)') trim(normname)//&
            '  norm file does not exist or invalid file name entered'
        GOTO 912 ! print error, usage options, and exit
    else
        write(*,'(A)') trim(normname)//&
            '  norm file found'
    endif
    inquire(file=normname2, exist=dir) ! case.normsodn[_X]
    if (.not.(dir)) then
        write(*,'(A)') trim(normname2)//&
            '  norm file does not exist or invalid file name entered'
        GOTO 912 ! print error, usage options, and exit
    else
        write(*,'(A)') trim(normname2)//&
            '  norm file found'
    endif
endif

!! Check if the number of folds argument is in correct format

i=INDEX(folds, ':') ! return the position of : in the string
if (i.eq.0) then
    write(*,*) 'folds = ', TRIM(folds)
    write(*,*) 'Unknown number of folds. See below or type: "fold2Bloch -h" '//&
        'for more information.'
    GOTO 912 ! print error, usage options, and exit
endif
read (folds(1:i-1), *, iostat=ios) Dp2s(1,1), Dp2s(1,2), Dp2s(1,3)
if (ios.ne.0) then
    write (*,*) 'Unable to read the Dp2s(1,:) matrix'
    write (*,*) 'Relevant input line = ', trim(folds)
    write (*,*) 'Parsed input section = ', folds(1:i-1)
    write (*,*) 'while expected 3 numerical values separated by space.'
    GOTO 912 ! print error, usage options, and exit
endif
folds=folds(i+1:)
i=INDEX(folds, ':')
if (i.eq.0) then
    write(*,*) 'folds = ', TRIM(folds)
    write(*,*) 'Unknown number of folds. See below or type: "fold2Bloch -h" '//&
        'for more information.'
    GOTO 912 ! print error, usage options, and exit
endif
read (folds(1:i-1),*, iostat=ios) Dp2s(2,1), Dp2s(2,2), Dp2s(2,3)
if (ios.ne.0) then
    write (*,*) 'Unable to read the Dp2s(2,:) matrix'
    write (*,*) 'Parsed input section = ', folds(1:i-1)
    write (*,*) 'while expected 3 numerical values separated by space.'
    GOTO 912 ! print error, usage options, and exit
endif
read(folds(i+1:),*, iostat=ios) Dp2s(3,1), Dp2s(3,2), Dp2s(3,3)
if (ios.ne.0) then
    write (*,*) 'Unable to read the Dp2s(3,:) matrix'
    write (*,*) 'Parsed input section = ', folds(i+1:)
    write (*,*) 'while expected 3 numerical values separated by space.'
    GOTO 912 ! print error, usage options, and exit
endif
! Evaluate the volume change prim -> supercell lattice
vscale = Dp2s(1,1)*Dp2s(2,2)*Dp2s(3,3) + Dp2s(1,2)*Dp2s(2,3)*Dp2s(3,1) &
    +Dp2s(1,3)*Dp2s(2,1)*Dp2s(3,2) - Dp2s(3,1)*Dp2s(2,2)*Dp2s(1,3) &
    -Dp2s(1,1)*Dp2s(3,2)*Dp2s(2,3) - Dp2s(2,1)*Dp2s(1,2)*Dp2s(3,3) ! det(Dp2s)
! Check volume chenge (not less than 0 and integer)
if (vscale.le.0) then
    write (*,*) 'The input matrix Dp2s that transforms the real space'
    write (*,*) 'primitive lattice vectors to a supercell leads to'
    write (*,*) 'the following volume change: ', vscale
    write (*,*) 'Parsed input section = ', folds(i+1:)
    write (*,*) 'FYI: Dp2s = ', Dp2s
    write (*,*) 'This means that the Dp2s matrix is not positively defined.'
    write (*,*) 'Also the volume change should be an integer number.'
    write (*,*) 'Please review your Dp2s input and rerun fold2Bloch.'
    GOTO 912 ! print error, usage options, and exit
else
    write (*,'(A)') 'Dp2s input matrix is successfully parsed.'
    write (*,'(A,3(I5,X),A)') '       | ', Dp2s(1,:), '|'
    write (*,'(A,3(I5,X),A)') 'Dp2s = | ', Dp2s(2,:), '|'
    write (*,'(A,3(I5,X),A)') '       | ', Dp2s(3,:), '|'
    write (*,'(A,I0)') 'The primitive to supercell volume scale is: ', &
        vscale
endif

!! Get case name

i=INDEX(vectorname, '.')
read (vectorname(1:i-1), *) casename

!! Get number of atoms from case.struct

inquire(file=trim(casename)//'.struct', exist=dir)
if (.not.(dir)) then ! case.struct is not present
    write(*,*) "     ", trim(casename)//'.struct file does not exist'
    ERROR STOP 1
else
    open(unit=unitstruct,file=trim(casename)//'.struct',status='old')
    ! read case.struct and return stru%nat
    call struct_read(unitstruct,&! <-- args in
        stru) ! --> args out
    close(unitstruct)
endif

!! Set output file name

dir=.TRUE.
i=1
endout=''
do while (dir)
    inquire(file=trim(casename)//'.f2b'//trim(endout), exist=dir)
    if (dir) then ! the file already exist
        write(endout, "('_', I0)") i ! append the name with next integer
        i=i+1
    end if
end do

!! Check if case.klist is available to determine number of k points for progress
!! calculation

inquire(file=trim(casename)//'.klist', exist=dir)
if (.not.(dir)) then
    write(*,*) trim("WARNING: ")//trim(casename)//'.klist could '//&
        'not be found. Progress can not be calculated!'
    write(*,*) "Continuing to process file..."
else
    open(unit=unitklist,file=trim(casename)//'.klist',status='old')
    nkpoints = line_count(unitklist) ! read the number of K points
    close(unitklist)
endif

!! Display starting message

write(*,'(A,A)')"FILE TO PROCESS: ", TRIM(vectorname)
if (lso) then
    write(*,'(A,A)')"                 ", TRIM(vectorname2)
    write(*,'(A,A)')"additional norm files: ", TRIM(normname)
    write(*,'(A,A)')"                       ", TRIM(normname2)
endif
write(*,'(A)') "/\/\/\ UNFOLDING VECTOR FILE /\/\/\"

!! Open vector file(s) and output file

open(unit=unitout,file=trim(casename)//".f2b"//trim(endout), &
    status='unknown',form='formatted')   
open(unit=unitvector, file=vectorname, status='old',form='unformatted')
if (lso) then ! SOC needs 2 vector files and 2 norm files
    unitvector2 = 21 ! vector file 2
    unitnorm = 22 ! norm file 1
    unitnorm2 = 23 ! norm file 2
    open(unit=unitvector2, file=vectorname2, status='old',form='unformatted')
    open(unit=unitnorm, file=normname, status='old',form='formatted')
    open(unit=unitnorm2, file=normname2, status='old',form='formatted')
endif

!! Read atom-specific info from the vector file
!! (dummy, values are not used to compute Bloch character)

allocate( E(LMAX,stru%nat) )
allocate( ELO(0:LOMAX,nloat,stru%nat) )
if (lso) then ! SOC
    allocate( E2(LMAX,stru%nat) )
    allocate( ELO2(0:LOMAX,nloat,stru%nat) )
endif
do jatom= 1,stru%nat
    read(unitvector) (E(jl,jatom),jl=1,LMAX)
    read(unitvector) ((ELO(jl,k,jatom),jl=0,LOMAX),k=1,nloat)
    if (lso) then ! SOC
        read(unitvector2) (E2(jl,jatom),jl=1,LMAX)
        read(unitvector2) ((ELO2(jl,k,jatom),jl=0,LOMAX),k=1,nloat)
    endif
enddo
deallocate(E,ELO) ! discard values

!! MAIN LOOP

nkcount=0
ios=0
bar='\'
do while (ios.eq.0)
    read(unitvector,end=888,err=911, iostat=ios) &
        KX, KY, KZ, KNAME, NV, NE, WEIGHT
    if (lso) then ! SOC
        read(unitvector2,end=888,err=911, iostat=ios2) &
            KX2, KY2, KZ2, KNAME2, NV2, NE2, WEIGHT2
        ! check consistency of vector files
        if ( (KX .ne. KX2) .or. (KY .ne. KY2) .or. (KZ .ne. KZ2) .or. &
                (NV .ne. NV2) .or. (NE .ne. NE2) .or. (ios .ne. ios2) ) then
            write(*,*) 'inconsistency detected between up/dn vector files'
            write(*,*) 'KX=', KX, 'KX2=', KX2
            write(*,*) 'KY=', KY, 'KY2=', KY2
            write(*,*) 'KZ=', KZ, 'KZ2=', KZ2
            write(*,*) 'NV=', NV, 'NV2=', NV2
            write(*,*) 'NE=', NE, 'NE2=', NE2
            write(*,*) 'IOSTAT specifier: ios=', ios, 'ios2=', ios2
            ERROR STOP 1
        endif
    endif
    write(*,'(a, f6.3, f6.3, f6.3)') 'Processing K-Point:', KX, KY, KZ
    nkcount=nkcount+1
    allocate( NKVal(3,vscale) )
    ! Convert 1 k point [KX, KY, KZ] in the supercell BZ into 'vscale' number 
    ! k points in the primitive BZ. The tolerance of 0.00001 is used to find 
    ! duplicate k points. 'vscale' new k points should be generated as a result 
    ! of unfolding.
    call NewK( (/KX, KY, KZ/), vscale, Dp2s, DBLE(0.00001), &! <-- args in
        NKVal) ! --> args out
    allocate( Vector(3,NV) )
    read(unitvector) (Vector(1,I),Vector(2,I),Vector(3,I),I=1,NV)
    do l=2, NV !Determine the number of local orbitals
        if ( (Vector(1,l).eq.Vector(1,1)) .and. (Vector(2,l).eq.Vector(2,1)) &
                .and. (Vector(3,l).eq.Vector(3,1)) ) then 
            ORB = NV-l+1
            exit ! loop
        endif
    enddo
    if (lso) then ! SOC
        allocate( Vector2(3,NV) )
        read(unitvector2) (Vector2(1,I),Vector2(2,I),Vector2(3,I),I=1,NV)
        do l=2, NV !Determine the number of local orbitals
            if ( (Vector2(1,l).eq.Vector2(1,1)) &
                    .and. (Vector2(2,l).eq.Vector2(2,1)) &
                    .and. (Vector2(3,l).eq.Vector2(3,1)) ) then 
                ORB2 = NV-l+1
                exit ! loop
            endif
        enddo
        ! check consistency of vector files
        if ( (ORB .ne. ORB2) ) then
            write(*,*) 'inconsistency detected between up/dn vector files'
            write(*,*) 'ORB=', ORB, 'ORB2=', ORB2
            ERROR STOP 1
        endif
    endif
    allocate(EIGVAL(NE))
    if (lso) then ! SOC
        allocate( EIGVAL2(NE), sonorm(NE), sonorm2(NE) )
        ! read norms from files for all eigenvalues
        CALL read_norms(unitnorm, unitnorm2, NE, &! <-- args in
            sonorm, sonorm2) ! ---> args out
    endif
    do jj = 1, NE
        read(unitvector) I, EIGVAL(I)
        if (i .ne. jj) then
            write(*,*) 'inconsistency in eigenvalue index'
            write(*,*) 'jj=', jj, 'I=', I
            ERROR STOP 1
        endif
        if (lso) then ! SOC
            read(unitvector2) I2, EIGVAL2(I2)
            ! check consistency of vector files
            if ( (EIGVAL(I) .ne. EIGVAL2(I2)) ) then
                write(*,*) 'inconsistency detected between up/dn vector files'
                write(*,*) 'I=', I, 'I2=', I2
                write(*,*) 'EIGVAL(I)=', EIGVAL(I), 'EIGVAL2(I2)=', EIGVAL2(I2)
                ERROR STOP 1
            endif
        endif
        allocate( Weights(vscale) )
        if (usecomplex) then
            allocate( CoefC(NV) )
            read(unitvector, ioStat=iocplx) CoefC(1:NV)
            if (iocplx.ne.0) then
                write(*,*) "ERROR: you either forgot [-r] (no complex "//&
                    "calculation) switch, or there are not enough"//&
                    "PW coefficients in the vector file."
                ERROR STOP 1
            endif
            call Sort(ks=(/KX, KY, KZ/), kp=NKVal, vscale=vscale, &!<-- args in
                Dp2s=Dp2s, toldk=DBLE(1e-6), G=Vector, &!<-- args in
                pwcoeffz=CoefC, NV=NV, Orb=Orb,&!<-- args in
                w=Weights) ! ---> args out
            deallocate(CoefC)
            if (lso) then ! SOC
                allocate( CoefC2(NV) )
                read(unitvector2, ioStat=iocplx) CoefC2(1:NV)
                if (iocplx.ne.0) then
                    write(*,*)
                    write(*,*) "ERROR: there are not enough PW coefficients."
                    ERROR STOP 1
                endif
                allocate( Weights2(vscale) )
                call Sort(ks=(/KX, KY, KZ/), kp=NKVal, vscale=vscale, &!<-- args in
                    Dp2s=Dp2s, toldk=DBLE(1e-6), G=Vector, &!<-- args in
                    pwcoeffz=CoefC2, NV=NV, Orb=Orb,&!<-- args in
                    w=Weights2) ! ---> args out
                deallocate(CoefC2)
                ! update Bloch spectral weights with the corresponding norms
                if ( (sonorm(jj) + sonorm2(jj)) < 0.9999 .or. &
                        (sonorm(jj) + sonorm2(jj)) > 1.0001 ) then
                    write(*,*) 'norms do not add up to 1'
                    write(*,*) 'jj=', jj
                    write(*,*) 'sonorm(jj)=', sonorm(jj), &
                        'sonorm2(jj)', sonorm2(jj)
                    ERROR STOP 1
                else ! norm check OK, proceed with weights merge
                    Weights = Weights*REAL(sonorm(jj))
                    Weights2 = Weights2*REAL(sonorm2(jj))
                    Weights = Weights + Weights2
                endif
                deallocate(Weights2)
            endif
        else
            allocate( Coef(NV) )
            read(unitvector, ioStat=iocplx) Coef(1:NV)
            if (iocplx.ne.0) then
                write(*,*)
                write(*,*) "ERROR in reading the vector file in real mode [-r]."
                ERROR STOP 1
            endif
            call Sort(ks=(/KX, KY, KZ/), kp=NKVal, vscale=vscale, &!<-- args in
                Dp2s=Dp2s, toldk=DBLE(1e-6), G=Vector, &!<-- args in
                pwcoeffr=Coef, NV=NV, Orb=Orb,&!<-- args in
                w=Weights) ! ---> args out
            deallocate(Coef)
        endif
        ! write results to a file
        do q=1, vscale
            write(unitout,'(5(f11.6))') NKval(1,q), NKVal(2,q),&
                NKVal(3,q), EIGVAL(jj), Weights(q)
        enddo
        deallocate(Weights)
    enddo ! loop ofver eigenvalues
    deallocate(EIGVAL,Vector, NKVal)
    if (lso) deallocate( EIGVAL2, Vector2, sonorm, sonorm2 ) ! SOC
    888  continue 
enddo ! end of vector file

!! Close I/O files

close(unitvector)
close(unitout)
if (lso) then ! SOC
    close(unitvector2)
    close(unitnorm)
    close(unitnorm2)
endif

!! Display finish message

write(*,'(A)') "\/\/\/ UNFOLDING FINISHED SUCCESSFULLY \/\/\/"
write(*,'(A,I0)') "Number of K points processed: ", nkcount
write(*,'(A,A)') "Data was written to: ", trim(casename)//".f2b"//trim(endout)
write(*,'(A)') "Data format: k_1, k_2, k_3, Eigenvalue(Ry), Bloch weight"

!! Compare number of k point found in the vector file with the number of 
!! point in case.klist file

if (.not.(dir)) then
    write(*,'(A)') trim(casename)//&
        '.klist could not be found for comparison.'
else
    if (nkcount.lt.nkpoints) then
        write(*,'(A,I0,A)') trim(casename)//'.klist file has ', &
            nkpoints-nkcount,' K point(s) more than the vector file.'
    elseif (nkcount.gt.nkpoints)  then
        write(*,'(A,I0,A)') trim(casename)//'.klist file has ', &
            nkcount-nkpoints,' K point(s) less than the vector file.' 
    else
        write(*,'(A)') trim(casename)//'.klist matches the vector file.'
    endif
endif
write(*,'(A)') 'If you have questions/suggestions/bugs, '//&
    'please communicate via the WIEN2k mailing list: '//&
    'http://www.wien2k.at/reg_user/mailing_list'
write(*,'(A)') 'If you found results useful and worth publishing, '//&
    'please consider citing these papers:'
write(*,'(A)') '[1] O. Rubel, A. Bokhanchuk, S. J. Ahmed, and E. Assmann '//&
    '"Unfolding the band structure of disordered solids: from bound '//&
    'states to high-mobility Kane fermions", Phys. Rev. B 90, 115202 (2014).'
write(*,'(A)') '[2] L.-W. Wang, L. Bellaiche, S.-H. Wei, and A. Zunger '//&
    '"Majority representation of alloy electronic states", Phys. Rev. '//&
    'Lett. 80, 4725 (1998).'

STOP ! Main part is successfully completed (exit code 0)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                           Help
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

910 CONTINUE
write(*,'(A)') 'Usage options:'
write(*,'(A)') 'fold2Bloch -h # get help'
write(*,'(A)') 'fold2Bloch -r '//&
    'case.vector[_1] "''P11 P12 P13:P21 P22 P23:P31 P32 P33''" # '//&
    'real calculation (inversion symm.) no SO'
write(*,'(A)') 'fold2Bloch -c '//&
    'case.vector[_1] "''P11 P12 P13:P21 P22 P23:P31 P32 P33''" # '//&
    'complex calc. (no inv. symm.) no SO'
write(*,'(A)') 'fold2Bloch '//&
    'case.vector[_1] "''P11 P12 P13:P21 P22 P23:P31 P32 P33''" '//&
    '# complex calc. implied no SO'
write(*,'(A)') 'fold2Bloch -so case.vectorso[_1] case.vectorsodn[_1] '//&
    'case.normsoup[_1] '//&
    'case.normsodn[_1] "''P11 P12 P13:P21 P22 P23:P31 P32 P33''" # '//&
    'spin-orbit without -sp'
    write(*,'(A)') 'fold2Bloch -so case.vectorsoup[_1] case.vectorsodn[_1] '//&
    'case.normsoup[_1] '//&
    'case.normsodn[_1] "''P11 P12 P13:P21 P22 P23:P31 P32 P33''" # '//&
    'spin-orbit with -sp'
    write(*,'(A)') ''
write(*,'(A)') 'Notes:'
write(*,'(A)') '(1) [P] matrix (internally called Dp2s) is used to '//&
    'transform primitive a_p to supercell a_s lattice '//&
    'vectors (same as in VESTA):'
write(*,'(A)') '       a_s(i) = sum_j a_p(j)*P(j,i)      i,j = 1, 2, 3'
write(*,'(A)') '(2) Use quotations to input the [P] matrix _exactly_ as '//&
    'shown in this help'
write(*,'(A)') '(3) Tutorials can be found at '//&
    'https://github.com/rubel75/fold2Bloch-Wien2k/wiki'
STOP ! legitimate exit

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                           Errors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! vector file
911 CONTINUE
ERROR STOP 'error detected when reading vector file!'

! input options
912 CONTINUE
write(*,'(A)') 'ERROR: Unable to recognize command line options. '//&
    'Possible options are:'
write(*,'(A)') 'fold2Bloch -h # get help'
write(*,'(A)') 'fold2Bloch -r '//&
    'case.vector[_1] "''P11 P12 P13:P21 P22 P23:P31 P32 P33''" # '//&
    'real calculation (inversion symm.) no SO'
write(*,'(A)') 'fold2Bloch -c '//&
    'case.vector[_1] "''P11 P12 P13:P21 P22 P23:P31 P32 P33''" # '//&
    'complex calc. (no inv. symm.) no SO'
write(*,'(A)') 'fold2Bloch '//&
    'case.vector[_1] "''P11 P12 P13:P21 P22 P23:P31 P32 P33''" '//&
    '# complex calc. implied no SO'
write(*,'(A)') 'fold2Bloch -so case.vectorso[_1] case.vectorsodn[_1] '//&
    'case.normsoup[_1] '//&
    'case.normsodn[_1] "''P11 P12 P13:P21 P22 P23:P31 P32 P33''" # '//&
    'spin-orbit without -sp'
    write(*,'(A)') 'fold2Bloch -so case.vectorsoup[_1] case.vectorsodn[_1] '//&
    'case.normsoup[_1] '//&
    'case.normsodn[_1] "''P11 P12 P13:P21 P22 P23:P31 P32 P33''" # '//&
    'spin-orbit with -sp'
    write(*,'(A)') ''
write(*,'(A)') 'Notes:'
write(*,'(A)') '(1) [P] matrix (internally called Dp2s) is used to '//&
    'transform primitive a_p to supercell a_s lattice '//&
    'vectors (same as in VESTA):'
write(*,'(A)') '       a_s(i) = sum_j a_p(j)*P(j,i)      i,j = 1, 2, 3'
write(*,'(A)') '(2) Use quotations to input the [P] matrix _exactly_ as '//&
    'shown in this help'
write(*,'(A)') '(3) Tutorials can be found at '//&
    'https://github.com/rubel75/fold2Bloch-Wien2k/wiki'
ERROR STOP

END PROGRAM fold2Bloch
