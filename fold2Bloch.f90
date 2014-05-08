!!!##################################################################################
!!! fold2Bloch.f90		
!!!TEST 
!!!    Unfold the data in .vector[_1] file.  Based on
!!!    join_vectorfiles.
!!!
!!!    Usage: fold2Bloch [-up/-dn] [-c] <case.vector[_1]> x:y:z (x,y, and z are number of folds)
!!!
!!! Copyright 2013 Elias Assmann, Anton Bokhanchuk
!!!###################################################################################

!!!#################################################################################
!!!Progress
!!!	Write a progress bar to the screen indicating percentage of work complete.
!!!Input: 	nkcount, nkpoints
!!!Output:	none
!!!#################################################################################

SUBROUTINE progress(symbol)
 implicit none
 character, intent(inout):: symbol
 character(len=17)::bar="? processing file"
 if (symbol.eq.'\') then
    symbol='|'
 elseif (symbol.eq.'|') then
    symbol='/'
 elseif (symbol.eq.'/') then
    symbol='-'
 elseif (symbol.eq.'-') then
    symbol='\'
 endif
 write(unit=bar(1:1),fmt="(a1)") symbol
 bar(2:1) = symbol
 ! print the progress bar.
 write(unit=6,fmt="(a1,a23,$)") char(13), bar
 return 
END SUBROUTINE progress


MODULE Unfold
Contains

integer FUNCTION line_count(ident)

 implicit none
 integer, intent(in) :: ident
 character(20) temp
 integer iostat
 logical ioEndStat

 line_count=0
 ioEndStat=.false.

 do while (.not. ioEndStat )
   read(ident,"(A20)", iostat=ioStat ) temp
   if( temp .ne. 'END') then 
      line_count = line_count + 1
   else 
      ioEndStat = .TRUE. 
   endif
 end do
 rewind(ident)
END FUNCTION line_count        

!!!#################################################################################
!!! NewK
!!!	Finds new K Values for each group, depending on the number of folds in 
!!!	x, y, and z directions.
!!!Input: 	KX, KY, KZ, FoldX, FoldY, FoldZ
!!!Output:	NKVal(array)
!!!#################################################################################

SUBROUTINE NewK (X, Y, Z, FX, FY, FZ, NKVal)
 implicit none
 real, allocatable :: NKVal (:,:)
 real :: field_x, field_y, field_z
 DOUBLE PRECISION :: TX, TY, TZ
 double precision, intent(in) :: X, Y, Z
 integer, intent(in) :: FX, FY, FZ
 integer :: loop, i, j, k, size 
 size=FX*FY*FZ
 allocate (NKVal(3,size))

 !Brillouin zone range
 field_x=0.5*FX
 field_y=0.5*FY
 field_z=0.5*FZ
 loop=1 

 do i=0, (FX-1)
    if ((X+i).gt.(field_x)) then
       TX=X-(field_x*2)+i
    else 
	TX=X+i
    endif
    do j=0, (FY-1)
       if ((Y+j).gt.(field_y)) then
          TY=Y-(field_y*2)+j
       else 
	   TY=Y+j
	endif
       do k=0,(FZ-1)
          if ((Z+k).gt.(field_z)) then
             TZ=Z-(field_z*2)+k
          else 
             TZ=Z+k
	   endif
	   NKVal(1,loop)=TX
           NKVal(2,loop)=TY
           NKVal(3,loop)=TZ
           loop=loop+1
       enddo
    enddo
 enddo   
END SUBROUTINE NewK

!!!#################################################################################
!!! Sort
!!!     Sorts energy coefficients into appropriate groups according to relative
!!!     location to Gamma point. Then calculates the weight of each group.
!!!Input:       FoldX, FoldY, FoldZ, Vectors, Coef, NV, Orb
!!!Output:      Weights(array)
!!!#################################################################################
SUBROUTINE Sort(FX, FY, FZ, Vector, Coef, NV, Orb, Weights)
 implicit none
 double precision, allocatable :: TGroup(:,:,:,:)
 double precision, allocatable, intent(in) :: Coef(:)
 real, allocatable, Intent(inout) ::  Weights(:)
 DOUBLE PRECISION, allocatable :: Sums(:)
 integer, allocatable, intent(in) :: Vector(:,:)
 integer, allocatable :: counter(:,:,:)
 real :: sumtot
 integer :: remainder_x, remainder_y, remainder_z, j, k, l, p, el
 integer, intent(in) :: FX, FY, FZ, NV, Orb
 allocate (TGroup(FX, FY, FZ, NV-Orb))
 allocate (Sums(FX*FY*FZ))
 allocate (Weights(FX*FY*FZ))
 allocate (counter(FX, FY, FZ))

 !Initiates the counter and TGroupC elements at 0
  do j=1,FX
   do k=1,FY
     do l=1,FZ
       counter(j,k,l)=0
       do p=1,NV-Orb
         TGroup(FX,FY,FZ,p)=0
       enddo
     enddo
   enddo
 enddo
 
 !Sorts the energy coeeficients
 do j=1, (NV-Orb)
        remainder_x=MODULO(Vector(1,j),FX)
        remainder_y=MODULO(Vector(2,j),FY)
        remainder_z=MODULO(Vector(3,j),FZ)
        counter(remainder_x+1, remainder_y+1, remainder_z+1)=counter(remainder_x+1, remainder_y+1, remainder_z+1)+1
 	 TGroup(remainder_x+1, remainder_y+1, remainder_z+1, counter(remainder_x+1, remainder_y+1, remainder_z+1))=Coef(j)
 enddo

 !Sums the squares  of all coefficients per group
 el=1
 do j=1, FX
        do k=1, FY
                do l=1, FZ
			if (counter(j, k, l).gt.0) then
			   do p=1, counter(j, k, l)
			     TGroup(j, k, l,p)=TGroup(j, k, l,p)*(TGroup(j, k, l,p))
			   enddo
			   Sums(el)=SUM(TGroup(j, k, l,1:counter(j, k, l)))
			   el=el+1
			else 
			   Sums(el)=0.0
			   el=el+1
		        endif
                enddo
        enddo
 enddo
 sumtot=SUM(Sums(:))
 do j=1, (FX*FY*FZ)
	Weights(j)=Sums(j)/sumtot
 enddo
 deallocate (TGroup, Sums)
END SUBROUTINE Sort

!!!#################################################################################
!!! SortC
!!!     Sorts energy coefficients (compex number form) into appropriate groups 
!!!	according to relative location to Gamma point. Then calculates the weight 
!!!	of each group.
!!!Input:       FoldX, FoldY, FoldZ, Vectors, CoefC, NV, Orb
!!!Output:      Weights(array)
!!!#################################################################################
SUBROUTINE SortC(FX, FY, FZ, Vector, CoefC, NV, Orb, Weights)
 implicit none
 double complex, allocatable :: TGroupC(:,:,:,:)
 double complex, allocatable, intent(in) :: CoefC(:)
 real, allocatable, Intent(inout) ::  Weights(:)
 DOUBLE PRECISION, allocatable :: Sums(:)
 integer, allocatable, intent(in) :: Vector(:,:)
 integer, allocatable :: counter(:,:,:)
 real :: sumtot
 integer :: remainder_x, remainder_y, remainder_z, j, k, l, p, el
 integer, intent(in) :: FX, FY, FZ, NV, Orb
 allocate (TGroupC(FX, FY, FZ, NV-Orb))
 allocate (Sums(FX*FY*FZ))
 allocate (Weights(FX*FY*FZ))
 allocate (counter(FX,FY,FZ))

 !Initiates the counter and TGroupC elements at 0
 do j=1,FX
   do k=1,FY
     do l=1,FZ
       counter(j,k,l)=0
       do p=1,NV-Orb
         TGroupC(j,k,l,p)=0.0
       enddo
     enddo
   enddo
 enddo

 !Sorts the energy coeeficients
 do j=1, (NV-Orb)
        remainder_x=MODULO(Vector(1,j), FX)
        remainder_y=MODULO(Vector(2,j), FY)
        remainder_z=MODULO(Vector(3,j), FZ)
        counter(remainder_x+1, remainder_y+1, remainder_z+1)=counter(remainder_x+1, remainder_y+1, remainder_z+1)+1
 	TGroupC(remainder_x+1, remainder_y+1, remainder_z+1, counter(remainder_x+1, remainder_y+1, remainder_z+1))=CoefC(j)
 enddo

 !Sums the squares  of all coefficients per group
 el=1
 do j=1, FX
        do k=1, FY
                do l=1, FZ
			if (counter(j, k, l).gt.0) then
			   do p=1, counter(j, k, l)
			     TGroupC(j, k, l,p)=TGroupC(j, k, l,p)*conjg(TGroupC(j, k, l,p))
			   enddo
			   Sums(el)=SUM(TGroupC(j, k, l,1:counter(j, k, l)))
			   el=el+1
			else 
			   Sums(el)=0.0
			   el=el+1
		        endif
                enddo
        enddo
 enddo
 sumtot=SUM(Sums)
 do j=1, (FX*FY*FZ)
	Weights(j)=Sums(j)/sumtot
 enddo 
 deallocate (TGroupC, Sums)
END SUBROUTINE SortC
END MODULE Unfold

!!!##################################################################################
!!! 					MAIN
!!!###################################################################################
PROGRAM Unfolder
  use structmod, only: struct, struct_read
  use Unfold

  implicit none

 character(70) seedname, vectorname
 character(20) :: folds
 integer num_args,argcount !command line input argument counter
 integer jatom,i,j,k,jl,jj,l,q
 integer lmax,lomax,nloat, field
 integer unitklist,unitstruct,unitvector, unittargetvector
 character(70)startmessage
 character(30), allocatable :: args(:)
 logical :: usecomplex, dir
 INTEGER            NE, NV, ORB, FoldX, FoldY, FoldZ, ios, iocplx, nkcount, nkpoints
 DOUBLE PRECISION   KX, KY, KZ, WEIGHT
 CHARACTER(10)      KNAME
 character :: bar
 real, allocatable :: Weights(:)
 real, allocatable :: NKVal(:,:)
 DOUBLE PRECISION, allocatable ::  EIGVAL(:)
 DOUBLE PRECISION, allocatable ::   E(:,:),ELO(:,:,:)
 INTEGER, allocatable :: Vector(:,:)
 DOUBLE PRECISION, allocatable ::  Coef(:)
 DOUBLE COMPLEX, allocatable ::  CoefC(:)
 type(struct) stru

  !default fileending: non spin-polarized
  startmessage = "     /\/\/\ UNFOLDING VECTOR FILE /\/\/\"
  unitvector = 11
  unittargetvector = 12
  unitklist = 13
  unitstruct = 14
  !parameters
  LMAX = 13
  LOMAX = 3
  nloat = 3
 
  write(*,*) '		***********************'
  write(*,*) '		** Fold2Bloch V 1.04 **'
  write(*,*) '		** Build May 6, 2014 **'
  write(*,*) '		***********************'

  !command line argument read-in
  num_args=command_argument_count()
  allocate (args(num_args))
  argcount=1
  usecomplex = .true.
  do j=1,num_args
     call get_command_argument(j,args(j))
     if (num_args.eq.3) then
        if (j.eq.1) then
           if (args(j).eq.'-r') then    
              usecomplex = .false.
              write(*,*) '     Regular (non-complex) calculation indicated'
           elseif (args(j).eq.'-c') then
              usecomplex = .true.
              write(*,*) '     Complex calculation indicated'
           else
              write(*,*) '     Unknown option. See below or type: "fold2Bloch -h" or "fold2Bloch --help" for more information.'
              write(*,*) '     Usage: fold2Bloch [-r/-c] case.vector[_1] x:y:z (folds)'
              stop    
           endif
	elseif (j.eq.2) then
           read(args(j),*) vectorname
       else
           read(args(j),*) folds
       endif
     elseif (num_args.eq.2) then
	if ((j.eq.1).and.(args(j).eq.'-r'.or.args(j).eq.'-c')) then
           write(*,*) '     Unknown number of folds. See below or type: "fold2Bloch -h" or "fold2Bloch --help" for more information.'
           write(*,*) '     Usage: fold2Bloch [-r/-c] case.vector[_1] x:y:z (folds)'
           stop
       elseif (j.eq.1) then
          write(*,*) '     Complex calculation assumed (default)' 
          read(args(j),*) vectorname
       else 
           read(args(j),*) folds
       endif
     elseif ((args(j).eq.'-h').or.(args(j).eq.'--help')) then
       write(*,*) '     Usage: fold2Bloch [-r] case.vector[_1] x:y:z (folds)'
       write(*,*) '     -r is an option indicating no oomplex number calculations use'
       write(*,*) '     case.vector[_1] is a vector file name (ex. atom.vector, atom.vector_1)'
       write(*,*) '     x:y:z integers, greater than 0, that represent a multiplicity in the corresponding directions used when constructing the supercell.'
       stop
     else  
       write(*,*) '     Unknown option. See below or type: "fold2Bloch -h" or "fold2Bloch --help" for more information.'
       write(*,*) '     Usage: fold2Bloch [-r/-c] case.vector[_1] x:y:z (folds)'
       stop 
     endif 
  enddo
 if (folds.eq.'') then
   write(*,*) '     Unknown number of folds. See below or type: "fold2Bloch -h" or "fold2Bloch --help" for more information.'
   write(*,*) '     Usage: fold2Bloch [-r/-c] case.vector[_1] x:y:z (folds)'
   stop
 endif
 inquire(file=vectorname, exist=dir)
 if (not(dir)) then
  write(*,50) trim(vectorname)//' CASE FILE DOES NOT EXIST OR NO FILENAME ENTERED'
  50 format(47A)
  write(*,*) '     Usage: fold2Bloch [-r/-c] case.vector[_1] x:y:z (folds)'
  stop
 endif

 i=INDEX(vectorname, '.')
 read (vectorname(1:i-1), *) seedname
 
 !Check if the number of folds argument is in correct format
 i=0
 i=INDEX(folds, ':')
 if (i.eq.0) then
   write(*,*) '     Unknown number of folds. See below or type: "fold2Bloch -h" or "fold2Bloch --help" for more information.'
   write(*,*) '     Usage: fold2Bloch [-r/-c] case.vector[_1] x:y:z (folds)'
   stop
 endif
 read (folds(1:i-1), *, iostat=ios) FoldX
 if ((ios.ne.0).or.(FoldX.le.0)) then
   write (*,*) '     Number of folds has to be a positive integer greater than 0'
   stop
 endif
 folds=folds(i+1:)
 i=0
 i=INDEX(folds, ':')
 if (i.eq.0) then
   write(*,*) '     Unknown number of folds. See below or type: "fold2Bloch -h" or "fold2Bloch --help" for more information.'
   write(*,*) '     Usage: fold2Bloch [-r/-c] case.vector[_1] x:y:z (folds)'
   stop
 endif
 read (folds(1:i-1),*, iostat=ios) FoldY
 if ((ios.ne.0).or.(FoldY.le.0)) then
   write (*,*) '     Number of folds has to be a positive integer greater than 0'
   stop
 endif
 read(folds(i+1:),*, iostat=ios) FoldZ
 if ((ios.ne.0).or.(FoldZ.le.0)) then
   write (*,*) '     Number of folds has to be a positive integer greater than 0'
   stop
 endif
 field=FoldX*FoldY*FoldZ

 !Check if .klist and .struct files exist
 inquire(file=trim(seedname)//'.struct', exist=dir)
 if (not(dir)) then
  write(*,*) "     ", trim(seedname)//'.struct file does not exist'
  stop
 else
  open(unit=unitstruct,file=trim(seedname)//'.struct',status='old')
  call struct_read(unitstruct, stru)
  close(unitstruct)
 endif

 write(*,35)"     FILE TO PROCESS: ", vectorname
 write(*,*)startmessage
 35 format(a23, a70)

 inquire(file=trim(seedname)//'.klist', exist=dir)
 if (not(dir)) then
   write(*,*) trim("     WARNING: ")//trim(seedname)//'.klist could not be found. Progress can not be calculated!'
   write(*,*) "      Continuing to process file..."
 else
   open(unit=unitklist,file=trim(seedname)//'.klist',status='old')
   nkpoints = line_count(unitklist) !Determines the number of K points
   close(unitklist)
 endif
 
 allocate( E(LMAX,stru%nat) )
 allocate( ELO(0:LOMAX,nloat,stru%nat) )
 
 open(unit=unittargetvector,file=trim(seedname)//".unfolded", &
      & status='unknown',form='formatted')   
 open(unit=unitvector, file=vectorname, status='old',form='unformatted')
 do jatom= 1,stru%nat
    read(unitvector) (E(jl,jatom),jl=1,LMAX)
    read(unitvector) ((ELO(jl,k,jatom),jl=0,LOMAX),k=1,nloat)
 enddo
 nkcount=0
 ios=0
 bar='\'
 do while (ios.eq.0)
    read(unitvector,end=888,err=888, iostat=ios) KX, KY, KZ, KNAME, NV, NE, WEIGHT
    nkcount=nkcount+1
    call progress(bar)
    !write(unittargetvector,*) '##############################################'
    !write(unittargetvector,*) KX, KY, KZ
    !write(unittargetvector,*) '##############################################'
    call NewK(KX, KY, KZ, FoldX, FoldY, FoldZ, NKVal)
    allocate( Vector(3,NV) )
    read(unitvector) (Vector(1,I),Vector(2,I),Vector(3,I),I=1,NV)
    !write(unittargetvector,50) (Vector(1,I),Vector(2,I),Vector(3,I),I=1,NV)
    !50 format(i4, i4, i4)
    !write(unittargetvector,*) NV
    do l=2, NV !Determine the number of local orbitals
       if ((Vector(1,l).eq.Vector(1,1)).and.(Vector(2,l).eq.Vector(2,1)).and.(Vector(3,l).eq.Vector(3,1))) then 
          ORB = NV-l+1
          exit
       endif
    enddo
    allocate(EIGVAL(NE))
    do jj = 1, NE
      allocate(Coef(NV),CoefC(NV))
      read(unitvector) I, EIGVAL(I)
      !write(unittargetvector, *) EIGVAL(jj)
      !write(unittargetvector,*) '##############################################'
      if (usecomplex) then
          read(unitvector, ioStat=iocplx) CoefC(1:NV)
          if (iocplx.ne.0) then
            write(*,*)
            write(*,*) "Ooops,  you either forgot [-r] (no complex calculation) switch, or there are not enough coefficients."
            stop
          endif
	   call SortC(FoldX, FoldY, FoldZ, Vector, CoefC, NV, Orb, Weights)
      else
          read(unitvector) Coef(1:NV)
          call Sort(FoldX, FoldY, FoldZ, Vector, Coef, NV, Orb, Weights)
      endif 
      do q=1, field !Writes results to a file
      write(unittargetvector,100) (NKval(1,q)/FoldX), (NKVal(2,q)/FoldY) ,(NKVal(3,q)/FoldZ), EIGVAL(jj), Weights(q)
	100 format(f11.6, f11.6, f11.6, f11.6, f11.6)
      enddo
     deallocate(Weights, CoefC, Coef)
    enddo
    deallocate(EIGVAL,Vector, NKVal)
 888  continue 
 enddo

 write(*,*)
 write(*,*) "     \/\/\/ UNFOLDING FINISHED SUCCESSFULLY \/\/\/"
 write(*,*) "     Number of K points processed:", nkcount
 write(*,*) "     Data was written to: ", trim(seedname)//".unfolded"
 write(*,*) "     Data format: KX, KY, KZ, Eigenvalue(Ry), Weight"

 if (not(dir)) then
   write(*,*) "     ", trim(seedname)//'.klist could not be found for comparison.'
 else
   if (nkcount.lt.nkpoints) then
     write(*,*) "     ", trim(seedname)//'.klist file has ',nkpoints-nkcount,' K point(s) more than the vector file.'
   elseif (nkcount.gt.nkpoints)  then
     write(*,*) "     ", trim(seedname)//'.klist file has ',nkcount-nkpoints,' K point(s) less than the vector file.' 
   else
     write(*,*) "     ", trim(seedname)//'.klist matches the vector file.'
   endif
 endif
 close(unitvector)
 close(unittargetvector)
END PROGRAM Unfolder
