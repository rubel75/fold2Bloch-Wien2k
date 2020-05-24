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
 integer, intent(in) :: FX, FY, FZ, NV, Orb
 double complex, allocatable :: TGroupC(:,:,:,:)
 double complex, intent(in) :: CoefC(NV)
 real, Intent(out) ::  Weights(FX*FY*FZ)
 DOUBLE PRECISION, allocatable :: Sums(:)
 integer, intent(in) :: Vector(3,NV)
 integer, allocatable :: counter(:,:,:)
 real :: sumtot
 integer :: remainder_x, remainder_y, remainder_z, j, k, l, p, el

 allocate (TGroupC(FX, FY, FZ, NV-Orb))
 allocate (Sums(FX*FY*FZ))
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
    counter(remainder_x+1, remainder_y+1, remainder_z+1) = &
        counter(remainder_x+1, remainder_y+1, remainder_z+1)+1
 	TGroupC(remainder_x+1, remainder_y+1, remainder_z+1, &
 	    counter(remainder_x+1, remainder_y+1, remainder_z+1))=CoefC(j)
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
