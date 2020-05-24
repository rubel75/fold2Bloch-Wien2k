!!!#################################################################################
!!! NewK
!!!	Finds new K Values for each group, depending on the number of folds in 
!!!	x, y, and z directions.
!!!Input: 	KX, KY, KZ, FoldX, FoldY, FoldZ
!!!Output:	NKVal(array)
!!!#################################################################################

SUBROUTINE NewK (X, Y, Z, FX, FY, FZ, NKVal)
 implicit none
 integer, intent(in) :: FX, FY, FZ
 real, intent(out) :: NKVal(3,FX*FY*FZ)
 real :: field_x, field_y, field_z
 DOUBLE PRECISION :: TX, TY, TZ
 double precision, intent(in) :: X, Y, Z

 integer :: loop, i, j, k, size 

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
