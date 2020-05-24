! Count number of k points in the case.klist file

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
