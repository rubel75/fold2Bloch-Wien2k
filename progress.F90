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
