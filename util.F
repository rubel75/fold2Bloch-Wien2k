!!! wien2wannier/util/util.f90
!!!
!!!    Collection of routines for the programs in util/ and woptic/
!!!
!!! Copyright 2010-2012 Philipp Wissgott
!!! Copyright 2013      Elias Assmann

!--------------- Assorted miscellanea                         ---------------
MODULE util
  implicit none

  integer, parameter :: BUFSZ = 256
  !! Default kinds
  integer, parameter :: DP = selected_real_kind(15,300) ! inherited from W90
  private
  public :: string, inverse3x3

  interface string
     module procedure int2str, real2str
  end interface string
contains
  character(len=10) function int2str(n)
    integer, intent(in) :: n
    write(int2str,"(I5)") n
  end function int2str

  character(len=15) function real2str(r)
    real(DP), intent(in) :: r
    write(real2str,"(E16.9)") r
  end function real2str

  subroutine inverse3x3(a, ainv)
    !inverse of th 3x3 matrix A
    implicit none

    real(DP), intent(in)  :: a(3,3)
    real(DP), intent(out) :: ainv(3,3)
    real(DP) :: det

    det = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) &
         +a(1,3)*a(2,1)*a(3,2) - a(3,1)*a(2,2)*a(1,3) &
         -a(1,1)*a(3,2)*a(2,3) - a(2,1)*a(1,2)*a(3,3)

    ainv(1,1) = (  a(2,2)*a(3,3) - a(2,3)*a(3,2) ) / det
    ainv(2,1) = (- a(2,1)*a(3,3) + a(2,3)*a(3,1) ) / det
    ainv(3,1) = (  a(2,1)*a(3,2) - a(2,2)*a(3,1) ) / det
    ainv(1,2) = (- a(1,2)*a(3,3) + a(1,3)*a(3,2) ) / det
    ainv(2,2) = (  a(1,1)*a(3,3) - a(1,3)*a(3,1) ) / det
    ainv(3,2) = (- a(1,1)*a(3,2) + a(1,2)*a(3,1) ) / det
    ainv(1,3) = (  a(1,2)*a(2,3) - a(1,3)*a(2,2) ) / det
    ainv(2,3) = (- a(1,1)*a(2,3) + a(1,3)*a(2,1) ) / det
    ainv(3,3) = (  a(1,1)*a(2,2) - a(1,2)*a(2,1) ) / det
  end subroutine inverse3x3
end MODULE util


!--------------- Helper procedures for command line interface ---------------
module clio
#include "fetcharg.h"
  use iso_fortran_env, only: ERROR_UNIT
  use util,            only: string
  
  implicit none
  private

  integer, parameter :: BUFSZ = 256

  public :: croak, carp, fetcharg_buf
#ifdef HAVE_VARLEN_STR
  public :: fetcharg
#endif

  argstr  :: progname
  logical :: got_progname = .false.
contains
  subroutine get_progname()
    integer s

    call fetcharg(0, progname, status=s)

    if (s /= 0) progname='WIEN2WANNIER'

    got_progname = .true.
  end subroutine get_progname

  subroutine croak(message, status)
    character(len=*), intent(in) :: message
    integer, intent(in), optional :: status

    integer            :: s=1

    if (present(status)) s=status

    if (.not. got_progname) call get_progname()
    
    write(ERROR_UNIT, '(A, ": ", A)') progname, message
    call exit(s)
  end subroutine croak

  subroutine carp(message)
    character(len=*), intent(in) :: message

    if (.not. got_progname) call get_progname()

    write(ERROR_UNIT, '(A, ": ", A)') progname, message
  end subroutine carp

  !! This elegant variable-length version of `fetcharg´ does not work
  !! on pre-4.8 gcc
#ifdef HAVE_VARLEN_STR
  subroutine fetcharg(i, str, message, status)
    integer,          intent(in)               :: i
    character(len=:), intent(out), allocatable :: str
    character(len=*), intent(in),  optional    :: message
    integer,          intent(out), optional    :: status

    integer :: s, l

    if (allocated(str)) deallocate(str)

    call get_command_argument(i, length=l, status=s)
    if (present(status)) then
       status = s
    elseif (s /= 0) then
       if (present(message)) then
          call croak(message)
       else
          call croak("FETCHARG: failed to get command argument #" // &
               &     trim(string(i)) // " length: " // trim(string(s)))
       end if
    end if

    allocate(character(len=l) :: str)
    !! zero-length arguments seem to need special treatment
    !! (hooray Fortran!)
    if (l==0) return

    call get_command_argument(i, value=str, status=s)
    if (present(status)) then
       status = s
    elseif (s /= 0) then
       if (present(message)) then
          call croak(message)
       else
          call croak("FETCHARG: failed to get command argument #" // &
               &     trim(string(i)) // ": " // trim(string(s)))
       end if
    end if
  end subroutine fetcharg
#endif

  subroutine fetcharg_buf(i, buf, message, status)
    integer,          intent(in)               :: i
    character(len=*), intent(out)              :: buf
    character(len=*), intent(in), optional     :: message
    integer,          intent(out), optional    :: status

    integer :: s

    !! here, a zero-length argument is okay
    !! (probably because len(buf) /= 0)
    call get_command_argument(i, value=buf, status=s)
    if (present(status)) then
       status = s
    elseif (s /= 0) then
       if (present(message)) then
          call croak(message)
       elseif (s < 0) then
          call croak("FETCHARG_BUF: buffer too small for command argument #" &
               &     // trim(string(i)))
       else
          call croak("FETCHARG_BUF: failed to get command argument #" // &
               &     trim(string(i)) // ": " // trim(string(s)))
       end if
    end if
  end subroutine fetcharg_buf
end module clio


!----------------- ‘struct’ type and associated procedures  -----------------
module structmod
  implicit none

  private
  public :: struct, struct_read

  type intptr
     integer, allocatable :: p(:)
  end type intptr

  type struct
     !! Structure to represent a ‘struct’ file
     !!
     !! The conventions here follow mostly those of ‘structeditor’

     character(len=80) :: title
     character(len= 4) :: lattic       ! lattice type
     integer           :: nat, nneq    ! number of noneq. and total atoms
     character(len= 4) :: mode         ! ‘RELA’ or ‘NREL’
     real(8)           :: a(3)         ! lattice constants
     real(8)           :: alpha(3)     ! angles
     real(8)           :: brlat(3,3)   ! Bravais lattice (row is a vector)
     real(8)           :: lat2car(3,3) ! lattice to cartesian transformation
                                       ! (different from brlat in some cases)
     logical           :: ortho
     real(8)           :: vol          ! u.c. volume

     ! positions (3 × nneq); local rotation matrices (3 × 3 × nat)
     real(8),           allocatable :: pos(:,:), locrot(:,:,:)
     ! neq2at(ineq) is the atom numer corresponding to ineq
     integer,           allocatable :: mult(:), isplit(:), npt(:), neq2at(:)
     ! at2neq(iat)%p is an array of the ‘ineq’s corresponding to iat
     type(intptr),      allocatable :: at2neq(:)
     character(len=10), allocatable :: aname(:)
     real(8),           allocatable :: r0(:), rmt(:), Z(:)
  end type struct
contains
  subroutine struct_read(lun, stru)
    use clio,  only: croak
    use util,  only: inverse3x3

    integer,      intent(in)  :: lun
    type(struct), intent(out) :: stru
    integer, parameter :: DP = selected_real_kind(15,300) ! inherited from W90
    real(DP), parameter :: PI         = 3.1415926535897932d0
    real(DP), parameter :: ORTHO_TEST = 1.d-12
    real(DP), parameter :: SQ3        = sqrt(3d0)
    integer :: iat, ineq, N, i
    ! “Bravais lattice” and “lattice->cartesian” transformation
    ! matrices [in some cases, the two are different], both for direct
    ! and reciprocal space
    real(8) :: br1_dir(3,3), br1_rec(3,3), br2_dir(3,3), br2_rec(3,3)
    ! abbreviations
    real(8) :: pia(3), alpha(3), cosab, cosac, cosbc, sinab, sinbc
    real(8) :: rvfac, wurzel

    read(lun, '(A)')          stru%title
    read(lun, '(A, 23X, I3)') stru%lattic, stru%nat
    read(lun, '(13X, A)')     stru%mode
    read(lun, '(6F10.6)')     stru%a, stru%alpha

    where (abs(stru%alpha) < 1e-5) stru%alpha = 90

    N = stru%nat
    allocate(stru%mult(N), stru%isplit(N), stru%aname (N), &
         &   stru%npt (N), stru%at2neq(N), stru%Z     (N), &
         &   stru%r0  (N), stru%rmt   (N), stru%locrot(3,3,N))

    countneq: do iat = 1, stru%nat
       read(lun,*) ! first position
       read(lun, '(15X, I2, 17X, I2)') stru%mult(iat), stru%isplit(iat)
       ! rest of positions
       do ineq = 2, stru%mult(iat)
          read(lun,*)
       end do
       read(lun, '(A, 5X, I5, 5X, F10.8, 5X, F10.5, 5X, F5.2)') &
            stru%aname(iat), stru%npt(iat), stru%r0(iat), stru%rmt(iat), &
            stru%Z(iat)
       read(lun,'(20X, 3F10.7)') (stru%locrot(1, i, iat), i=1,3)
       read(lun,'(20X, 3F10.7)') (stru%locrot(2, i, iat), i=1,3)
       read(lun,'(20X, 3F10.7)') (stru%locrot(3, i, iat), i=1,3)
    end do countneq

    stru%nneq = sum(stru%mult)
    N = stru%nneq
    allocate(stru%pos(3, N), stru%neq2at(N))

    N=0
    do iat = 1, stru%nat
       stru%neq2at(N+1:N+stru%mult(iat)) = iat

       allocate(stru%at2neq(iat)%p(stru%mult(iat)))
       stru%at2neq(iat)%p = iat

       N = N+stru%mult(iat)
    end do

    rewind(lun)
    read(lun,*) ! title
    read(lun,*) ! lattic
    read(lun,*) ! mode
    read(lun,*) ! uc
    N=0
    readpos: do iat = 1, stru%nat
       N=N+1
       read(lun, '(9X, 3(3X, F10.8))') stru%pos(:, N)
       read(lun,*) ! mult
       do ineq = 2, stru%mult(iat)
          N=N+1
          read(lun, '(9X, 3(3X, F10.8))') stru%pos(:, N)
       end do
       read(lun,*) ! aname
       read(lun,*) ! locrot
       read(lun,*) !
       read(lun,*) !
    end do readpos

    !! Now compute lattice vectors.  This is copied from
    !! SRC_structeditor/module.f (subroutine latgen_struct)
    pia = 2*PI/stru%a
    alpha = stru%alpha*PI/180

    cosab = cos(alpha(3)); sinab = sin(alpha(3))
    cosac = cos(alpha(2))
    cosbc = cos(alpha(1)); sinbc = sin(alpha(1))

    br1_rec = 0; br1_dir = 0; br2_rec = 0; br2_dir = 0

    lattic: select case (stru%lattic(1:1))
    case ('H')
       br1_rec(1,1) = pia(1)*2/SQ3
       br1_rec(1,2) = pia(1)  /SQ3
       br1_rec(2,2) = pia(2)
       br1_rec(3,3) = pia(3)

       br2_rec = br1_rec

       rvfac = 2/SQ3
       stru%ortho = .false.

    case ('S', 'P')             ! what's ‘S’??
       wurzel = sinbc**2 - cosac**2 - cosab**2 + 2*cosbc*cosac*cosab

       br1_rec(1,1) = pia(1) * sinbc/wurzel
       br1_rec(1,2) = pia(2) * (-cosab + cosbc*cosac)/(sinbc*wurzel)
       br1_rec(1,3) = pia(3) * (-cosac + cosbc*cosab)/(sinbc*wurzel)
       br1_rec(2,2) = pia(2) / sinbc
       br1_rec(2,3) =-pia(3) * cosbc/sinbc
       br1_rec(3,3) = pia(3)

       br2_rec = br1_rec

       rvfac = 1/wurzel

       stru%ortho = all(abs(alpha - PI/2) <= ortho_test)

    case ('F')
       br1_rec(1,1) = pia(1)
       br1_rec(2,2) = pia(2)
       br1_rec(3,3) = pia(3)

       br2_rec(1,:) = pia(1) * (/-1, 1, 1 /)
       br2_rec(2,:) = pia(2) * (/ 1,-1, 1 /)
       br2_rec(3,:) = pia(3) * (/ 1, 1,-1 /)

       rvfac = 4
       stru%ortho = .true.
    case ('B')
       br1_rec(1,1) = pia(1)
       br1_rec(2,2) = pia(2)
       br1_rec(3,3) = pia(3)

       br2_rec(1,:) = pia(1) * (/ 0, 1, 1 /)
       br2_rec(2,:) = pia(2) * (/ 1, 0, 1 /)
       br2_rec(3,:) = pia(3) * (/ 1, 1, 0 /)

       rvfac = 2
       stru%ortho = .true.

    case ('R')
       br1_rec(1, :) = pia(1) * (/  1, 1, -2 /)/SQ3
       br1_rec(2, :) = pia(2) * (/ -1, 1,  0 /)
       br1_rec(3, :) = pia(3)

       br2_rec = br1_rec

       rvfac = 6/SQ3
       stru%ortho = .false.

    case ('C')
       ! “defaults”, to be changed in nonorthogonal XZ-case
       br1_rec(1,1) = pia(1)
       br1_rec(2,2) = pia(2)
       br1_rec(3,3) = pia(3)

       rvfac = 2
       stru%ortho = .true.

       C: select case (stru%lattic(2:3))
       case ('XY')
          br2_rec(1,:) = pia(1) * (/ 1, 1, 0 /)
          br2_rec(2,:) = pia(2) * (/-1, 1, 0 /)
          br2_rec(3,:) = pia(3) * (/ 0, 0, 1 /)

       case ('YZ')
          br2_rec(1,:) = pia(1) * (/ 1, 0, 0 /)
          br2_rec(2,:) = pia(2) * (/ 0, 1, 1 /)
          br2_rec(3,:) = pia(3) * (/ 0,-1, 1 /)

       case ('XZ')
          ortho: if (abs(alpha(3) - PI/2) <= 0.0001) then
             br2_rec(1,:) = pia(1) * (/ 1, 0, 1 /)
             br2_rec(2,:) = pia(2) * (/ 0, 1, 0 /)
             br2_rec(3,:) = pia(3) * (/-1, 0, 1 /)
          else                  ! CXZ monoclinic case
             br1_rec(1,1) = pia(1) / sinab
             br1_rec(1,2) =-pia(2) / sinab * cosab

             br2_rec(1,1) = pia(1) / sinab
             br2_rec(1,2) =-pia(2) / sinab * cosab
             br2_rec(1,3) = pia(1) / sinab
             br2_rec(2,:) = pia(2) * (/ 0, 1, 0 /)
             br2_rec(3,:) = pia(3) * (/-1, 0, 1 /)

             rvfac = 2/sinab
             stru%ortho = .false.
          end if ortho
       end select C

    case default
       rvfac = 0                ! silence warning
       call croak('unknown lattice type `' // trim(stru%lattic) // "'")
    end select lattic

    call inverse3x3(br1_rec, br1_dir); br1_dir = br1_dir*2*PI
    call inverse3x3(br2_rec, br2_dir); br2_dir = br2_dir*2*PI

    stru%vol = product(stru%a) / rvfac

    stru%lat2car = br1_dir
    stru%brlat   = br2_dir
  end subroutine struct_read
end module structmod
