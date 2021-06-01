
subroutine read_coor !(Ifile,step)
  use input_parameter, only : Natom, Nbeads, Nstep, Nfile, FileName, Nstart, atom, FIbinary
  use calc_parameter,  only : r !r(xyz,Natom_Nbead,Nstep)
  implicit none
  integer :: Ifile, step
  integer :: i, j, k, check

  step = 0
  do Ifile = 1, Nfile
  select case(FIbinary)
    case(.False.)
      call read_coor_format(Ifile,step)
    case(.True.)
      atom(:) = ""
      call read_coor_binary(Ifile,step)
  end select
  end do

  print '(a,/)', " *****END reading coordinate*****"

  return
contains

subroutine read_coor_binary(Ifile,step)
  integer, intent(in) :: Ifile
  integer, intent(inout) :: step
  real(8) :: dummy(3)

  ! --- Reading binary file ---
!  open(21, file=trim(FileName(Ifile)), status='old', err=902)
  open(21, file=trim(FileName(Ifile)), form='unformatted', access='stream', status='old',err=902)

  ! --- skip reading coor file ---
    do k = 1, Nstart(Ifile)-1
      do j = 1, Nbeads
        do i = 1, Natom
          read(21,end=911) dummy(:)
        end do
      end do
    end do
  ! --- End skip reading coor file ---

    do k = Nstart(Ifile), Nstep(Ifile)
      step = step+1
      do j = 1, Nbeads
        do i = 1, Natom
          read(21,end=911) r(:,i,j,step)
        end do
      end do
    end do
  close(21)
  ! --- Reading binary file ---

  return
  902 print *, "ERROR!!: There is no coordinate file"; stop
  911 print *, "ERROR!!: Reading line exceed the coor lines"; stop
end subroutine read_coor_binary



subroutine read_coor_format(Ifile,step)
  integer, intent(in) :: Ifile
  integer, intent(inout) :: step

  ! --- Reading formated file ---
  open(21, file=trim(FileName(Ifile)), status='old', err=902)
    read(21,*) check
      if (check /= Natom * Nbeads) then; goto 904; end if
    rewind(21)

  ! --- skip reading coor file ---
    do k = 1, Nstart(Ifile)-1
      read(21,'()',end=911)
      read(21,'()',end=911)
      do j = 1, Nbeads
        do i = 1, Natom
          read(21,'()',end=911)
        end do
      end do
    end do
  ! --- End skip reading coor file ---

    do k = Nstart(Ifile), Nstep(Ifile)
      read(21, '()',end=911)
      read(21, '()',end=911)
      step = step+1
      do j = 1, Nbeads
        do i = 1, Natom
          read(21,*,end=911) atom(i), r(:,i,j,step)
        end do
      end do
    end do
  close(21)
  ! --- Reading formated file ---

  return
  902 print *, "ERROR!!: There is no coordinate file"; stop
  904 print *, "ERROR!!: The number of atom or beads is wrong"; stop
  911 print *, "ERROR!!: Reading line exceed the coor lines"; stop
end subroutine read_coor_format


end subroutine read_coor



