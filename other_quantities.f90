module mod_other_quantities
  use input_parameter
  implicit none
  character(len=:), allocatable :: path_charge

contains

  subroutine charge_analysis
    call read_charge
  end subroutine charge_analysis


  subroutine read_charge
    integer :: i, j, k, Uinp
    character(len=128) :: line
!    character(:), allocatable :: input_file
!
!    print '(" *****START reading parameters*****")'
!    block ! Reading input file
!      integer :: Nargu, leng
!      Nargu = command_argument_count()
!      if ( Nargu == 0) then
!        print '(a)',   "   There is no argument"
!        print '(a,/)', '   Reading from "input.dat"'
!        allocate(character(9) :: input_file)
!        write(input_file,'(a)') "input.dat"
!      else
!        call get_command_argument(1, length=leng)
!          allocate(character(leng) :: input_file)
!          call get_command_argument(1, input_file)
!        print '(a,a,/)', "   Reading from ", '"'//input_file//'"'
!      end if
!    end block ! End Reading input file

! --- Rereading input file for charge analysis ---
    print '(" *****reading input file for charge analysis*****")'
    open(newunit=Uinp,file=input_file,status='old', err=900)
      do
        read(Uinp,'(a)',end=901) line
        if ( index(trim(line), "# charge analysis") > 0 ) exit
      end do
    close(Uinp)
! --- Rereading input file for charge analysis ---
  return
  900 print *, 'ERROR!!: There is no "# job type"'; stop
  901 print *, 'ERROR!!: There is no "# charge analysis"'; stop
  end subroutine read_charge

end module mod_other_quantities
! you can change text to line




