module mod_other_quantities
  use input_parameter
  implicit none
  real(8), allocatable :: charge(:,:), dipole(:,:,:)
  integer, private :: Uinp, ierr

contains

  subroutine  other_quantities
    select case(jobtype)
      case(61)
        call charge_analysis
      case(62)
        call dipole_analysis
    end select
  end subroutine other_quantities

  subroutine charge_analysis
    call read_charge
  end subroutine charge_analysis

  subroutine dipole_analysis
    integer :: i,j,k,Istep
    real(8) :: abs_dipole
    allocate(dipole(3,Nbeads,TNstep))

    open(newunit=Uinp, file=other_path,status='old',iostat=ierr)
      if ( ierr > 0 ) then
        print *, 'Check the path : ', other_path
        stop 'ERROR!!: There is no "dipole.dat"'
      end if

      do k = 1, Nstart(1)-1
        read(Uinp, '()')
        do j = 1, Nbeads
          read(Uinp,*) dipole(:,j,k)
        end do
      end do

      Istep = 0
      do k = Nstart(1), Nstep(1)
        Istep = Istep + 1
        read(Uinp, '()')
        do j = 1, Nbeads
          read(Uinp,*) dipole(:,j,Istep)
        end do
      end do
    close(Uinp)

    abs_dipole = 0.0d0
    do Istep = 1, TNstep
      do j = 1, Nbeads
        abs_dipole = abs_dipole + dsqrt( dot_product(dipole(:,j,Istep),dipole(:,j,Istep)) )
      end do
    end do
    abs_dipole = abs_dipole/dble(TNstep*Nbeads)
    print '(a)', "The absolute dipole moment"
    print *, abs_dipole, " D "

    return
  end subroutine dipole_analysis


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




