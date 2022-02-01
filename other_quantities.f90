module mod_other_quantities
  use input_parameter
  use calc_parameter, only: data_beads, data_step
  use calc_histogram1D
  implicit none
  real(8), allocatable :: charge(:,:,:), dipole(:,:,:)
  integer, private :: Uinp, ierr
  integer :: atom1, atom2

contains

  subroutine  other_quantities
!    atom1 = other_atom1
!    atom2 = other_atom2
    atom1 = atom_num(1,1)
    atom2 = atom_num(2,1)

    select case(jobtype)
      case(61:62)
        call charge_analysis
      case(63)
        call dipole_analysis
    end select
  end subroutine other_quantities

  subroutine charge_analysis
    integer :: Istep, i, j, k
!    real(8) :: average(Natom)
    print *, other_path
    allocate(charge(Natom,Nbeads,TNstep))

    open(newunit=Uinp, file=other_path, status='old', iostat=ierr)
      if ( ierr > 0 ) then
        print *, 'Check the path : ', other_path
        stop 'ERROR!!: There is no "charge.dat"'
      end if

      read(Uinp,'()')
      do Istep = 1, Nstart(1)-1
        read(Uinp,'()')
        do j = 1, Nbeads
          read(Uinp,'()')
        end do
      end do

      Istep = 0
      do k = Nstart(1), Nstep(1)
        Istep = Istep + 1
        read(Uinp,'()')
        do j = 1, Nbeads
          read(Uinp,*) charge(:,j,Istep)
        end do
      end do
    close(Uinp)

    if ( jobtype == 62 ) then
    block
      integer :: Ounit
      character(len=128) :: name_out
      if ( save_beads .eqv. .True. ) then
        open(newunit=Ounit, file=FNameBinary1, form='unformatted', access='stream', status='replace')
          do Istep = 1, TNstep
            do j = 1, Nbeads
              write(Ounit) charge(atom1,j,Istep)
            end do
          end do
        close(Ounit)
      end if
      write(*,*) "***** atomic charge of ", atom1, "is saved *****"
      write(*,*) "***** in ", FNameBinary1, " *****"
      data_beads = charge(atom1,:,:)
      name_out = "hist_charge.out"
      call calc_1Dhist(out_hist_ex=name_out)


      open(newunit=Ounit,file="step_charge.out",status='replace')
        do k = 1, TNstep
          if (mod(k,graph_step) == 0 ) write(Ounit,'(I7,F10.5)') k, sum(charge(atom1,:,k))/dble(Nbeads)
        end do
      close(Ounit)
    end block
    end if

!    average(:) = 0.0d0
!    do i = 1, Natom
!      do j = 1, Nbeads
!        average(i) = average(i) + sum(charge(i,j,:))
!      end do
!      average(i) = average(i) / dble(Nbeads)
!    end do
!    do i = 1, Natom
!      print '(I3,F12.7)', i, average(i)
!    end do

    do i = 1, Natom
      print '(I3,F12.7)', i, sum(charge(i,:,:)) / dble(TNstep*Nbeads)
    end do

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
        ! abs_dipole = abs_dipole + dsqrt( dot_product(dipole(:,j,Istep),dipole(:,j,Istep)) )
        data_beads(j,Istep) = dsqrt( dot_product(dipole(:,j,Istep),dipole(:,j,Istep)) )
      end do
    end do
    ! abs_dipole = abs_dipole/dble(TNstep*Nbeads)
    abs_dipole = sum(data_beads)/dble(TNstep*Nbeads)
    print '(a)', "The absolute dipole moment"
    print *, abs_dipole, " D "

block
  integer :: Ounit
  if ( save_beads .eqv. .True. ) then
    open(newunit=Ounit,file=FNameBinary1, form='unformatted', access='stream', status='replace')
      do Istep = 1, TNstep
        do i = 1, Nbeads
          write(Ounit) data_beads(i,Istep)
        end do
      end do
    close(Ounit)
  end if
end block

    return
  end subroutine dipole_analysis


!  subroutine read_charge
!    integer :: i, j, k, Uinp
!    character(len=128) :: line
!
!! --- Rereading input file for charge analysis ---
!    print '(" *****reading input file for charge analysis*****")'
!    open(newunit=Uinp,file=input_file,status='old', err=900)
!      do
!        read(Uinp,'(a)',end=901) line
!        if ( index(trim(line), "# charge analysis") > 0 ) exit
!      end do
!    close(Uinp)
!! --- Rereading input file for charge analysis ---
!  return
!  900 print *, 'ERROR!!: There is no "# job type"'; stop
!  901 print *, 'ERROR!!: There is no "# charge analysis"'; stop
!  end subroutine read_charge

end module mod_other_quantities
! you can change text to line




