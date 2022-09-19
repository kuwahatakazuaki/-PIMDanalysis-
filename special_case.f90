module mod_special_case
  use input_parameter
  use calc_parameter, only: data_beads, data_step
  use calc_histogram1D
  implicit none
!  real(8), allocatable :: charge(:,:,:), dipole(:,:,:), hfcc(:,:,:)
  integer, private :: Uinp, ierr
  integer :: atom1, atom2, atom3, atom4

contains

  subroutine special_case
    atom1 = atom_num(1,1)
    atom2 = atom_num(2,1)
    atom3 = atom_num(3,1)
    atom4 = atom_num(4,1)

    select case(jobtype)
      case(91)
        call out_plane  ! atom2-atom1-atom3 -> atom1-atom4
    end select
  end subroutine special_case

! +++++++++++++++++++++++
! +++ Start out_plane +++
! +++++++++++++++++++++++
  subroutine out_plane
    integer :: Istep, i, j, k
    character(len=:), allocatable :: out_name
    write(out_name,'(a,I0,"-",a,I0)') trim(atom(atom2)), atom2, trim(atom(atom1)), atom1


!    integer :: Istep, i, j, k
!    allocate(charge(Natom,Nbeads,TNstep))
!
!    open(newunit=Uinp, file=other_path, status='old', iostat=ierr)
!      if ( ierr > 0 ) then
!        print *, 'Check the path : ', other_path
!        stop 'ERROR!!: There is no "charge.dat"'
!      end if
!
!      read(Uinp,'()')
!      do i = 1, Nstart(1)-1
!        read(Uinp,'()')
!        do j = 1, Nbeads
!          read(Uinp,'()')
!        end do
!      end do
!
!      Istep = 0
!      do k = Nstart(1), Nstep(1)
!        Istep = Istep + 1
!        read(Uinp,'()')
!        do j = 1, Nbeads
!          read(Uinp,*) charge(:,j,Istep)
!        end do
!      end do
!    close(Uinp)
!
!    if ( jobtype == 62 ) then
!    block
!      integer :: Ounit
!      character(len=128) :: name_out
!      if ( save_beads .eqv. .True. ) then
!        open(newunit=Ounit, file=FNameBinary1, form='unformatted', access='stream', status='replace')
!          do Istep = 1, TNstep
!            do j = 1, Nbeads
!              write(Ounit) charge(atom1,j,Istep)
!            end do
!          end do
!        close(Ounit)
!      end if
!      write(*,*) "***** atomic charge of ", atom1, "is saved *****"
!      write(*,*) "***** in ", FNameBinary1, " *****"
!      data_beads = charge(atom1,:,:)
!      name_out = "hist_charge.out"
!      call calc_1Dhist(out_hist_ex=name_out)
!
!
!      open(newunit=Ounit,file="step_charge.out",status='replace')
!        do k = 1, TNstep
!          if (mod(k,graph_step) == 0 ) write(Ounit,'(I7,F10.5)') k, sum(charge(atom1,:,k))/dble(Nbeads)
!        end do
!      close(Ounit)
!    end block
!    end if
!
!    do i = 1, Natom
!      print '(I3,F12.7)', i, sum(charge(i,:,:)) / dble(TNstep*Nbeads)
!    end do
  end subroutine out_plane
! +++++++++++++++++++++
! +++ End out_plane +++
! +++++++++++++++++++++

end module mod_special_case




