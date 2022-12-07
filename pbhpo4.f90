
subroutine pbhpo4
  use input_parameter, &
        only: jobtype, Nbeads, TNstep, atom, atom_multi, hist_max, &
              hist_min, Nfile, atom_num, save_beads, FNameBinary1, lattice
  use calc_parameter
  use calc_histogram1D
  use utility
  implicit none
  integer :: i, j, k
  integer :: UnitAtom
  real(8) :: lat_inv(3,3)
  real(8) :: rij(3)
  real(8), allocatable :: s(:,:,:,:)

  print *, "Nunit is ", Nunit
  allocate(s(3,Natom,Nbeads,TNstep))
  do i = 1, Nunit
  end do

!  select case(jobtype)
!    case(191)
!    call distribution_oho
!  end select
contains

  subroutine distribution_oho
  integer :: Ndist = 6
  real(8) :: roh(3)
  real(8), allocatable :: bond(:,:,:)
  allocate(bond(Nbeads,TNstep,Ndist))
  do k = 1, TNstep
  end do
  end subroutine distribution_oho

end subroutine pbhpo4





