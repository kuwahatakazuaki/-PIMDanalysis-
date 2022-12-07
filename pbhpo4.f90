
subroutine pbhpo4
  use input_parameter, &
        only: jobtype, Nbeads, TNstep, atom, atom_multi, hist_max, &
              hist_min, Nfile, atom_num, save_beads, FNameBinary1, lattice
  use calc_parameter
  use calc_histogram1D
  use utility
  implicit none
  integer :: i, j, k, atom1, atom2
  real(8) :: lat_inv(3,3)
  real(8) :: rij(3), sij(3), dij
  real(8), allocatable :: s(:,:,:,:), unit_bead(:,:,:), unit_step(:,:)
  real(8) :: data_max, data_min, data_ave, data_dev, data_err

  call get_inv_mat(lattice,lat_inv,3)

  allocate(s(3,Natom,Nbeads,TNstep))
  allocate(unit_bead(Nunit*2,Nbeads,TNstep))
  allocate(unit_step(Nunit*2,TNstep))

  do k = 1, TNstep
    do j = 1, Nbeads
      do i = 1, Natom
        s(:,i,j,k) = matmul(r(:,i,j,k), lat_inv(:,:))
      end do
    end do
  end do

  do i = 0,2
    call calc_unit_bond(19+i,22+i,i+1)
  end do
  do i = 0,1
    call calc_unit_bond(26+i,28+i,i+4)
  end do
  call calc_unit_bond(25,30,6)

  do k = 1, TNstep
    do i = 1,Nunit*2
      unit_step(i,k) = sum(unit_bead(i,:,k))/dble(Nbeads)
    end do
  end do


  do k = 1, TNstep
    if (mod(k,graph_step)==0) then
      print 9998, k, unit_step(:,k)
    end if
  end do


!  select case(jobtype)
!    case(191)
!    call distribution_oho
!  end select
9999  format(3F13.7)
9998  format(i0,6F11.6)
contains

  subroutine distribution_oho
  integer :: Ndist = 6
  real(8) :: roh(3)
  real(8), allocatable :: bond(:,:,:)
  allocate(bond(Nbeads,TNstep,Ndist))
  do k = 1, TNstep
  end do
  end subroutine distribution_oho

  subroutine calc_1Dhistogram_unit
    real(8) :: Dhist
    real(8), allocatable :: hist(:,:)
    allocate(hist(Nhist,2))

    print '(" ***** START calculation 1D histgram ******")'

    data_max = maxval(unit_bead)
    data_min = minval(unit_bead)
    data_ave = sum(unit_bead)/size(unit_bead)

    if ( hist_min(1) == 0.0d0 .and. hist_max(1) == 0.0d0 ) then
      print *, "   Using the margin parameter"
      hist_min(1) = data_min - hist_margin
      hist_max(1) = data_max + hist_margin
    else
      print *, "   Using the hist_X_min and hist_X_max"
    end if

    Dhist = (hist_max(1) - hist_min(1)) / dble(Nhist)

    print '("    X range max =", F13.6)', hist_max(1)
    print '("    X range min =", F13.6)', hist_min(1)
    print '("    Number hist =", I8)',    Nhist
    print '("    Delta hist  =", F13.6)', Dhist
  end subroutine calc_1Dhistogram_unit

  subroutine calc_unit_bond(atom1,atom2,Iunit)
    integer, intent(in) :: atom1, atom2, Iunit
    integer :: i,j,k
    do k = 1, TNstep
      do j = 1, Nbeads
        sij(:) = s(:,atom1,j,k) - s(:,atom2,j,k)
        sij(:) = sij(:) - nint(sij(:))
        rij(:) = matmul(sij(:),lattice(:,:))
        dij = dsqrt(sum(rij(:)*rij(:)))
        unit_bead(Iunit,j,k) = dij
      end do
    end do
  end subroutine calc_unit_bond

end subroutine pbhpo4



