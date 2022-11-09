subroutine rotation
  use input_parameter,  only: atom, atom_num, TNstep, save_beads, Nbeads, Natom, &
      FNameBinary1, graph_step, Nstart, Nstep, weight, r_ref, jobtype, label
  use calc_parameter,   only: r
  use utility,          only: calc_deviation, calc_cumulative, get_rot_mat
  implicit none
  real(8), parameter :: pi = 4.0d0*atan(1.0d0)
  integer, parameter :: N = 4
  integer, parameter :: liwork = 5*N+3,  lwork = 2*N*N+6*N+1
  real(8) :: work(lwork)
  integer :: iwork(liwork)
  integer :: i, j, k, xyz, info
  real(8), allocatable :: rnew(:,:,:,:)
  real(8) :: rot(3,3), qua(4), rave(3), Tweight
  real(8) :: matA(4,4), matB(4,4), eigval(4)

  allocate(rnew(3,Natom,Nbeads,TNstep))
  Tweight = sum(weight(:))
  do xyz = 1, 3
    rave(xyz) = dot_product(r_ref(xyz,:), weight(:)) / (Tweight)
  end do
  r_ref(:,:) = r_ref(:,:) - spread(rave(:), dim=2, ncopies=Natom)

  do k = 1, TNstep
    rave(:) = 0.0d0
    do j = 1, Nbeads
      do xyz = 1, 3
        rave(xyz) = rave(xyz) + dot_product(r(xyz,:,j,k), weight(:))
      end do
    end do
    rave(:) = rave(:) / (Tweight * dble(Nbeads))
    r(:,:,:,k) = r(:,:,:,k) - spread( spread(rave(:),dim=2,ncopies=Natom),dim=3,ncopies=Nbeads )
  end do

  do k = 1, TNstep
    matB(:,:) = 0.0d0
    do j = 1, Nbeads
      do i = 1, Natom
        matA(:,:) = make_matA(r_ref(:,i)+r(:,i,j,k),r_ref(:,i)-r(:,i,j,k))
        matB(:,:) = matB(:,:) + weight(i)*matmul(transpose(matA),matA)
      end do
    end do
    matB(:,:) = matB(:,:) / (Tweight*dble(Nbeads))

    call dsyevd('V', 'U', N, matB, N, eigval, work, lwork, iwork, liwork, info)
 !   print *, eigval(:)    ! ??????
    qua(:) = matB(:,1)
    rot(:,:) = get_rot_mat(qua(:))
    do j = 1, Nbeads
      do i = 1, Natom
        rnew(:,i,j,k) = matmul(rot(:,:), r(:,i,j,k))
      end do
    end do
  end do

  select case(jobtype)
    case(75)
      call save_movie
!    case(76)
!      call save_cube
  end select

contains

  subroutine save_movie
    integer :: Uout
    open(newunit=Uout,file='vmd.xyz',status='replace')
      !do k = 1, TNstep, graph_step
      do k = 1, TNstep
        write(Uout,'(I10)') Natom*Nbeads
        write(Uout,'(I10)') k
        do j = 1, Nbeads
          do i = 1, Natom
            write(Uout,9999) label(i), rnew(:,i,j,k)
          end do
        end do
      end do
    close(Uout)
9999 format(a,4F11.7)
  end subroutine save_movie

  function make_matA(x,y) result(mat)
    real(8) :: x(3), y(3)
    real(8) :: mat(4,4)
    mat(1,:) = [0.d0,  -y(1), -y(2), -y(3)]
    mat(2,:) = [y(1),   0.d0, -x(3),  x(2)]
    mat(3,:) = [y(2),   x(3),  0.d0, -x(1)]
    mat(4,:) = [y(3),  -x(2),  x(1),  0.d0]
  end function make_matA

!  integer :: i, j, k, Ifile ! i=atom, j=beads, k=step, l=hist, Ifile=file
!  integer :: step
!  character(len=128) :: out_name
!  real(8) :: data_max, data_min, data_ave, data_dev, data_err
!  real(8) :: e(3), r2(3), rt(3)
!  integer :: atom1, atom2, atom3, atom4
!
!!  write(bond_name, '(a,I0,"-",a,I0)')   trim(atom(atom_num(1,1))), atom_num(1,1), trim(atom(atom_num(2,1))), atom_num(2,1)
!!  write(out_bond, '("bond_",a,".out")') trim(bond_name)
!!  write(out_cumulative, '("cumu_",a,".out")') trim(bond_name)
!
!  step = 0
!  do Ifile = 1, Nfile
!    atom1 = atom_num(1,Ifile)
!    atom2 = atom_num(2,Ifile)
!    atom3 = atom_num(3,Ifile)
!    atom4 = atom_num(4,Ifile)
!
!    do k = Nstart(Ifile), Nstep(Ifile)
!    step = step + 1
!    do j = 1, Nbeads
!      e(:) = r(:,atom4,j,step) - r(:,atom3,j,step)
!      e(:) = e(:) / dsqrt( dot_product(e(:),e(:)))
!
!      r2(:) = r(:,atom2,j,step) - r(:,atom1,j,step)
!      rt(:) = r2(:) - dot_product(r2(:),e(:)) * e(:)
!
!      data_beads(j,step) = dsqrt( dot_product(rt(:),rt(:)) )
!    end do
!    data_step(step) = sum(data_beads(:,step)) / dble(Nbeads)
!    end do
!  end do
!
!  data_max = maxval(data_beads)
!  data_min = minval(data_beads)
!  data_ave = sum(data_beads)/size(data_beads)
!  call calc_deviation(data_dev, data_err)
!
!block
!  integer :: Ounit
!  if ( save_beads .eqv. .True. ) then
!    open(newunit=Ounit,file=FNameBinary1, form='unformatted', access='stream', status='replace')
!      do step = 1, TNstep
!        do i = 1, Nbeads
!          write(Ounit) data_beads(i,step)
!        end do
!      end do
!    close(Ounit)
!  end if
!end block

!  open(Usave, file=out_bond, status='replace')
!    do Ifile = 1, Nfile
!      write (Usave,'(" # From file",I0, " : "a)') Ifile,trim(bond_name)
!    end do
!    write(Usave,'(a,F13.6)') " # Maximum bond  = ", data_max
!    write(Usave,'(a,F13.6)') " # Minimum bond  = ", data_min
!    write(Usave,'(a,F13.6)') " # Average bond  = ", data_ave
!    write(Usave,'(a,F13.6)') " # St. deviation = ", data_dev
!    write(Usave,'(a,F13.6)') " # St. error     = ", data_err
!    do k = 1, TNstep
!      if (mod(k,graph_step) == 0) then
!        write(Usave,'(I7,F10.5)') k, data_step(k)
!      end if
!    end do
!  close(Usave)
!
!  write(Uprint,*) "*****START calculating bond length*****"
!  do Ifile = 1, Nfile
!    write (Uprint,'("    From file",I0, " : "a)') Ifile, trim(bond_name)
!  end do
!  write(Uprint, '("    Maximum bond =", F13.6)') data_max
!  write(Uprint, '("    Minimum bond =", F13.6)') data_min
!  write(Uprint, '("    Average bond =", F13.6)') data_ave
!  write(Uprint, '("    St. deviation=", F13.6)') data_dev
!  write(Uprint, '("    St. error    =", F13.6)') data_err
!  write(Uprint,*) ""
!  call calc_1Dhist(out_hist_ex=out_hist)
end subroutine rotation

