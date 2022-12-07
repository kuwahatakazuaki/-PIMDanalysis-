subroutine rotation
  use input_parameter,  only: atom, atom_num, TNstep, save_beads, Nbeads, Natom, &
      FNameBinary1, graph_step, Nstart, Nstep, weight, r_ref, jobtype, label, &
      muon => atom_density, Nhyd, hyd
  use calc_parameter,   only: r
  use utility,          only: calc_deviation, calc_cumulative, get_rot_mat, lowerchr
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
    case(76)
      call save_cube
  end select

contains

  subroutine save_cube
    integer :: Uout
    integer, parameter :: Ndiv = 20
    real(8), parameter :: Ledge = 10.0d0
    real(8), parameter :: Bohr2Angs = 0.529177249
    real(8), parameter :: Angs2Bohr = 1.8897259886
    real(8), parameter :: margine = 1d-2
    real(8) :: grid(Ndiv,Ndiv,Ndiv)
    real(8) :: Lmin(3), Lmax(3)
    real(8) :: dL(3), base_vec(3,3)
    integer, allocatable :: coun(:,:,:)
    integer :: cx,cy,cz

    rnew(:,:,:,:) = rnew(:,:,:,:) * Angs2Bohr
    do i = 1, 3
    Lmin(i) = minval(rnew(i,muon,:,:)) - margine
    Lmax(i) = maxval(rnew(i,muon,:,:)) + margine
    end do
    dL(:) = (Lmax(:) - Lmin(:)) / dble(Ndiv)

print *, 'dL : ', dL(:)
print *, 'Lmin:', Lmin(:)
    base_vec(:,:) = 0.0d0
    do i = 1, 3
      base_vec(i,i) = dL(i)
    end do

    allocate(coun(3,Nbeads,TNstep))
    !coun(:,:,:,:) = (rnew(:,:,:,:)-spread( spread( spread(Lmin(:),dim=2,ncopies=Natom),dim=3,ncopies=Nbeads),dim=4,ncopies=TNstep))
    do k = 1, TNstep
      do j = 1, Nbeads
        coun(:,j,k) = int( ( rnew(:,muon,j,k)-Lmin(:) ) / dL(:) ) + 1
      end do
    end do
!do k = 1, TNstep
!  do j = 1, Nbeads
!    print *, i,j,k, coun(:,j,k)
!  end do
!end do
    grid(:,:,:) = 0.0d0
    do k = 1, TNstep
      do j = 1, Nbeads
          cx = coun(1,j,k)
          cy = coun(2,j,k)
          cz = coun(3,j,k)
          grid(cx,cy,cz) = grid(cx,cy,cz) + 1.0d0
      end do
    end do
    open(newunit=Uout,file='hyd.cube',status='replace')
      write(Uout,*) "commnet"
      write(Uout,*) "commnet"
      write(Uout,9999) Natom-Nhyd, Lmin(:)
      do i = 1, 3
        write(Uout,9999) Ndiv, base_vec(i,:)
      end do
      j = 1
      do i = 1, Natom-Nhyd
        if ( i == hyd(j) ) then
          j = j + 1
          cycle
        end if
        write(Uout,9999) atom2num(trim(label(i))), dble(i), &
                         [sum(rnew(1,i,:,:)),sum(rnew(2,i,:,:)),sum(rnew(3,i,:,:))]/dble(TNstep*Nbeads)
      end do
      do i = 1, Ndiv
        do j = 1, Ndiv
          do k = 1, Ndiv
            write(Uout,'(E13.5)',advance='no') grid(i,j,k)
            if ( mod(k,6) == 0 ) write(Uout,*)
          end do
          write(Uout,*)
        end do
      end do
    close(Uout)

  9998  format(I5,4F12.6)
  9999  format(I5,4F12.6)
  end subroutine save_cube

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

  function atom2num(cha) result(num)
    character(*) :: cha
    integer :: num
    cha = lowerchr(cha)
    num = 0
    if     ( trim(cha) == 'h' ) then
      num = 1
    elseif ( trim(cha) == 'li' ) then
      num = 3
    elseif ( trim(cha) == 'b' ) then
      num = 5
    elseif ( trim(cha) == 'c' ) then
      num = 6
    elseif ( trim(cha) == 'n' ) then
      num = 7
    elseif ( trim(cha) == 'o' ) then
      num = 8
    elseif ( trim(cha) == 'f' ) then
      num = 9
    else
      stop 'ERROR!! "atom2num" cannot chage '
    end if
  end function

  character(len=2) function itoc(i)
    integer :: i
    select case(i)
      case(1)
        itoc = 'H'
      case(3)
        itoc = 'Li'
      case(5)
        itoc = 'B'
      case(6)
        itoc = 'C'
      case(7)
        itoc = 'N'
      case(8)
        itoc = 'O'
      case(9)
        itoc = 'F'
    end select
  end function itoc

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

