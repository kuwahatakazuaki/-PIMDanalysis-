
module calc_histogram2D
  use input_parameter, &
      only: Natom, Nbeads, TNstep, Nhist, Nfile, Nbond, &
            jobtype, atom, atom_num, out_hist, Usave, &
            hist_min, hist_max, FNameBinary1, FNameBinary2, &
            hist_margin
  use calc_parameter, only: r, data_beads
  implicit none
  real(8), save :: hist2D_min(2), hist2D_max(2), Dhist(2)
  real(8), save :: hist2D_ave(2), bond_dev2D(2)
  real(8), save, allocatable :: hist_data(:,:), hist_axis(:,:)
  real(8), save, allocatable :: hist2D_beads(:,:,:)
  real(8), allocatable  :: hist_data1D(:)
  integer, private :: i, j, k, l, Ifile, step! i=atom, j=beads, k=step, l=hist, Ifile=file
contains

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_2Dhist
    integer :: temp1(2,Nfile), temp2(2,Nfile)
    character(len=128) name_axis(2)
    character :: axis(2) = (/"X","Y"/)

    print '(" ***START calculation 2D histgram***")'

    select case(jobtype)
      case(21)
        write(name_axis(1),'(a,I0,"-",a,I0)') trim(atom(atom_num(1,1))),atom_num(1,1), trim(atom(atom_num(2,1))),atom_num(2,1)
        write(name_axis(2),'(a,I0,"-",a,I0)') trim(atom(atom_num(3,1))),atom_num(3,1), trim(atom(atom_num(4,1))),atom_num(4,1)
        write(out_hist,'("hist_",a,"and",a,".out")') trim(name_axis(1)), trim(name_axis(2))
      case(22)
        write(name_axis(1),'(a,I0,"-",a,I0)') trim(atom(atom_num(1,1))),atom_num(1,1), trim(atom(atom_num(2,1))),atom_num(2,1)
        write(name_axis(2),'(a,I0,"-",a,I0,"-",a,I0)') &
          trim(atom(atom_num(3,1))),atom_num(3,1), trim(atom(atom_num(4,1))),atom_num(4,1),trim(atom(atom_num(5,1))),atom_num(5,1)
        write(out_hist,'("hist_",a,"and",a,".out")') trim(name_axis(1)), trim(name_axis(2))
      case default
        stop 'ERROR!!! wrong "Job type" option'
    end select

    allocate(hist_data(Nhist,Nhist))
    allocate(hist_axis(Nhist,2))
    allocate(hist2D_beads(Nbeads,TNstep,2))

    if (jobtype == 21) then ! 2D hist for bond and bond
      temp1(1,:)=atom_num(1,:); temp1(2,:)=atom_num(3,:)
      temp2(1,:)=atom_num(2,:); temp2(2,:)=atom_num(4,:)
      do i = 1, 2 ! 1:X-axis, 2:Y-axis
        step=0
        do Ifile = 1, Nfile
          call calc_bond_sub(Ifile,temp1(i,Ifile),temp2(i,Ifile),step)
        end do
        hist2D_beads(:,:,i) = data_beads(:,:) ! data_beads(Nbeads,TNstep)
        hist2D_max(i) = maxval(data_beads)
        hist2D_min(i) = minval(data_beads)
        hist2D_ave(i) = sum(data_beads)/size(data_beads)
      end do
    elseif (jobtype == 22) then ! 2D hist for bond and angle
      ! calculation of bond length
      step=0
      do Ifile = 1, Nfile
        call calc_bond_sub(Ifile,atom_num(1,Ifile),atom_num(2,Ifile),step)
      end do
      hist2D_beads(:,:,1) = data_beads(:,:) ! data_beads(Nbeads,TNstep)
      hist2D_max(1) = maxval(data_beads)
      hist2D_min(1) = minval(data_beads)
      hist2D_ave(1) = sum(data_beads)/size(data_beads)

      ! calculation of angle
      step = 0
      do Ifile = 1, Nfile
        call calc_angle_sub(Ifile,atom_num(3,Ifile),atom_num(4,Ifile),atom_num(5,Ifile),step)
      end do
      hist2D_beads(:,:,2) = data_beads(:,:) ! data_beads(Nbeads,TNstep)
      hist2D_max(2) = maxval(data_beads)
      hist2D_min(2) = minval(data_beads)
      hist2D_ave(2) = sum(data_beads)/size(data_beads)
    endif

    call calc_2Dhist_sub(hist2D_beads)
  ! --- START Calculating 2D histgram --- !

    open(Usave, file=trim(out_hist), status='replace')
      do i = 1, 2
        write(Usave,'(" # ",a," axis is ",a)') axis(i),trim(name_axis(i))
        write(Usave,'(" # Maximum hist  ", F13.6)')  hist2D_max(i)
        write(Usave,'(" # Minimum hist  ", F13.6)')  hist2D_min(i)
        write(Usave,'(" # Average hist  ", F13.6)')  hist2D_ave(i)
      end do

      write(Usave,*) "# Hist parameter is as follows"
      do i = 1, 2
        write(Usave, '(" # Range max  ",a,F13.6)') axis(i), hist_max(i)
        write(Usave, '(" # Range min  ",a,F13.6)') axis(i), hist_min(i)
        write(Usave, '(" # Delta hist ",a,F13.6)') axis(i), Dhist(i)
      end do
      write(Usave, '(" # Number hist   ", I6)')     Nhist

      do i = 1, Nhist
        do j = 1, Nhist
          write(Usave,'(2F11.5, E14.5)') hist_axis(i,1), hist_axis(j,2), hist_data(i,j)
        end do
        write(Usave,*) ""
      end do
    close(Usave)

    do i = 1, 2
      print '("    Range max  ",a,F13.6)', axis(i), hist_max(i)
      print '("    Range min  ",a,F13.6)', axis(i), hist_min(i)
      print '("    Delta hist ",a,F13.6)', axis(i),    Dhist(i)
    end do
    print '("    Data is saved in ",a)', '"'//trim(out_hist)//'"'
    print '(" *****END 2D histgram*****",/)'

    deallocate(hist_data)
    deallocate(hist_axis)
    deallocate(hist2D_beads)

  ! --- END Calculating 2D histgram --- !
  end subroutine calc_2Dhist
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_2Dhist_sub(hist2D_beads)
  !  implicit none
    integer :: coun(2)
    real(8), intent(in) :: hist2D_beads(Nbeads,TNstep,2)
    character :: axis(2) = (/"X","Y"/)

    do i = 1, 2
      if ( hist_min(i) == 0.0d0 .and. hist_max(i) == 0.0d0 ) then
        print '(a,a)', "    Using the margin parameter to ", axis(i)
        hist_min(i) = hist2D_min(i) - hist_margin
        hist_max(i) = hist2D_max(i) + hist_margin
      else
        print '("    Using the hist_min and hist_max")'
      end if

      Dhist(i) = (hist_max(i) - hist_min(i)) / dble(Nhist)
    end do

    do i = 1, 2
      do l = 1, Nhist
        hist_axis(l,i) = hist_min(i) + Dhist(i) * dble(l)
      end do
    end do

! +++++++ HERE ++++++++
! +++++++ hist_data1D ++++++
    hist_data(:,:) = 0.0d0
    do j = 1, TNstep
      do k = 1, Nbeads
        do i = 1, 2

          do l = 1, Nhist
            if ( hist_axis(l,i) >= hist2D_beads(k,j,i) ) then
              coun(i) = l
              exit
            else
              cycle
            end if
          end do

        end do
        hist_data(coun(1),coun(2)) = hist_data(coun(1),coun(2)) + 1.0d0
      end do
    end do

    do i = 1, 2
      hist_axis(:,i) = hist_axis(:,i) - 0.5d0 * Dhist(i)
    end do
    hist_data(:,:) = hist_data(:,:) / (dble(TNstep*Nbeads)*Dhist(1)*Dhist(2))
  end subroutine calc_2Dhist_sub
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine external_2Dhits
    implicit none
    character :: axis(2) = (/"X","Y"/)
    integer :: Nunit

    block ! Check the number of line between Finput1 and Finput2
      integer :: Nline1, Nline2, Nline
      real(8) :: dummyD
      open(newunit=Nunit, file=FNameBinary1, form='unformatted', access='stream', status='old', err=901)
        Nline = 0
        do
          read(Nunit, end=800) dummyD
          Nline = Nline + 1
        end do
        800 continue
      close(Nunit)
      Nline1 = Nline

      open(newunit=Nunit, file=FNameBinary2, form='unformatted', access='stream', status='old', err=902)
        Nline = 0
        do
          read(Nunit, end=801) dummyD
          Nline = Nline + 1
        end do
        801 continue
      close(Nunit)
      Nline2 = Nline

      if ( Nline1 /= Nline2 ) then
        print *, "The number of lines are different between ", FNameBinary1, " and ", FNameBinary2
        stop "ERROR!!"
      end if
    end block

    allocate(hist_data(Nhist,Nhist))
    allocate(hist_data1D(Nhist))
    allocate(hist_axis(Nhist,2))
    allocate(hist2D_beads(Nbeads,TNstep,2))
    open(newunit=Nunit, file=FNameBinary1, form='unformatted', access='stream', status='old', err=901)
      do j = 1, TNstep
        do i = 1, Nbeads
          read(Nunit) hist2D_beads(i,j,1)
        end do
      end do
    close(Nunit)

    open(newunit=Nunit, file=FNameBinary2, form='unformatted', access='stream', status='old', err=901)
      do j = 1, TNstep
        do i = 1, Nbeads
          read(Nunit) hist2D_beads(i,j,2)
        end do
      end do
    close(Nunit)

    do i = 1, 2
      hist2D_max(i) = maxval(hist2D_beads(:,:,i))
      hist2D_min(i) = minval(hist2D_beads(:,:,i))
    end do
    call calc_2Dhist_sub(hist2D_beads)

    out_hist='hist_2D_external.out'
    open(Usave, file=trim(out_hist), status='replace')
      do i = 1, 2
        write(Usave,'(" # Maximum hist  ", F13.6)')  hist2D_max(i)
        write(Usave,'(" # Minimum hist  ", F13.6)')  hist2D_min(i)
        write(Usave,'(" # Average hist  ", F13.6)')  hist2D_ave(i)
      end do

      write(Usave,*) "# Hist parameter is as follows"
      do i = 1, 2
        write(Usave, '(" # Range max  ",a,F13.6)') axis(i), hist_max(i)
        write(Usave, '(" # Range min  ",a,F13.6)') axis(i), hist_min(i)
        write(Usave, '(" # Delta hist ",a,F13.6)') axis(i), Dhist(i)
      end do
      write(Usave, '(" # Number hist   ", I6)')     Nhist

      do i = 1, Nhist
        do j = 1, Nhist
          write(Usave,'(2F11.5, E14.5)') hist_axis(i,1), hist_axis(j,2), hist_data(i,j)
        end do
        write(Usave,*) ""
      end do
    close(Usave)

    do i = 1, 2
      print '("    Range max  ",a,F13.6)', axis(i), hist_max(i)
      print '("    Range min  ",a,F13.6)', axis(i), hist_min(i)
      print '("    Delta hist ",a,F13.6)', axis(i),    Dhist(i)
    end do
    print '("    Data is saved in ",a)', '"'//trim(out_hist)//'"'
    print '(" *****END 2D histgram*****",/)'

  return
  901 stop "ERROR!! Tere is no binary input1 file"
  902 stop "ERROR!! Tere is no binary input2 file"
  end subroutine external_2Dhits
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module calc_histogram2D


