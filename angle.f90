subroutine calc_angle
  use input_parameter,  only: atom, atom_num, out_hist, TNstep, Nfile, save_beads, FNameBinary1
  use calc_parameter,   only: data_step, data_beads
  use calc_histogram1D, only: calc_1Dhist
  use utility,          only: calc_deviation, calc_cumulative
  integer :: Ifile, step, k ! i=atom, j=beads, k=step
  character(len=128) :: out_angle, angle_name, out_cumulative
  real(8) :: data_max, data_min, data_ave, data_dev, data_err
!  integer :: atom1(Nfile),atom2(Nfile),atom3(Nfile),atom4(Nfile),atom5(Nfile)
!  atom1(:) = atom_num(1,:)
!  atom2(:) = atom_num(2,:)
!  atom3(:) = atom_num(3,:)
!  atom4(:) = atom_num(4,:)
!  atom5(:) = atom_num(5,:)

  write(angle_name,'(a,I0,"-",a,I0,"-",a,I0)') &
    trim(atom(atom_num(1,1))),atom_num(1,1), trim(atom(atom_num(2,1))),atom_num(2,1), trim(atom(atom_num(3,1))),atom_num(3,1)
  write(out_angle,'("angle_",a,".out")') trim(angle_name)
  write(out_cumulative, '("cumu_",a,".out")') trim(angle_name)
!  write(out_hist, '("hist_angle_",a,".out")') trim(angle_name)

  if ( trim(out_hist) == "0") write(out_hist, '("hist_",a,".out")') trim(angle_name)

    step = 0
    do Ifile = 1, Nfile
!      call calc_angle_sub(Ifile,atom1(Ifile),atom2(Ifile),atom3(Ifile),step)
      call calc_angle_sub(Ifile,atom_num(1,Ifile),atom_num(2,Ifile),atom_num(3,Ifile),step)
    end do
    data_max = maxval(data_beads)
    data_min = minval(data_beads)
    data_ave = sum(data_beads)/size(data_beads)
    call calc_deviation(data_dev, data_err)

block
  integer :: Ounit
  if ( save_beads .eqv. .True. ) then
    open(newunit=Ounit,file=FNameBinary1, form='unformatted', access='stream', status='replace')
      do step = 1, TNstep
        do i = 1, Nbeads
          write(Ounit) data_beads(i,step)
        end do
      end do
    close(Ounit)
  end if
end block

  open(23,file=out_angle,status='replace')
    write(23, '("# Angel Histgram of ",a)') trim(angle_name)
    write(23, '("# Maximum angle = ", F13.6)') data_max
    write(23, '("# Minimum angle = ", F13.6)') data_min
    write(23, '("# Average angle = ", F13.6)') data_ave
    write(23, '("# Deviation     = ", F13.6)') data_dev
    write(23, '("# St. erro      = ", F13.6)') data_err
    do k = 1, TNstep
      if (mod(k,10) == 0) then
        write(23,'(I7,F10.5)') k, data_step(k)
      end if
    end do
  close(23)

  print '("    Maximum angle = ", F13.6)', data_max
  print '("    Minimum angle = ", F13.6)', data_min
  print '("    Average angle = ", F13.6)', data_ave
  print '("    Deviation     = ", F13.6)', data_dev
  print '("    St. error     = ", F13.6)', data_err
  print '("    Data is saved in ",a,/)', '"'//trim(out_angle)//'"'
!  call calc_1Dhist ! you need "data_beads"
  call calc_1Dhist(out_hist_ex=out_hist) ! you need "data_beads"
  call calc_cumulative(out_cumulative)
end subroutine calc_angle

subroutine calc_angle_sub(Ifile,atom_temp1,atom_temp2,atom_temp3,step)
  use input_parameter,only: Nstart, Nstep, Nbeads
  use calc_parameter, only: r, data_beads, data_step
  use utility,        only: norm
  implicit none
  integer, intent(in) :: Ifile, atom_temp1, atom_temp2, atom_temp3
  integer, intent(inout) :: step
  integer :: j, k
  real(8) :: pi = atan(1.0d0)*4.0d0
  real(8) :: vec(3,2)

  do k = Nstart(Ifile), Nstep(Ifile)
    step = step+1
    do j = 1, Nbeads
      vec(:,1) = r(:,atom_temp1,j,step) - r(:,atom_temp2,j,step)
      vec(:,2) = r(:,atom_temp3,j,step) - r(:,atom_temp2,j,step)

      vec(:,1) = vec(:,1)/ norm(vec(:,1))
      vec(:,2) = vec(:,2)/ norm(vec(:,2))

      data_beads(j,step) = dacos(dot_product(vec(:,1),vec(:,2)))
      data_beads(j,step) = data_beads(j,step) *180.0/pi

!      print *, data_beads(j,step)
    end do
    data_step(step) = sum(data_beads(:,step))/dble(Nbeads)
  end do
end subroutine calc_angle_sub

