subroutine calc_dihedral
!  use input_parameter, only: TNstep, atom, atom_num, graph_step, out_hist
!  If you set only keyword in input_parameter, you cannot use calc_dihedral_sub subroutine
!  I don't know why
  use input_parameter
  use calc_parameter,  only: data_step, data_beads
  use calc_histogram1D, only: calc_1Dhist
  use utility, only: calc_deviation, calc_cumulative, reblock_step
  implicit none
  integer :: i, k ! r(xyz,i=atom,j=beads,k=step)
  integer :: Ifile, step
  character(len=128) :: out_angle, angle_name, out_cumulative
  real(8) :: data_max, data_min, data_ave, data_dev

  write(angle_name,'(a,I0,"-",a,I0,"-",a,I0,"-",a,I0)') &
    trim(atom(atom_num(1,1))),atom_num(1,1), trim(atom(atom_num(2,1))),atom_num(2,1), &
    trim(atom(atom_num(3,1))),atom_num(3,1), trim(atom(atom_num(4,1))),atom_num(4,1)
  write(out_angle,'("dihe_",a,".out")') trim(angle_name)
  write(out_cumulative, '("cumu_",a,".out")') trim(angle_name)
  if ( trim(out_hist) == "0") write(out_hist, '("hist_dihe_",a,".out")') trim(angle_name)

    step = 0
    do Ifile = 1, Nfile
      call calc_dihedral_sub(Ifile, &
          atom_num(1,Ifile),atom_num(2,Ifile), &
          atom_num(3,Ifile),atom_num(4,Ifile), step)
    end do
    data_max = maxval(data_beads)
    data_min = minval(data_beads)
    data_ave = sum(data_beads)/size(data_beads)
    call calc_deviation(data_dev)

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
    do k = 1, TNstep
      if (mod(k,graph_step) == 0) then
        write(23,'(I7,F10.5)') k, data_step(k)
      end if
    end do
  close(23)

  print '("    Maximum angle = ", F13.6)', data_max
  print '("    Minimum angle = ", F13.6)', data_min
  print '("    Average angle = ", F13.6)', data_ave
  print '("    Deviation     = ", F13.6)', data_dev
  print '("    Data is saved in ",a,/)', '"'//trim(out_angle)//'"'
  call calc_1Dhist(out_hist_ex=out_hist)
  call reblock_step()
!  call calc_cumulative(out_cumulative)
end subroutine calc_dihedral


subroutine calc_dihedral_sub(Ifile, Tatom1, Tatom2, Tatom3, Tatom4, step)
  use input_parameter
  use calc_parameter, only: r, data_beads, data_step
  use utility, only: norm, outer_product, pi
  implicit none
  integer, intent(in) :: Ifile, Tatom1, Tatom2, Tatom3, Tatom4
  integer, intent(inout) :: step
  integer :: j, k ! i=atom, j=beads, k=step
!  real(8) :: pi = atan(1.0d0)*4.0d0
  real(8) :: r12(3), r23(3), r34(3)
  real(8) :: vec123(3), vec234(3)

  do k = Nstart(Ifile), Nstep(Ifile)
    step = step+1
    do j = 1, Nbeads
      r12(:) = r(:,Tatom1,j,step) - r(:,Tatom2,j,step)
      r23(:) = r(:,Tatom2,j,step) - r(:,Tatom3,j,step)
      r34(:) = r(:,Tatom3,j,step) - r(:,Tatom4,j,step)

      vec123(:) = outer_product(r12(:),r23(:))
      vec234(:) = outer_product(r23(:),r34(:))

      vec123(:) = vec123(:) / norm(vec123(:))
      vec234(:) = vec234(:) / norm(vec234(:))

      data_beads(j,step) = dacos( dot_product(vec123(:), vec234(:)) )
      data_beads(j,step) = data_beads(j,step) *180.0/pi
    end do
    data_step(step) = sum(data_beads(:,step))/dble(Nbeads)
  end do
end subroutine calc_dihedral_sub

