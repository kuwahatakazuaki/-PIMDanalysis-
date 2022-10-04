module mod_special_case
  use input_parameter
  use calc_parameter, only: data_beads, data_step
  use calc_histogram1D
  use utility
  implicit none
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
    integer :: Istep, Ifile, i, j, k
    character(len=128) :: out_name
    real(8) :: r21(3), r31(3), r41(3), rt(3)
    real(8) :: data_max, data_min, data_ave, data_dev, data_err
!    character(len=:), allocatable :: out_name
    write(out_name,'("outplane_",a,I0,"-",a,I0,"-"a,I0,"to",a,I0)') &
            & trim(atom(atom2)), atom2, trim(atom(atom1)), atom1, &
            & trim(atom(atom3)), atom3, trim(atom(atom4)), atom4
    if ( trim(out_hist) == "0") write(out_hist, '("hist_",a,".out")') trim(out_name)

    Istep = 0
    do Ifile = 1, Nfile
      do k = Nstart(Ifile), Nstep(Ifile)
        Istep = Istep + 1
        do j = 1, Nbeads
          r21(:) = r(:,atom2,j,Istep) - r(:,atom1,j,Istep)
          r31(:) = r(:,atom3,j,Istep) - r(:,atom1,j,Istep)
          r41(:) = r(:,atom4,j,Istep) - r(:,atom1,j,Istep)

          r21(:) = r21(:) / norm( r21(:) )
          r31(:) = r31(:) / norm( r31(:) )
          r41(:) = r41(:) / norm( r41(:) )

          rt(:) = outer_product(r21(:),r31(:))
          data_beads(j,Istep) = ( pi - acos( dot_product(r41(:),rt(:)) ) ) * 180.0 / pi
        end do
        data_step(Istep) = sum(data_beads(:,Istep))/dble(Nbeads)
      end do
    end do
    data_max = maxval(data_beads)
    data_min = minval(data_beads)
    data_ave = sum(data_beads)/size(data_beads)
    call calc_deviation(data_dev, data_err)

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

    open(Usave, file=trim(out_name)//'.out', status='replace')
      do Ifile = 1, Nfile
        write (Usave,'(" # From file",I0, " : "a)') Ifile,trim(out_name)
      end do
      write(Usave,'(a,F13.6)') " # Maximum bond  = ", data_max
      write(Usave,'(a,F13.6)') " # Minimum bond  = ", data_min
      write(Usave,'(a,F13.6)') " # Average bond  = ", data_ave
      write(Usave,'(a,F13.6)') " # St. deviation = ", data_dev
      write(Usave,'(a,F13.6)') " # St. error     = ", data_err
      do k = 1, TNstep
        if (mod(k,graph_step) == 0) then
          write(Usave,'(I7,F10.5)') k, data_step(k)
        end if
      end do
    close(Usave)

    write(Uprint,*) "*****START calculating bond length*****"
    do Ifile = 1, Nfile
      write (Uprint,'("    From file",I0, " : "a)') Ifile, trim(out_name)
    end do
    write(Uprint, '("    Maximum bond =", F13.6)') data_max
    write(Uprint, '("    Minimum bond =", F13.6)') data_min
    write(Uprint, '("    Average bond =", F13.6)') data_ave
    write(Uprint, '("    St. deviation=", F13.6)') data_dev
    write(Uprint, '("    St. error    =", F13.6)') data_err
    write(Uprint,*) ""
    call calc_1Dhist(out_hist_ex=out_hist) ! you need "data_beads"

  end subroutine out_plane
! +++++++++++++++++++++
! +++ End out_plane +++
! +++++++++++++++++++++

end module mod_special_case




