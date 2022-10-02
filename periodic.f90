module mod_periodic
  use input_parameter
  use calc_parameter, only: data_beads, data_step
  use calc_histogram1D
  implicit none
  integer, private :: Uinp, ierr
  integer :: atom1, atom2, atom3, atom4

contains

  subroutine periodic
    atom1 = atom_num(1,1)
    atom2 = atom_num(2,1)
    atom3 = atom_num(3,1)
    atom4 = atom_num(4,1)

    select case(jobtype)
      case(81)
        call RDF1
!      case(82)
!        call RDF2
    end select
  end subroutine periodic

! ++++++++++++++++++
! +++ Start RDF1 +++
! ++++++++++++++++++
  subroutine RDF1
    integer :: Istep, Ifile, i, j, k
    integer :: Uout
    character(len=128) :: out_name
    real(8) :: r21(3), r31(3), r41(3)
!    character(len=:), allocatable :: out_name
    write(out_name,'("outplane_",a,I0,"-",a,I0,"-"a,I0,"->",a,I0)') &
            & trim(atom(atom2)), atom2, trim(atom(atom1)), atom1, &
            & trim(atom(atom3)), atom3, trim(atom(atom4)), atom4
    print *, out_name

    Istep = 0
    do Ifile = 1, Nfile
      do k = Nstart(Ifile), Nstep(Ifile)
        Istep = Istep + 1
        do j = 1, Nbeads
          r21(:) = r(:, atom2,j,Istep) - r(:,atom1,j,Istep)
          r31(:) = r(:, atom3,j,Istep) - r(:,atom1,j,Istep)
          r41(:) = r(:, atom4,j,Istep) - r(:,atom1,j,Istep)
        end do
      end do
    end do

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

    open(newunit=Uout, file=out_name, status='replace')
    close(Uout)

  end subroutine RDF1
! ++++++++++++++++
! +++ End RDF1 +++
! ++++++++++++++++

end module mod_periodic




