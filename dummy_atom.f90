subroutine dummy_atom
use input_parameter, &
    only: Natom, Nbeads, TNstep, atom, jobtype, atom_num, &
          definition_dummy, atom_dummy, save_beads, FNameBinary1
use calc_parameter, only: r, data_beads
implicit none
integer :: i, j, k, step
integer :: Natom_back
real(8) :: r_new(3,Natom+1,Nbeads,TNstep)  ! r_new(:,i,j,k)
real(8) :: r_X(3)        ! dummy atom
real(8) :: r_temp(3,2)  ! r_temp(xyz,num)
character(len=2) :: atom_back(Natom)


! --- 1: middle point, 2: difference ---
select case(definition_dummy)
  case(1)  ! middle point
    do k = 1, TNstep
      do j = 1, Nbeads
        do i = 1, Natom
          r_new(:,i,j,k) = r(:,i,j,k)
        end do
          r_temp(:,1) = r(:,atom_dummy(1),j,k)
          r_temp(:,2) = r(:,atom_dummy(2),j,k)
          r_X(:) = 0.5 * (r_temp(:,1) + r_temp(:,2))
          r_new(:,Natom+1,j,k) = r_X(:)
      end do
    end do
  case(2)  ! difference
    stop 'Not yet'
  case default
    stop 'ERROR!!! wrong "definition_dummy" option'
end select

atom_back(:) = atom(:)
deallocate(r)
deallocate(atom)
allocate(r(3,Natom+1,Nbeads,TNstep))
allocate(atom(Natom+1))
r(:,:,:,:) = r_new(:,:,:,:)
atom(1:Natom) = atom_back(:)
atom(Natom+1) = "X"


Natom_back = Natom
Natom = Natom + 1

select case(jobtype)
  case(41)
    atom_num(2,:) = Natom
    call calc_bond
  case(42)
    atom_num(3,:) = Natom
    call calc_angle
  case(43)
    atom_num(3,:) = Natom
    call calc_dihedral
  case default
    stop 'ERROR!!! wrong "definition_dummy" option'
end select

Natom = Natom_back

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

return

contains




end subroutine dummy_atom

