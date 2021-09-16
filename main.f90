! Analysing PIMD data
! last modified 2020/05/24
! program structure
! all of the data are stored in the 'data_beads', and this is analyzed according to the 'jobtype' option

program analysis
use input_parameter
use calc_parameter
use calc_centoroid
use calc_histogram1D
use calc_histogram2D
use mod_other_quantities
implicit none

! +++ Reading the input file +++
call read_input

! r(:,i,j,k) = r(xyz,atom,beads,step)
allocate(r(3,Natom,Nbeads,TNstep))
allocate(atom(Natom))
allocate(data_step(TNstep), source=0.0d0)
allocate(data_beads(Nbeads,TNstep), source=0.0d0)

! +++ Reading coordinate +++
select case(jobtype)
  case(60:69)
    print *, "Skip reading coor.xyz"
  case default
    call read_coor
end select

! Choose "job type"
!  1 : -1D histgram of bond         (atom1-atom2)
!  2 : -Angle histgram              (atom1-atom2-atom3)
!  3 : -Dihedral angle              (atom1-atom2-atom3-atom4)
!  9 : -1D histgram from External
! 11 : -Multi bond calc all
! 12 : -Multi bond sort
! 13 : -Multi bond diff              (atom1-atom2  -  atom3-atom4)
! 14 : -Multi bond diff              (atom1-atom2  +  atom3-atom4)
! 21 : -2D histogram_bond            (atom1-atom2 and atom3-atom4)
! 22 : -2D histogram_angle           (atom1-atom2 and atom3-atom4-atom5)
! 29 : -2D histogram from External   (Normal)
! 31 : -1D histogram for Centroid
! 32 : -2D histogram for Centroid
! 33 : -Angle histgram for Centroid
! 41 : -Dummy atom (X) for bond      (atom1-atomX)
! 42 : -Dummy atom (X) for angle     (atom1-atom2-atomX)
! 43 : -Dummy atom (X) for dihedral  (atom1-atom2-atomX-atom4)
! 51 : -Beads expansion   (all atoms)
! 52 : -Beads expansion   (export binary of atom1)
! 53 : -Beads expansion   (atom1 projected to atom2-atom3)
! 61 : -charge_analysis   (all atoms)
! 62 : -charge_analysis   (atom1)
! 63 : -dipole_analysis
! 71 : -projection        (atom1-atom2  T  atom3-atom4)
!!! 91 : -Specific purpose (Dihedral of NH4+(H2O))

select case(jobtype)
  case(1)
    call calc_bond
  case(2)
    call calc_angle
  case(3)
    call calc_dihedral
  case(9)
    call external_1Dhist
  case(11:14)
    call multi_bond
  case(21:22)
    call calc_2Dhist
  case(27:29)
    call external_2Dhits
  case(31:39)
    call calc_cent
  case(41:49)
    call dummy_atom
  case(51:54)
    call beads_expansion
  case(61:63)
    call other_quantities
  case(71)
    call projection
  case default
    stop 'ERROR!!! wrong "Job type" option'
end select

deallocate(r)
deallocate(atom)
deallocate(data_step)
deallocate(data_beads)
print '(a,/)', " Normal termination"
stop
end program analysis
