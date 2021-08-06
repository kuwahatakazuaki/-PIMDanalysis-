module input_parameter
  implicit none
  integer, save :: Natom, Nbeads, TNstep, Nhist, Nfile, Nbond
  integer, save :: jobtype
  integer, save :: definition_dummy, atom_dummy(2)
  integer, save :: graph_step = 10
  integer, save :: Usave = 20, Uprint = 6  ! setting for output ! output = 0
  integer, save, allocatable :: Nstep(:), Nstart(:) !, atom1(:), atom2(:), atom3(:), atom4(:), atom5(:)
  integer, save, allocatable :: atom_num(:,:) ! atom_num(Number,Nfile)
  integer, save, allocatable :: atom_multi(:,:)
  real(8), save :: hist_min(2)=0.0d0, hist_max(2)=0.0d0, hist_margin !, hist_max
  character(len=128) :: out_hist = "0"
  character(len=2), save, allocatable :: atom(:) ! the element of atom
  character(len=128), save, allocatable :: FileName(:)
  logical :: FIbinary
  logical :: Lfolding = .False.
  logical :: save_beads = .False.
  character(len=:), allocatable :: FNameBinary1, FNameBinary2
  character(:), allocatable :: input_file
  integer :: other_type
  character(:), allocatable :: other_path
  integer :: umbrella_type, umbrella_atom1, umbrella_atom2, umbrella_atom3
end module input_parameter

module calc_parameter
  real(8), save, allocatable :: r(:,:,:,:) ! r(xyz,i=atom,j=beads,k=step)
  real(8), save, allocatable :: data_step(:), data_beads(:,:) ! step(TNstep), beads(Nbeads,TNstep)
end module calc_parameter


