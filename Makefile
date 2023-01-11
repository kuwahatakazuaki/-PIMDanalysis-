program = run.exe
# +++ gfortran +++
fc = gfortran
fcopt =  -O2 -pipe -lblas -llapack
#fcopt =  -Wall -O3 -fbacktrace -fbounds-check -lblas -llapack 
#fcopt = -Wall -O3
# +++ End gfortran +++

# +++ ifort +++
#fc = ifort
#fcopt =  -CB -traceback -fpe0
#fcopt =  -warn all -traceback
# +++ End ifort +++
objs = \
parameters.o       \
utility.o          \
hist1D.o           \
hist2D.o           \
other_quantities.o  \
special_case.o      \
beads_expansion.o  \
periodic.o         \
angle.o            \
bond.o             \
cent.o             \
main.o             \
multi_bond.o       \
read_coor.o        \
read_inp.o         \
dummy_atom.o       \
dihedral.o         \
projection.o         \
rotation.o         \
pbhpo4.o           \

module =              \
calc_centoroid.mod    \
calc_histogram1d.mod  \
calc_histogram2d.mod  \
calc_parameter.mod    \
input_parameter.mod   \
utility.mod           \
mod_other_quantities.mod \
mod_special_case.mod    \
mod_periodic.mod    \

%.mod : %.f90 %.o
	@true

$(program): $(objs)
	@echo
	$(fc) $(objs) -o $@ $(fcopt)
	cp $@ ../
	@echo -e '\e[34m Noraml termination!!!\e[m\n'
#	$(fc) $(fcopt) $(objs) -o $@

%.o : %.f90
	@echo
	@echo ' << Compiling >>' '"'$<'"'
	$(fc) $(fcopt) -c $< -o $@


clean:
	rm -f *.o *.mod $(program)

install: $(objs)
	$(fc) $(objs) -o $(program) $(fcopt)
	cp $(program) /Users/kuwahatakazuaki/Program/bin/PIMDanalysis
# 	cp $(program) /Users/kuwahatakazuaki/PIMD/Analysis/Program/PIMDanalysis

#.SUFFIXES: .f90 .o
#.f90.o:
#	@echo
#	@echo ' << Compiling >>' '"'$<'"'

#	$(fc) $(fcopt) -c $< -o $@
