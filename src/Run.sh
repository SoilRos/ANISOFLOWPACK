make

np=2	# Number of processors to be used
InputDir=

mpiexec -np $np ./ANISOFLOPACK.exe -Input_Dir $InputDir