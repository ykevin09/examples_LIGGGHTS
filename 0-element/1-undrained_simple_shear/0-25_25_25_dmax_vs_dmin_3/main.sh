start=`date +%s`
#
#mpirun -np 2 /Users/myang/Documents/1-Programming/LIGGGHTS-PUBLIC/src/lmp_auto < 0-iso_comp.lmp
#mpirun -np 2 /Users/myang/Documents/1-Programming/LIGGGHTS-PUBLIC/src/lmp_auto < 1-iso_comp.lmp
mpirun -np 2 /Users/myang/Documents/1-Programming/LIGGGHTS-PUBLIC/src/lmp_auto < 2-iso_comp_bp.lmp
#mpirun -np 2 /Users/myang/Documents/1-Programming/LIGGGHTS-PUBLIC/src/lmp_auto < 3-iso_comp_bp.lmp
#mpirun -np 2 /Users/myang/Documents/1-Programming/LIGGGHTS-PUBLIC/src/lmp_auto < 4-undrained_simple_shear_bp.lmp
#mpirun -np 2 /Users/myang/Documents/1-Programming/LIGGGHTS-PUBLIC/src/lmp_auto < 5-undrained_simple_shear_bp.lmp
end=`date +%s`

runtime=$((end-start))

echo "The running time is ${runtime}s"
echo "The running time is ${runtime}s"
echo "The running time is ${runtime}s"
echo "The running time is ${runtime}s"
echo "The running time is ${runtime}s"
