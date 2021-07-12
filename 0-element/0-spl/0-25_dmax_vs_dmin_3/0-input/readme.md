# How to create a sample file for LIGGGHTS?

To construct a sample, LIGGGHTS has the capability of inserting particles with various sizes. But based on the work of Mutabaruka et al. (2019), we can take advantange of the cpp code __pse_43.cpp__ to automatically generate a particle assembly following a specified particle size distribution (PSD). The procedures to achive this goal are:

1. Go to folder _data/spl/_, and compile __pse_43.cpp__ to generate the executable __pse_43__:

   > clang++ pse_43.cpp -o pse_43

   Then open the file __inputParameters.dat__ to make the corresponding changes according to Mutabaruka et al. (2019). Then run

   > ./pse_43 inputParameters.dat

   to generate a list of files. One of them named __sample.spl__ is the sample file.

2. Get out of the folder _spl_, and open the python document __main_psd.py__. Run the function __cvd_to_psd()__ in the __main()__ function and it will allow us to compare the sample PSD with the experimental one (the experiment corresponds to Ottawa-F65 sand in this example). Run the function __spl_to_data()__ and it will generate the sample file compatible with LIGGGHTS, namely __Ottawa-F65.spl__ in the current folder.

3. To run triaxial test in tri-periodic boundary conditions, one can neglect the part of how to generate mesh-plate file since there is no need using it. In simple shear test, however, one way of doing this test is to introduce mesh-plate files as the top and bottom walls. In macOS or linux, one can install __gmsh__ to create the mesh. Now go to the folder _mesh_, and modify the file __plane.geo__ to have a size comparable with the dimensions of Ottawa-F65.spl (one can open __Ottawa-F65.spl__ to have a sense of these dimensions). Then use the command in the terminal:

   > gmsh -2 plane.geo -o plane.stl -format stl

   to generate the output __plane.stl__ that can be used for LIGGGHTS.