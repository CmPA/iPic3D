Author: Alexander E. Vapirev, KU Leuven, Afdeling Plasma-astrofysica

Note  : All *-2D scripts and folders are for the Parsek2D case. They can be easily merged with a simple 2D/3D-switch option in one script/code but are left separate for clarity (and due to laziness).

Note  : If the gnuplot generating scripts are used, then a little knowledge of gnuplot is helpful.

I. Data processing and extraction from iPIC3D simulation VirtualSatellite* data files
-------------------------------------------------------------------------------------

The CPP* folders contain the c++ code for various cases of data extraction from the VirtualTraces satellite files in iPIC3D. Each c++ code can be used as a standalone code, but for user convenience the shell scripts have been created to govern a more flexible input/output. Once the data is extracted, a set of shell gnuplot script generators is provided for further plotting and data anlysis in the name of the brilliant science one is doing.

1. process-virtual-satellite-files-point.sh

The script governs the c++ code 'cpp-code-point/process-virtual-sat-point.cpp'. The c++ code extracts data for a point satellite which is closest to a set of user defined coordinates (x,y,z). The c++ code searches through all the VirtualSattelite* files and finds the best match. Then it extracts the data and also computes some extra quantities (see the c++ code). The out put is written in a text file with a header explanation and columns of data. Each column corresponds to a specific quantity and each row is the simulation step/cycle.

2. process-virtual-satellite-files-block-single-var.sh

This script extracts the data for a specified trace in X,Y, or Z direction. The 'block-single-var' in the names stems from the fact that each quantity is written a separate file in the form of a block matrix - one format which gnuplot likes when plotting 3D data. In that case the column numbers are X, and row numbers are Y provided (as it is in this case) that the grid has constant dx and dy. In the gnuplot script generator later there is an option to scale X and Y to whatever is needed.

3. process-virtual-satellite-files-block.sh

This script extracts the data for a specified trace in X,Y, or Z direction - the format in each row is {X0 Y0 Z00_1 ... Z00_n} .. {Xi Yj Zij_1 ... Zij_n} (n quantities are extracted in a single block). A plotting script for 'process-virtual-satellite-files-block.sh' has not been created since the respective extracted data file is huge and requires enormous computing resources to plot it. Code is not kept up-to-date.



II. Gnuplot scripts
-------------------

After extracting the data, one can use the generate-gnuplot* shell scripts to plot it.

II.A.

1. generate-gnuplot-code-point.sh

Generates a gnuplot script to plot data processed with 'process-virtual-satellite-files-point.sh' from (I.1).

2. generate-gnuplot-code-block-single-var.sh

Generates a gnuplot script to plot data processed with 'process-virtual-satellite-files-block-single-var.sh' from (I.2).

II.B.

Sometimes an iPIC3D run is split into 2 runs, the second run being a continuation of the first one via a 'restart' file. In case one would like to plot quantities continuously from first run to second run, the following script are provided.

1. generate-gnuplot-code-stitch-point.sh

Stitch quantities from first and second run for a 'point' satellite data probe.

2. generate-gnuplot-code-stitch-block-single-var.sh

Stitch quantities from first and second run for a 'trace' satellite data probe.

II.C.

Script generators are developed for plotting the FFT after the VirtualSatellite data has been extracted. The gnuplot generator scripts first utilize a c++ code which computes the FFT of the data. Then a gnuplot script is generated and plots the data from the output FFT data file.

1. generate-gnuplot-code-point-fft.sh

Compiles and executes a c++ code to compute the FFT form the output data obtained by 'process-virtual-satellite-files-point.sh'. Then it plots the FFTs.


