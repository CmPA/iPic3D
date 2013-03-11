Plotting script library to extract data from already existing VTK iPic3D files and performing FFT in space and time.


Main script - it generates some supplementary files and runs the *.cpp code:

fft2d-slice-of-vtk.sh

c++ code - there are some explanations on the FFT too:

fft2d-slice-of-vtk.cpp

Output files:

--- suplementary (automatically generated)
indices.txt
user_inputs.h
plot-fft2d-slice-of-vtk.gnuplot

--- data for plotting
output_2D_slice_VTK.txt
output_fft_2D_slice_VTK.txt

