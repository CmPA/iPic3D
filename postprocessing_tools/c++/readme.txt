parsek-hdf-to-vtk-verX.sh is a bahs script which generates a c++ hdf-to-vtk
code based on the original idea of Stefano Markidis. Ver2 supports HDF5 1.6.5
to 1.6.7. Ver3 supports HDF5 1.8.8 and 1.8.9.

Open the bash file with an editor and modify the inputs. You need to know
what has been written out in the proc****.hdf files. Trying to extract
inexisting quantities may give only a warning, but also crash the script.

The c++ code will automatically include the user's choice of quantities to be
extracted and then the execution is automatic, unless otherwise specified.

It is a long process to extract everything, so on should use a shell screen:

Type in terminal: screen -dR some_screen_name_here

Then execute the bash script: ./parsek-hdf-to-vtk-verX.sh

Exit the screen by holding ctrl and pressing a and d.

View the screen by typing: screen -dR the_screen_name

To see active screens type: screen -list
