How-to install and run iPic3D on your (local) Linux machine.
Doesn't need sudo or root permissions, everything could be installed locally.


1. Install the MPICH2 (for code development) or OpenMPI (for performance runs).

You may choose to use the official "openmpi" library from your operating system's repository.
If you do so, you may skip over to step 2.

For example, this version should work:
http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/1.3.2p1/mpich2-1.3.2p1.tar.gz
During the installation of MPICH follow all the instructions in the README file.

Newer versions of OpenMPI should be fine by default in a linux distro.


2. Install HDF5 library.

You may choose to use the official "libhdf5-openmpi" package from your operating system's repository.
Please note that you also need the "libhdf5-openmpi-dev" package for the required include files.
If you do so, you may skip over to step 3.

The current version working properly with iPIC3D is 1.8.9
http://www.hdfgroup.org/ftp/HDF5/releases/

$ cd <top HDF5 source code directory>
$ ./configure --prefix=<location for HDF5 software>
$ make
$ make check
$ make install


3. Install H5hut library (optional).

If you want to use also the optional H5hut library, you may follow the same steps as for the HDF5 library as described in step 2.
You may download this library from:
https://amas.psi.ch/H5hut/
Build and installation are linked on this website.


4. Set environment variables:

$ PATH=mpich2_install_path/bin/:$PATH
$ PATH=hdf5_install_path/lib/:$PATH
$ LD_LIBRARY_PATH=hdf5_install_path/lib/:$LD_LIBRARY_PATH

Note: these paths can also go at the beginning of the makefile.

where mpich2_install_path is the path to the installed MPICH2, and hdf5_install_path is the path to the installed HDF5.


5. Download the repository:

$ git clone https://github.com/IWF-Graz/iPic3D.git iPic3D_folder

where "iPic3D_folder" is the desired folder for the code.


6. Compiling the code

$ cd iPic3D_folder

Edit the "Makefile" and check the documentation in the "USER-DEFINED SECTION".

Finally, build everything and compile in parallel:

$ make -j


7. Running iPic3D

$ mpiexec -n 4 ./iPic3D inputfiles/GEM.inp
Or for some MPI libraries that do not have "mpiexec":
$ mpirun -n 4 ./iPic3D inputfiles/GEM.inp
