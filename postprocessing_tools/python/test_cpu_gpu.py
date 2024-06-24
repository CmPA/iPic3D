"""
Created on Tue May 7 12:50 2024

@author: Pranab JD
"""

import os
import glob, h5py
import numpy as np

from datetime import datetime

startTime = datetime.now()

###* =================================================================== *###

dir_ref = "../build_Intel/data/"
dir_data = "./data/"
time_cycle = "cycle_48"

print("Comparing error at time ", time_cycle, "\n")

num_hdf_files = len([name for name in os.listdir(dir_ref) if os.path.isfile(os.path.join(dir_data, name))])

Jx_0_error = np.zeros(num_hdf_files); Jy_0_error = np.zeros(num_hdf_files); Jz_0_error = np.zeros(num_hdf_files); rho_0_error = np.zeros(num_hdf_files)
Jx_1_error = np.zeros(num_hdf_files); Jy_1_error = np.zeros(num_hdf_files); Jz_1_error = np.zeros(num_hdf_files); rho_1_error = np.zeros(num_hdf_files)
Jx_2_error = np.zeros(num_hdf_files); Jy_2_error = np.zeros(num_hdf_files); Jz_2_error = np.zeros(num_hdf_files); rho_2_error = np.zeros(num_hdf_files)
Jx_3_error = np.zeros(num_hdf_files); Jy_3_error = np.zeros(num_hdf_files); Jz_3_error = np.zeros(num_hdf_files); rho_3_error = np.zeros(num_hdf_files)

Ex_error = np.zeros(num_hdf_files); Ey_error = np.zeros(num_hdf_files); Ez_error = np.zeros(num_hdf_files)
Bx_error = np.zeros(num_hdf_files); By_error = np.zeros(num_hdf_files); Bz_error = np.zeros(num_hdf_files)

###? Iterate through all hdf files (one for each processor)
for ii in range(num_hdf_files):

    ###? Reference dataset
    for file in glob.glob(dir_ref + "/proc" + str(ii) + ".hdf"):

        data = h5py.File(file, "r")
        moments = data.get("moments")
        field = data.get("fields")

        ###? Moments (Jx, Jy, Jz, rho)
        sp_0_Jx = moments.get("species_0/Jx"); sp_0_Jy = moments.get("species_0/Jy"); 
        sp_0_Jz = moments.get("species_0/Jz"); sp_0_rho = moments.get("species_0/rho")

        sp_1_Jx = moments.get("species_1/Jx"); sp_1_Jy = moments.get("species_1/Jy"); 
        sp_1_Jz = moments.get("species_1/Jz"); sp_1_rho = moments.get("species_1/rho")

        sp_2_Jx = moments.get("species_2/Jx"); sp_2_Jy = moments.get("species_2/Jy"); 
        sp_2_Jz = moments.get("species_2/Jz"); sp_2_rho = moments.get("species_2/rho")

        sp_3_Jx = moments.get("species_3/Jx"); sp_3_Jy = moments.get("species_3/Jy"); 
        sp_3_Jz = moments.get("species_3/Jz"); sp_3_rho = moments.get("species_3/rho")

        sp_0_Jx_ref = np.array(sp_0_Jx.get(time_cycle)); sp_0_Jy_ref = np.array(sp_0_Jy.get(time_cycle));
        sp_0_Jz_ref = np.array(sp_0_Jz.get(time_cycle)); sp_0_rho_ref = np.array(sp_0_rho.get(time_cycle));

        sp_1_Jx_ref = np.array(sp_1_Jx.get(time_cycle)); sp_1_Jy_ref = np.array(sp_1_Jy.get(time_cycle));
        sp_1_Jz_ref = np.array(sp_1_Jz.get(time_cycle)); sp_1_rho_ref = np.array(sp_1_rho.get(time_cycle));

        sp_2_Jx_ref = np.array(sp_2_Jx.get(time_cycle)); sp_2_Jy_ref = np.array(sp_2_Jy.get(time_cycle));
        sp_2_Jz_ref = np.array(sp_2_Jz.get(time_cycle)); sp_2_rho_ref = np.array(sp_2_rho.get(time_cycle));

        sp_3_Jx_ref = np.array(sp_3_Jx.get(time_cycle)); sp_3_Jy_ref = np.array(sp_3_Jy.get(time_cycle));
        sp_3_Jz_ref = np.array(sp_3_Jz.get(time_cycle)); sp_3_rho_ref = np.array(sp_3_rho.get(time_cycle));
        
        ###? Fields (Ex, Ey, Ez, Bx, By, Bz)
        Bx = field.get("Bx"); By = field.get("By"); Bz = field.get("Bz")
        Ex = field.get("Ex"); Ey = field.get("Ey"); Ez = field.get("Ez")

        Bx_ref = np.array(Bx.get(time_cycle)); By_ref = np.array(By.get(time_cycle)); Bz_ref = np.array(Bz.get(time_cycle))
        Ex_ref = np.array(Ex.get(time_cycle)); Ey_ref = np.array(Ey.get(time_cycle)); Ez_ref = np.array(Ez.get(time_cycle))

    ###? Modified dataset
    for file in glob.glob(dir_data + "/proc" + str(ii) + ".hdf"):

        data = h5py.File(file, "r")
        moments = data.get("moments")
        field = data.get("fields")
        
        ###? Moments (Jx, Jy, Jz, rho)
        sp_0_Jx = moments.get("species_0/Jx"); sp_0_Jy = moments.get("species_0/Jy"); 
        sp_0_Jz = moments.get("species_0/Jz"); sp_0_rho = moments.get("species_0/rho")

        sp_1_Jx = moments.get("species_1/Jx"); sp_1_Jy = moments.get("species_1/Jy"); 
        sp_1_Jz = moments.get("species_1/Jz"); sp_1_rho = moments.get("species_1/rho")

        sp_2_Jx = moments.get("species_2/Jx"); sp_2_Jy = moments.get("species_2/Jy"); 
        sp_2_Jz = moments.get("species_2/Jz"); sp_2_rho = moments.get("species_2/rho")

        sp_3_Jx = moments.get("species_3/Jx"); sp_3_Jy = moments.get("species_3/Jy"); 
        sp_3_Jz = moments.get("species_3/Jz"); sp_3_rho = moments.get("species_3/rho")

        sp_0_Jx_data = np.array(sp_0_Jx.get(time_cycle)); sp_0_Jy_data = np.array(sp_0_Jy.get(time_cycle));
        sp_0_Jz_data = np.array(sp_0_Jz.get(time_cycle)); sp_0_rho_data = np.array(sp_0_rho.get(time_cycle));

        sp_1_Jx_data = np.array(sp_1_Jx.get(time_cycle)); sp_1_Jy_data = np.array(sp_1_Jy.get(time_cycle));
        sp_1_Jz_data = np.array(sp_1_Jz.get(time_cycle)); sp_1_rho_data = np.array(sp_1_rho.get(time_cycle));

        sp_2_Jx_data = np.array(sp_2_Jx.get(time_cycle)); sp_2_Jy_data = np.array(sp_2_Jy.get(time_cycle));
        sp_2_Jz_data = np.array(sp_2_Jz.get(time_cycle)); sp_2_rho_data = np.array(sp_2_rho.get(time_cycle));

        sp_3_Jx_data = np.array(sp_3_Jx.get(time_cycle)); sp_3_Jy_data = np.array(sp_3_Jy.get(time_cycle));
        sp_3_Jz_data = np.array(sp_3_Jz.get(time_cycle)); sp_3_rho_data = np.array(sp_3_rho.get(time_cycle));
        
        ###? Fields (Ex, Ey, Ez, Bx, By, Bz)
        Bx1 = field.get("Bx"); By = field.get("By"); Bz = field.get("Bz")
        Ex = field.get("Ex"); Ey = field.get("Ey"); Ez = field.get("Ez")

        Bx_data = np.array(Bx1.get(time_cycle)); By_data = np.array(By.get(time_cycle)); Bz_data = np.array(Bz.get(time_cycle))
        Ex_data = np.array(Ex.get(time_cycle)); Ey_data = np.array(Ey.get(time_cycle)); Ez_data = np.array(Ez.get(time_cycle))

    ###? Compute error
    Jx_0_diff = np.mean(abs(sp_0_Jx_ref - sp_0_Jx_data))/np.linalg.norm(sp_0_Jx_ref)
    Jy_0_diff = np.mean(abs(sp_0_Jy_ref - sp_0_Jy_data))/np.linalg.norm(sp_0_Jy_ref)
    Jz_0_diff = np.mean(abs(sp_0_Jz_ref - sp_0_Jz_data))/np.linalg.norm(sp_0_Jy_ref)
    rho_0_diff = np.mean(abs(sp_0_rho_ref - sp_0_rho_data))/np.linalg.norm(sp_0_rho_ref)

    Jx_1_diff = np.mean(abs(sp_1_Jx_ref - sp_1_Jx_data))/np.linalg.norm(sp_1_Jx_ref)
    Jy_1_diff = np.mean(abs(sp_1_Jy_ref - sp_1_Jy_data))/np.linalg.norm(sp_1_Jy_ref)
    Jz_1_diff = np.mean(abs(sp_1_Jz_ref - sp_1_Jz_data))/np.linalg.norm(sp_1_Jy_ref)
    rho_1_diff = np.mean(abs(sp_1_rho_ref - sp_1_rho_data))/np.linalg.norm(sp_1_rho_ref)

    Jx_2_diff = np.mean(abs(sp_2_Jx_ref - sp_2_Jx_data))/np.linalg.norm(sp_2_Jx_ref)
    Jy_2_diff = np.mean(abs(sp_2_Jy_ref - sp_2_Jy_data))/np.linalg.norm(sp_2_Jy_ref)
    Jz_2_diff = np.mean(abs(sp_2_Jz_ref - sp_2_Jz_data))/np.linalg.norm(sp_2_Jy_ref)
    rho_2_diff = np.mean(abs(sp_2_rho_ref - sp_2_rho_data))/np.linalg.norm(sp_2_rho_ref)

    Jx_3_diff = np.mean(abs(sp_3_Jx_ref - sp_3_Jx_data))/np.linalg.norm(sp_3_Jx_ref)
    Jy_3_diff = np.mean(abs(sp_3_Jy_ref - sp_3_Jy_data))/np.linalg.norm(sp_3_Jy_ref)
    Jz_3_diff = np.mean(abs(sp_3_Jz_ref - sp_3_Jz_data))/np.linalg.norm(sp_3_Jy_ref)
    rho_3_diff = np.mean(abs(sp_3_rho_ref - sp_3_rho_data))/np.linalg.norm(sp_3_rho_ref)

    Bx_diff = np.mean(abs(Bx_ref - Bx_data))/np.linalg.norm(Bx_ref)
    By_diff = np.mean(abs(By_ref - By_data))/np.linalg.norm(By_ref)
    Bz_diff = np.mean(abs(Bz_ref - Bz_data))/np.linalg.norm(Bz_ref)

    Ex_diff = np.mean(abs(Ex_ref - Ex_data))/np.linalg.norm(Ex_ref)
    Ey_diff = np.mean(abs(Ey_ref - Ey_data))/np.linalg.norm(Ey_ref)
    Ez_diff = np.mean(abs(Ez_ref - Ez_data))/np.linalg.norm(Ez_ref)

    Jx_0_error[ii] = Jx_0_diff; Jy_0_error[ii] = Jy_0_diff; Jz_0_error[ii] = Jz_0_diff; rho_0_error[ii] = rho_0_diff;
    Jx_1_error[ii] = Jx_1_diff; Jy_1_error[ii] = Jy_1_diff; Jz_1_error[ii] = Jz_1_diff; rho_1_error[ii] = rho_1_diff;
    Jx_2_error[ii] = Jx_2_diff; Jy_2_error[ii] = Jy_2_diff; Jz_2_error[ii] = Jz_2_diff; rho_2_error[ii] = rho_2_diff;
    Jx_3_error[ii] = Jx_3_diff; Jy_3_error[ii] = Jy_3_diff; Jz_3_error[ii] = Jz_3_diff; rho_3_error[ii] = rho_3_diff;

    Bx_error[ii] = Bx_diff; By_error[ii] = By_diff; Bz_error[ii] = Bz_diff
    Ex_error[ii] = Ex_diff; Ey_error[ii] = Ey_diff; Ez_error[ii] = Ez_diff


###* =================================================================== *###

print("Error in J(x, y, z) for species 0: ", np.mean(Jx_0_error), ", ", np.mean(Jy_0_error), ", ", np.mean(Jz_0_error))
print("Error in J(x, y, z) for species 1: ", np.mean(Jx_1_error), ", ", np.mean(Jy_1_error), ", ", np.mean(Jz_1_error))
print("Error in J(x, y, z) for species 2: ", np.mean(Jx_2_error), ", ", np.mean(Jy_2_error), ", ", np.mean(Jz_2_error))
print("Error in J(x, y, z) for species 3: ", np.mean(Jx_3_error), ", ", np.mean(Jy_3_error), ", ", np.mean(Jz_3_error))
print()

print("Error in density for species 0, 1, 2, and 3: ", np.mean(rho_0_error), ", ", np.mean(rho_1_error), ", ", np.mean(rho_2_error), ", ", np.mean(rho_3_error))
print()

print("Error in B(x, y, z): ", np.mean(Bx_error), ", ", np.mean(By_error), ", ", np.mean(Bz_error))
print("Error in E(x, y, z): ", np.mean(Ex_error), ", ", np.mean(Ey_error), ", ", np.mean(Ez_error))


print()
print("Complete .....", "Time Elapsed = ", datetime.now() - startTime)