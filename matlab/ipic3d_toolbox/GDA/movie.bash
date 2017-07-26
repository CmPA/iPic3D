#!/bin/bash
echo "inizio mat movie " |mail -s "from buteo on $OF"  valsusa@gmail.com

unset DISPLAY

echo "starting matlab"
#matlab 
matlab > matlab_movie.out  2>&1 << EOF
Ncyc_max=48000
%Ncyc_max=62000
Ncyc_ini=48000
%dir='/data1/gianni/HRmaha3Dopen2/data1/'
dir='/data2/gianni/gda/HRmaha3D3/'
!./clean
%fraciz=1-1./12. %Ygsm=8
%fraciz=.01
fraciz=.6 %Ygsm=4
%fraciz=1./6. %Ygsm=-7old
%fraciz=1./20. %Ygsm=-8odl

Ygsm = 7 

energy_flows
%spectrum_k

%movie_maha3D_slip

%movie_maha3D_rel_drift

%movie_maha3D
%movie_maha3D_vi
%movie_maha3D_p %modified to do only last frame with fiel lines
%movie_maha3D_pi
%!mkdir Y
%!cp *.png Y
%!tar cvf Y.tar Y

%maxp_rate
%maxp
%!tar cvf ELE *.png
%!rm *.png

%maxp_i
%!tar cvf ION *.png
%!rm *.png
%maxp_trace
%!tar cvf TRACE *.png
%!rm *.png

%maxT_agyro
%maxp_energy
%maxp_energy_ion
%maxT_energy_ion
%maxTpar_energy
%maxp_slip_i
%maxp_slip
%!tar cvf SLIP *.png
%!rm *.png
exit

EOF

#convert -density 300 Bn_combo*.eps Bn_combo.gif
#convert  Bn_combo*.png Bn_combo.gif
#convert  Bl_combo*.png Bl_combo.gif
#convert  Bm_combo*.png Bm_combo.gif
#convert  Bn0*.png Bn.gif
#convert  Bl0*.png Bl.gif
#convert  Bm0*.png Bm.gif

#convert Tepar_combo*.png Tepar_combo.gif
#convert Tepar0*.png Tepar.gif
#convert Teper1_combo*.png Teper1_combo.gif
#convert Teper10*.png Teper1.gif
#convert Teper2_combo*.png Teper2_combo.gif
#convert Teper20*.png Teper2.gif

#convert Tipar_combo*.png Tipar_combo.gif
#convert Tipar0*.png Tipar.gif
#convert Tiper1_combo*.png Tiper1_combo.gif
#convert Tiper10*.png Tiper1.gif
#convert Tiper2_combo*.png Tiper2_combo.gif
#convert Tiper20*.png Tiper2.gif

#convert deltay*.png deltaycode.gif
#convert Zgsm*.png Zgsm.gif
#convert Bz_xy*.png Bz_xygsm.gif
#convert Eyrate*.png Eyrate_xygsm.gif


echo "fine mat movie " |mail -s "from buteo on $OF"  valsusa@gmail.com
