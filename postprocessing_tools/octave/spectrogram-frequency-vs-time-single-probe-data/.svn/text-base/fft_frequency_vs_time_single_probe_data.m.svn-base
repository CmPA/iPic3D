#                           The script uses already extracted data recorded by a single satellite via the postprocessing scripts for virtual probe files
#
#                           Author: Alexander E. Vapirev
#
#                           The FFT contains information between 0 and fs, where fs is the sample frequency.
#                           The sampling frequency is the frequency rate between the adjacent samples.
#                           It is calculated as fs=1/delta_T. In the case of iPIC3D delta_T=wci*dt
#                           So if the frequency axis is plotted as (0:(M-1))/M it is actually f/fs what is plotted over X-axis.
#                           If we multiply by 2Pi/delta_T, e.g. [(0:(M-1))/M]*[2Pi/delta_T] we get it in radians/wci.
#                           If we multiply by 2Pi/dt, e.g. [(0:(M-1))/M]*[2Pi/dt] we get it in radians.


# load and concatenate the extracted data - files of type sat_output_point.txt are matrices (see the postprocessing script for more info)

# eliminate the records after the restart has taken place (after record 10000 in this case) and select column 6 (Ey here) - note that the first index is zero -> means we select the 7th column ---> dataselection1
# then read the data extracted after the restart up the where ever is desired ---> dataselection2
dataselection1=dlmread('/home/alexander/Desktop/data/tred46/output-data-point-x24y7.5z5/sat_output_point.txt',SEP=' ',[1,6,10000,6]);
dataselection2=dlmread('/home/alexander/Desktop/data/tred46.2/output-data-point-x24y7.5z5/sat_output_point.txt',SEP=' ',[1,6,6000,6]);

save dataselection1.txt dataselection1;
save dataselection2.txt dataselection2;

# concatenate them vertically - horizontal concatenation is B=[A1,A2]
dataselection=[dataselection1;dataselection2];

A=dataselection;

# get the number of data points in A - size reports number of rows and cols so we take the bigger, i.e., rows using max
Npoints = max(size(A))

# number of FFT points - change the number if desired by changing the denominator
NFFT = floor(Npoints/1)

# next in iPic3D code units
b0 = 0.0097
dt = 0.125
wci = b0
delta_t = wci*dt

# sampling frequency - fs=1/delta_t [delta_t in code units = wci*dt]
Fs = 2.0*3.1415927/((NFFT/2+1)*delta_t)

# plot the spectrogram and then test the result

 x = A;
 step=ceil(20*Fs/1);          # one spectral slice every XXX Fs
 window=ceil(100*Fs/1); # 100 Fs data window

 ## test of automatic plot
 [S, f, t] = specgram(x);
 specgram(x, 2^nextpow2(window), Fs, window, window-step);
 disp("shows a diagonal from bottom left to top right");
 input("press enter:","s");

 ## test of returned values
 S = specgram(x, 2^nextpow2(window), Fs, window, window-step);
 imagesc(20*log10(flipud(abs(S))));
 disp("same again, but this time using returned value");
