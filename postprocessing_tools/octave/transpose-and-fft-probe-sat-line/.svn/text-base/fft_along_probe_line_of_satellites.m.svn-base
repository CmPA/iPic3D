% Author: Alexnader E. Vapirev

% script uses already extracted data along a line of satellites via the postprocessing scripts for virtual probe files

% load and concatenate the matrices

% matrices in these files are time along x and satellites/coordinate along y - they are horizontal matrices for gnuplot splot format

matrix1=load("/home/alexander/Desktop/data/tred46/output-data-block-single-var-y7.5z5/output_filename_Ey.txt");
matrix2=load("/home/alexander/Desktop/data/tred46.2/output-data-block-single-var-y7.5z5/output_filename_Ey.txt");

F1=matrix1(:,1:10000); % eliminate the records after the restart has taken place
F2=matrix2;

% concatenate them horizontally - vertical concatenation is B=[A1;A2]
F=[F1,F2];

%transpose to vertical matrix
A=F';

% perform FFT along the space dimension

dimension1 = min(size(A))  % space - number of cols/satellite probes
dimension2 = max(size(A))  % time  - number of rows/time steps

NFFT = dimension1;

%w = hamming(NFFT);
%avg_A = A;
%for i = 2:dimension1
%       avg_A(i,1:dimension2)=A(i,1:dimension2)*.01+avg_A(i-1,1:dimension2)*.99;
%end

B=[];

for j = 1:dimension2
       row_A = A(j,1:dimension1);                       % get one row of A
       FFT_row_A = fft2(row_A);                         % FFT that row
       abs_FFT_row_A = 2*abs(FFT_row_A(1:NFFT/2+1));    % get absolute of the FFT result
       B = [B;abs_FFT_row_A];                           % concatenate the resulted rows in a 2D matrix
end

B_transposed=B';

% output the file in a horizontal matrix again - the first 5 lines contain brief description of the data and then the matrix itself
save fft_B.txt B_transposed;

%pcolor(log10(B_transposed));
%colorbar('log10');
