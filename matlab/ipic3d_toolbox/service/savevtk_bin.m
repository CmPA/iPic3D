function savevtk_bin(array, filename,label,dx,dy,dz,Ox,Oy,Oz)
if (nargin<7)
Ox=-dx/2;
Oy=-dy/2;
Oz=-dz/2;
end
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.
    [nx, ny, nz] = size(array);
    fid = fopen(filename, 'w');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'BINARY\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
    fprintf(fid, 'ORIGIN  %d  %d  %d\n',Ox,Oy,Oz);
    fprintf(fid, 'SPACING %d   %d   %d\n',dx,dy,dz);
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, ['SCALARS ' label ' float\n']);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fwrite(fid,reshape(array,1,nx*ny*nz),'float','b');
    fclose(fid)
return

