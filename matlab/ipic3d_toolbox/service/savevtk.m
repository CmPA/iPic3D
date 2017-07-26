function savevtk(array, filename,label,dx,dy,dz)
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.
    [nx, ny, nz] = size(array);
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
    fprintf(fid, '\n');
    fprintf(fid, 'ORIGIN  %d  %d  %d\n',-dx/2,-dy/2,-dz/2);
    fprintf(fid, 'SPACING %d   %d   %d\n',dx,dy,dz);
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, ['SCALARS ' label ' double\n']);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '\n');
    for a=1:nz
        for b=1:ny
            for c=1:nx
                fprintf(fid, '%d ', array(c,b,a));
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
return

Vector data:
