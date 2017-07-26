function saveTHOR(array, timelog, filename,Ox,Oy,Oz,Lx,Ly,Lz,Bx,By,Bz)
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.
    [nx, ny, nz] = size(array);
    fid = fopen(filename, 'wt');
    fprintf(fid, '# iPic3D data from run tred60 in format for the THOR mission\n');
    fprintf(fid, '# Local Magnetic Field [Bx, By, Bz]:  %d   %d   %d\n',Bx,By,Bz);
    fprintf(fid, 'TIME: %d \n',timelog);
    fprintf(fid, 'DATA_SIZE:    %d   %d   %d\n', nx, ny, nz);
    fprintf(fid, 'DATA_FORMAT: DOUBLE\n');
    fprintf(fid, 'DATA_COMPONENTS: 1\n');
    fprintf(fid, 'VARIABLE: phasespace_density\n');
    fprintf(fid, 'BRICK_ORIGIN: %d   %d   %d\n',Ox,Oy,Oz);
    fprintf(fid, 'BRICK_SIZE: %d   %d   %d\n',Lx,Ly,Lz);
    fprintf(fid, '\n');
    for a=1:nz
        for b=1:ny
            for c=1:nx
                fprintf(fid, '%d \n', array(c,b,a));
            end
        end
    end
    fclose(fid);
return
