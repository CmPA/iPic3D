

function savevtkvector(X, Y, Z, filename,label,dx,dy,dz,ox,oy,oz)
%  savevtkvector Save a 3-D vector array in VTK format
%  savevtkvector(X,Y,Z,filename) saves a 3-D vector of any size to
%  filename in VTK format. X, Y and Z should be arrays of the same
%  size, each storing speeds in the a single Cartesian directions.
    if ((size(X) ~= size(Y)) | (size(X) ~= size(Z)))
        fprint('Error: velocity arrays of unequal size\n'); return;
    end
    [nx, ny, nz] = size(X);
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'BINARY\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
    fprintf(fid, 'ORIGIN  %d  %d  %d\n',ox,oy,oz);
    fprintf(fid, 'SPACING %d   %d   %d\n',dx,dy,dz);
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, ['VECTORS ' label ' float\n']);
    V=zeros(3*nx*ny*nz,1);
    V(1:3:end)=reshape(X,nx*ny*nz,1);
    V(2:3:end)=reshape(Y,nx*ny*nz,1);
    V(3:3:end)=reshape(Z,nx*ny*nz,1);
     fwrite(fid,V,'float','b');
%     for a=1:nz
%         for b=1:ny
%             for c=1:nx
%                     fwrite(fid,X(c,b,a),'float','b');
%                     fwrite(fid,Y(c,b,a),'float','b');
%                     fwrite(fid,Z(c,b,a),'float','b');
%             end
%         end
%     end
    fclose(fid);
return
