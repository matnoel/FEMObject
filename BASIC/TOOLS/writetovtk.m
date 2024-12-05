function writetovtk(matrix, filename)
% writetovtk(matrix, filename)
%
% Writes a 3D matrix of doubles as a VTK file. View with paraview.
% The matrix must be 3D.

% Get the matrix dimensions.
[N M O] = size(matrix);

% Open the file.
fid = fopen(filename, 'w');
if fid == -1
    error('Cannot open file for writing.');
    end

% New line.
    nl = sprintf('\n'); % Stupid matlab doesn't interpret \n normally.

% Write the file header.
    fwrite(fid, ['# vtk DataFile Version 2.0' nl 'Volume example' nl 'ASCII' nl ...
        'DATASET STRUCTURED_POINTS' nl 'DIMENSIONS ' ...
        num2str(N) ' ' num2str(M) ' ' num2str(O) nl 'ASPECT_RATIO 1 1 1' nl ...
        'ORIGIN 0 0 0' nl 'POINT_DATA ' ...
        num2str(N*M*O) nl 'SCALARS volume_scalars double 1' nl 'LOOKUP_TABLE default' nl]);

    for z = 1:O
        v = matrix(:, :, z);
        fwrite(fid, num2str(v(:)'));
        fwrite(fid, nl);

% Display progress.
        disp([num2str(round(100*z/O)) '%']);
    end
    fclose(fid);
end
