function S = cast2matlab_paramelem(file,S)
% function S = cast2matlab_paramelem(file,S)

fid = fopen(file);
if fid==-1
    error(['Le fichier ' file ' n''existe pas']);
end

A = fscanf(fid,'%f %f %f %f %f',[1,5]);
nbnode = A(1);
nbelem = A(2);
node = zeros(nbnode,4);

node(:,[1,2:4]) = fscanf(fid,'%f %f %f %f',[4,nbnode])';

for p=1:nbelem
    temp1 = fscanf(fid,'%f %f',[2,1])';
    typeelemp = fscanf(fid,'%s',1);
    switch typeelemp
        case 'line'
            elem(p).type = 'SEG2';
            elem(p).connec = fscanf(fid,'%f %f ',[2,1])';
        case 'tri'
            elem(p).type = 'TRI3';
            elem(p).connec = fscanf(fid,'%f %f %f',[3,1])';
        case 'quad'
            elem(p).type = 'QUA4';
            elem(p).connec = fscanf(fid,'%f %f %f %f',[4,1])';
        otherwise
            error('type element non reconnu');
    end
end

nbparam = fscanf(fid,'%f ',1);
if ~isempty(nbparam)
    for j=1:nbparam
        B = fscanf(fid,'%f',1)';
    end
    for j=1:nbparam
        S.model.paramname{S.model.nbparam+j} = fscanf(fid,'%s',1);
        if S.model.paramname{S.model.nbparam+j}(end)==','
            S.model.paramname{S.model.nbparam+j} = S.model.paramname{S.model.nbparam+j}(1:end-1);
        else
            fscanf(fid,'%s',1);
        end
    end

    for p=1:nbelem
        temp = fscanf(fid,'%f',[1,1])';
        for j=1:nbparam
            temp = fscanf(fid,'%f',[1,1])';
            elem(p).param{S.model.nbparam+j}{1} = S.model.paramname{S.model.nbparam+j};
            elem(p).param{S.model.nbparam+j}{2} = temp; 
            switch S.model.paramname{S.model.nbparam+j}
                case 'MATT'
                    S.elem(p).mattype=  temp;
            end
        end
    end
end

S.model.nbparam = S.model.nbparam+nbparam;

fclose(fid);

end
