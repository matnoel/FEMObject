function [ch,paramname]=lirechamp(fid,nbsupp,paramchoix)

nbparam=fscanf(fid,'%f ',1);
for j=1:nbparam
    B=fscanf(fid,'%f',1)';
end
ch=zeros(nbsupp,nbparam); %ch = cell(1,nbparam);
for j=1:nbparam
    name = fscanf(fid,'%s',1);
    if name(end)==','
        name=name(1:end-1);
    else
        fscanf(fid,'%s',1);
    end
    paramname{j}=name;
%ch{j}=zeros(nbsupp,1);
end

for p=1:nbsupp
    temp=fscanf(fid,'%f',[1,1])';
    for j=1:nbparam
        temp=fscanf(fid,'%f',[1,1])';
        ch(p,j)=temp;%ch{j}(p)=temp;
    end    
end


if nargin==3
    ok=0 ;
    for j=1:nbparam
        if strcmp(paramname{j},paramchoix) 
            paramname = paramname{j};
            ch=ch(:,j);
            ok=1;
            break
        end
    end
    if ok==0
        error('le champ n''existe pas');
    end
end
