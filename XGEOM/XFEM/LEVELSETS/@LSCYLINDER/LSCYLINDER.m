function ls=LSCYLINDER(dim,varargin)
% function ls=LSCYLINDER(1,Px,Vx,r,normtype)
% function ls=LSCYLINDER(2,Px,Py,Vx,Vy,r,normtype)
% function ls=LSCYLINDER(3,Px,Py,Pz,Vx,Vy,Vz,r,normtype)

ls=struct();

param = cell(8,2);
param(:,1)={'Px','Py','Pz','Vx','Vy','Vz','r','normtype'};
rep = [[1:dim],3+[1:dim],7,8];
param = param(rep,:);

if length(varargin)<=length(rep)-1
    varargin{end+1}=2;
end

param(:,2)=varargin;

lsp = LEVELSET(param,varargin{:});
ls=class(ls,'LSCYLINDER',lsp);
