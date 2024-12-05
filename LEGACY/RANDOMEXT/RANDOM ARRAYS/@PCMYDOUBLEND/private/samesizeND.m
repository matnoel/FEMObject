function [u,v]=samesizeND(u,v)

utemp = u;
vtemp = v;

if isa(u,'PCMYDOUBLEND')
    u=double(u.V);
end
if isa(v,'PCMYDOUBLEND')
    v=double(v.V);
end

if ndims(u)>2 | ndims(v)>2
if isa(utemp,'MYDOUBLEND') 
    if size(u,vtemp.stodim)~=1
     error('mauvaise taille du MYDOUBLEND')
    end
end
if isa(vtemp,'MYDOUBLEND') 
    if size(v,utemp.stodim)~=1
        error('mauvaise taille du MYDOUBLEND')
    end  
end

repunu = find(sizeND(u,1:ndims(v)-2)==1);
repmatu = ones(1,ndims(v)-2);
repmatu(repunu)=sizeND(v,repunu);
if isa(utemp,'PCMYDOUBLEND') 
  repmatu(utemp.stodim-2)=1;
else
  repmatu(vtemp.stodim-2)=1;    
end
    
u = repmat(u,[1,1,repmatu]);

repunv = find(sizeND(v,1:ndims(u)-2)==1);
repmatv = ones(1,ndims(u)-2);
repmatv(repunv)=sizeND(u,repunv);
if isa(vtemp,'PCMYDOUBLEND') 
  repmatv(vtemp.stodim-2)=1;
else
  repmatv(utemp.stodim-2)=1;    
end
v = repmat(v,[1,1,repmatv]);
end


if isa(utemp,'PCMYDOUBLEND')
    utemp.V=MYDOUBLEND(u);
    u=utemp;
else
    u = MYDOUBLEND(u);
end
if isa(vtemp,'PCMYDOUBLEND')
    vtemp.V=MYDOUBLEND(v);
    v=vtemp;
else
    v = MYDOUBLEND(v);
end
