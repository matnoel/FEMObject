function [connecin1in2,connecin1out2,connecout1in2,connecout1out2,xlnodeplus]=...
    bilsdivide_oneelem(ls1,ls2,varargin)

connecin1in2=zeros(0,3);
connecin1out2=zeros(0,3);
connecout1in2=zeros(0,3);
connecout1out2=zeros(0,3);
xlnodeplus = zeros(0,2);

in1 = all(ls1<=0);
in2 = all(ls2<=0);
out1 = all(ls1>=0);
out2 = all(ls2>=0);

if in1
    [connecin1in2,connecin1out2,xlnodeplus] = lsdivide_oneelem(ls2,varargin{:});
elseif out1
    [connecout1in2,connecout1out2,xlnodeplus] = lsdivide_oneelem(ls2,varargin{:});
elseif in2
    [connecin1in2,connecout1in2,xlnodeplus] = lsdivide_oneelem(ls1,varargin{:});
elseif out2
    [connecin1out2,connecout1out2,xlnodeplus] = lsdivide_oneelem(ls1,varargin{:});
else
    
[connecin2,connecout2,xlnodeplus]=lsdivide_oneelem(ls2);

xlnodetotal = [nodelocalcoordtri3();xlnodeplus];

subls1=Ntri3(xlnodetotal)*ls1;

connecin1in2=zeros(0,3);
connecin1out2=zeros(0,3);
connecout1in2=zeros(0,3);
connecout1out2=zeros(0,3);

numnodetotal = size(xlnodetotal,1);

for k=1:size(connecin2,1)

conneck = connecin2(k,:);

subls1k = subls1(conneck);
if all(subls1k>=0)
cout = conneck;cin=zeros(0,3);xlnodeplus=zeros(0,2);
elseif all(subls1k<=0)
cin = conneck;cout=zeros(0,3);xlnodeplus=zeros(0,2);  
else
[cin,cout,xlnodeplus]=lsdivide_oneelem(subls1(conneck),xlnodetotal(conneck,:));
numnodekplus = [conneck,numnodetotal+[1:size(xlnodeplus,1)]];
numnodetotal = numnodetotal+size(xlnodeplus,1);
cout=numnodekplus(cout);
cin=numnodekplus(cin);
end
connecin1in2 = [connecin1in2 ; cin];
connecout1in2 = [connecout1in2 ; cout];
xlnodetotal = [xlnodetotal;xlnodeplus];    
end

for k=1:size(connecout2,1)
conneck = connecout2(k,:);
subls1k = subls1(conneck);
if all(subls1k>=0)
cout = conneck;cin=zeros(0,3);xlnodeplus=zeros(0,2);
elseif all(subls1k<=0)
cin = conneck;cout=zeros(0,3);xlnodeplus=zeros(0,2);  
else
[cin,cout,xlnodeplus]=lsdivide_oneelem(subls1(conneck),xlnodetotal(conneck,:));
numnodekplus = [conneck,numnodetotal+[1:size(xlnodeplus,1)]];
numnodetotal = numnodetotal+size(xlnodeplus,1);
cout=numnodekplus(cout);
cin=numnodekplus(cin);
end
connecin1out2 = [connecin1out2 ; cin];
connecout1out2 = [connecout1out2 ; cout];
xlnodetotal = [xlnodetotal;xlnodeplus];    
end

xlnodeplus = xlnodetotal(4:end,:);

end

if nargin==3
xlnodeplus = Ntri3(xlnodeplus)*varargin{3};    
end



