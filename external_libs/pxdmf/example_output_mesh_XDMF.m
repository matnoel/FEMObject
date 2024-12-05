%    By Felipe Bordeu - GEM 2010
%    Felipe.Bordeu@ec-nantes.fr
%    
%    Version date : 29/09/2010

clear all
close all
divisions = [1 1 3];

%the nodes
nodes = [0.0    0.0    0.0
1.0    0.0    0.0
0.0    1.0    0.0
1.0    1.0    0.0
0.0    0.0    1.0
1.0    0.0    1.0
0.0    1.0    1.0
1.0    1.0    1.0
0.0    0.0    2.0
1.0    0.0    2.0
0.0    1.0    2.0
1.0    1.0    2.0
];

% wedges 
elements = [0 1 4 2 3 6
4 5 9 6 7 11];

% dep
nodes_fields = cell(1,1);
nodes_fields{1,1} = [    0.9501    0.9218    0.1389
0.2311    0.7382    0.2028
0.6068    0.1763    0.1987
0.4860    0.4057    0.6038
0.8913    0.9355    0.2722
0.7621    0.9169    0.1988
0.4565    0.4103    0.0153
0.0185    0.8936    0.7468
0.8214    0.0579    0.4451
0.4447    0.3529    0.9318
0.6154    0.8132    0.4660
0.7919    0.0099    0.4186
];
for i=1:1
    nodes_fields{i,1} = rand(12,3);
end
% temp
nodes_fields{1,2} = [0 1 2 3 4 5 6 4 5 6 7 8 9 10 11]';
nodes_fields{2,2} = [0 2 4 6 8 10 12 8 10 12 14 16 18 20 22]';

nodes_fields_names= cell(1,2);
nodes_fields_names{1}= 'dep';
nodes_fields_names{2}= 'temp';

mesh_fields = cell(1,1) ;
mesh_fields{1,1} = [0 1]';
mesh_fields{2,1} = [2 3]';
mesh_fields_names= cell(1,1);
mesh_fields_names{1}= 'P';


a = zeros(1,7);
%binary and compressed
filename= 'tests/binary_comp.xdmf';
a(1)=output_mesh_XDMF(filename, divisions, nodes, elements ,nodes_fields, mesh_fields, nodes_fields_names, mesh_fields_names,'HDF5','gzip');

%binnary
filename= 'tests/binary.xdmf';
a(2)=output_mesh_XDMF(filename, divisions, nodes, elements ,nodes_fields, mesh_fields, nodes_fields_names, mesh_fields_names,'HDF5');

%ascii
filename= 'tests/ASCII.xdmf';
a(3)=output_mesh_XDMF(filename, divisions, nodes, elements ,nodes_fields, mesh_fields, nodes_fields_names, mesh_fields_names);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%binary and compressed
filename= 'tests/binary_comp_flip.xdmf';
a(4)=output_mesh_XDMF(filename, divisions, nodes, elements ,nodes_fields, mesh_fields, nodes_fields_names, mesh_fields_names,'HDF5','gzip','flipTimeSpace');

%binnary
filename= 'tests/binary_flip.xdmf';
a(5)=output_mesh_XDMF(filename, divisions, nodes, elements ,nodes_fields, mesh_fields, nodes_fields_names, mesh_fields_names,'HDF5','flipTimeSpace');

%ascii
filename= 'tests/ASCII_flip.xdmf';
a(6)=output_mesh_XDMF(filename, divisions, nodes, elements ,nodes_fields, mesh_fields, nodes_fields_names, mesh_fields_names,'flipTimeSpace');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'tests/ASCII_rectilinear';
a(7)=output_mesh_XDMF('test',[1 1 1], [[0 0 1];[1 2 3]], [2 3 6], cell(0,0),cell(0,0),cell(0,0),cell(0,0),'rectilinear');

if sum(a==0) == size(a,2)
    disp('run without error')
else
    disp('run with error')
end
