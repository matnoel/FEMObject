clc
%%clear all
close all

[SUCCESS,MESSAGE,MESSAGEID] = mkdir('tests');

%nodes is a cell contaning the nodes in each dimention
%the nodes
nodes = cell(3,1);
nodes{1} = [0.0    0.0    0.0
1.0    0.0    0.0
2.0    0.0    0.0
3.0    0.0    0.0
4.0    0.0    0.0
5.0    0.0    0.0
6.0    0.0    0.0
7.0    0.0    0.0
8.0    0.0    0.0
9.0    0.0    0.0
10.0    0.0    0.0
11.0    0.0    0.0];

nodes{2} = [0.0    0.0    0.0
1.0    0.0    0.0
2.0    0.0    0.0
3.0    0.0    0.0
4.0    0.0    0.0
5.0    0.0    0.0
6.0    0.0    0.0
7.0    0.0    0.0
8.0    0.0    0.0
9.0    0.0    0.0
10.0    0.0    0.0
11.0    0.0    0.0
];


nodes{3} = [0.0    0.0    0.0
1.0    0.0    0.0
2.0    0.0    0.0
3.0    0.0    0.0
4.0    0.0    0.0
5.0    0.0    0.0
6.0    0.0    0.0
7.0    0.0    0.0
8.0    0.0    0.0
9.0    0.0    0.0
10.0    0.0    0.0
11.0    0.0    0.0
];

%elements is a cell contaning the elements in each dimention
cells = cell(3,1);
cells{1} = [1 2
2 3
3 4
4 5 
5 6 
6 7
7 8
8 9
9 10
10 11
11 12];

cells{2} = [1 2
2 3
3 4
4 5 
5 6 
6 7
7 8
8 9
9 10
10 11
11 12];

cells{3} = [1 
2
3
4
5 
6 
7
8
9
10
11
12];

%names is a cell contaning the name of each dimention (the number of name determine the size of the dimention)
names = cell(3,1);
name = cell(1);
name{1} = 'x';
names{1} = name;
name{1} = 'y';
names{2} = name;
name{1} = 'z';
names{3} = name;

% three fields (dep_x, dep_y, dep_z)
nodes_fields = cell(3,3);
nodes_fields{1,1} =   [ -1   0.2311   0.6068   0.4860   0.8913   0.7621   0.4565   0.0185   0.8214   0.4447 0.6154   0.7919   
0   0.2311   0.6068   0.4860   0.8913   0.7621   0.4565   0.0185   0.8214   0.4447 0.6154   0.7919  ];
nodes_fields{2,1} =   [ 1 1 1 1 1 1 1 1 1 1 1 1 
2 1 1 4 1 7 1 1 0 1 -1 1];
nodes_fields{3,1} =   [ 3 1 1 1 1 1 1 1 1 1 1 1 
4 1 1 1 1 1 1 1 1 1 1 1];	

nodes_fields{1,2} =   [ 5 0.7382 0.1763 0.4057 0.9355 0.9169 0.4103 0.8936 0.0579 0.3529 0.8132  0.0099  ];	
nodes_fields{2,2} =   [ 6 1 1 1 1 1 1 1 1 1 1 1 ];
nodes_fields{3,2} =   [ 7 1 1 1 1 1 1 1 1 1 1 1 ];

for i=1:3
    nodes_fields{i,3} = rand(1,12);
end

% so three fields names (dep_x, dep_y, dep_z)
nodes_fields_names= cell(1,2);
nodes_fields_names{1}= 'dep_x';
nodes_fields_names{2}= 'dep_y';
nodes_fields_names{3}= 'dep_z';

% like the nodes_fields 
cell_fields = cell(3,1) ;
cell_fields{1,1} = [0 1 2 3 4 5 6 7 8 9 10 ];
cell_fields{2,1} = [1 2 3 4 5 6 7 8 9 10 11 ];
cell_fields{3,1} = [2 3 4 5 6 7 8 9 10 11 12 13];
cell_fields_names= cell(1,1);
cell_fields_names{1}= 'P';


%binary and compressed

filename= 'tests/binary_comp.pxdmf';
output_mesh_PXDMF(filename, nodes, cells, names, nodes_fields, cell_fields, nodes_fields_names, cell_fields_names, 'HDF5', 'gzip', 'from1');

filename= 'tests/binary.pxdmf';
output_mesh_PXDMF(filename, nodes, cells, names, nodes_fields, cell_fields, nodes_fields_names, cell_fields_names, 'HDF5', 'from1');

filename= 'tests/Ascii.pxdmf';
output_mesh_PXDMF(filename, nodes, cells, names, nodes_fields, cell_fields, nodes_fields_names, cell_fields_names, 'from1');

%also work for a struct with the data
data.filename = 'tests/Ascii_struct';
data.nodes = nodes;
data.cells = cells;
data.names = names;
data.nodes_fields = nodes_fields;
data.cell_fields = cell_fields;
data.nodes_fields_names =nodes_fields_names;
data.cell_fields_names = cell_fields_names;
data.flags = cell(1);
data.flags{1} = 'from1';

output_mesh_PXDMF_withdata(data);

data_Ascii = input_mesh_PXDMF_withdata('tests/Ascii.pxdmf', 'From1');
data_bin = input_mesh_PXDMF_withdata('tests/binary.pxdmf', 'From1');
data_compress = input_mesh_PXDMF_withdata('tests/binary_comp', 'From1');

% to check the integrity of the function 

result = zeros(21,1);

result(1)= min(min(cell2mat(data.nodes)	         == cell2mat(data_Ascii.nodes)));
result(2)= min(min(cell2mat(data.cells(1:2)) 	     == cell2mat(data_Ascii.cells(1:2))));
result(3)= min(min(cell2mat(data.cells(3:3)) 	     == cell2mat(data_Ascii.cells(3:3))));
result(4)= min(min(cell2mat(data.nodes_fields(1))   == cell2mat(data_Ascii.nodes_fields(1))));
result(5)= min(min(cell2mat(data.nodes_fields(2:3)) == cell2mat(data_Ascii.nodes_fields(2:3))));
result(6)= min(min(cell2mat(data.cell_fields(1:2))  == cell2mat(data_Ascii.cell_fields(1:2))));
result(7)= min(min(cell2mat(data.cell_fields(3:3))  == cell2mat(data_Ascii.cell_fields(3:3))));

result(8)= min(min(cell2mat(data.nodes)	         == cell2mat(data_bin.nodes)));
result(9)= min(min(cell2mat(data.cells(1:2)) 	     == cell2mat(data_bin.cells(1:2))));
result(10)= min(min(cell2mat(data.cells(3:3)) 	     == cell2mat(data_bin.cells(3:3))));
result(11)= min(min(cell2mat(data.nodes_fields(1))   == cell2mat(data_bin.nodes_fields(1))));
result(12)= min(min(cell2mat(data.nodes_fields(2:3)) == cell2mat(data_bin.nodes_fields(2:3))));      
result(13)= min(min(cell2mat(data.cell_fields(1:2))  == cell2mat(data_bin.cell_fields(1:2))));  
result(14)= min(min(cell2mat(data.cell_fields(3:3))  == cell2mat(data_bin.cell_fields(3:3))));  

result(15)= min(min(cell2mat(data.nodes)	         == cell2mat(data_compress.nodes)));
result(16)= min(min(cell2mat(data.cells(1:2)) 	     == cell2mat(data_compress.cells(1:2)))); 			     
result(17)= min(min(cell2mat(data.cells(3:3)) 	     == cell2mat(data_compress.cells(3:3))));		     
result(18)= min(min(cell2mat(data.nodes_fields(1))   == cell2mat(data_compress.nodes_fields(1))));   
result(19)= min(min(cell2mat(data.nodes_fields(2:3)) == cell2mat(data_compress.nodes_fields(2:3))));      
result(20)= min(min(cell2mat(data.cell_fields(1:2))  == cell2mat(data_compress.cell_fields(1:2)))); 
result(21)= min(min(cell2mat(data.cell_fields(3:3))  == cell2mat(data_compress.cell_fields(3:3))));

if (min(result) == 0 )
    disp('Error in the I/O funtions')
else
    disp('Everything looks fine.')
end

