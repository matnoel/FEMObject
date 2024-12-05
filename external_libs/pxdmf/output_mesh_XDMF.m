function status = output_mesh_XDMF(filename, divisions, nodes, elements ,nodes_fields, elements_fields, nodes_fields_names, elements_fields_names, varargin)
%
% Function to write a XDMF file:
%
%    output_mesh_XDMF(filename, divisions, nodes, elements ,nodes_fields, elements_fields, nodes_fields_names, elements_fields_names, <options>, ...)
%
%    filename              : file name Note : it will add .xdmf if the filename doesn't have it
%    division              : size 3 vector with the number of divisions. Ex. [2 1 3]
%    nodes                 : matrix (number of nodes, 2 or 3) with the position of the nodes (only nodes used by the elements are printed out)
%                            If z coordinate is absent z =0 is assumed 
%    elements              : connectivity matrix. size(number_of_elements, number_of_nodes_per_element) 
%                           elements supported  3 Triangle, 4 Quadrilateral,  6 Wedges ,  8 Hexahedrons 
%                     NOTE : index start from 0 NOT 1 (see option 'from1')
%    nodes_fields          : a cell with the nodes fields (TimeSteps, fields). 1 column per field, 1 row for each timeStep 
%                           Ex. nodes_fields{1,1} = [1 2 3 4]' (one field (with 4 nodes), and one timestep ). 
%    elements_fields       : a cell with the elements fields (TimeSteps, fields). 1 column per field, 1 row for each timeStep. 
%                           Ex. elements_fields{1,1} = [[1 0 0 ];[.2 .5 .6 ]] (one field (2 elements), and one timestep).
%                            Note : one line for each element. 
%    nodes_fields_names    : a cell with the names of each node field. Ex. nodes_fields_names{1,1} = "Temp"
%    elements_fields_names : a cell with the names of each element field. Ex. elements_fields_names{1,1} = "Damage"
%    options : 'HDF5' for binary output. NOTE : a new file is created filename.h5  (default ASCII)
%              'gzip' option can be used with 'HDF5' compress the information. (default : not compressed)
%              'z=...' force the z coordinate to a value different of 0 (only for 2D meshes). (default z = 0)
%              'from1' the connectivity matrix numbering start from 1 (matlab type) (default : connectivity star form 0)
%              'flipTimeSpace' to write the time inside every domain. (default : the domains are written inside every Timestep) 
%              'rectilinear' to create a rectilinear mesh, so nodes is the min and the max. element = number of elements in each direction 
%    status = 0 if the file was successfully written
%    status = -1 if error
%
%    By Felipe Bordeu - GEM 2010
%    Felipe.Bordeu@ec-nantes.fr
%    
%    Version date : 12/11/2010

try

    opt.HDF5 = 0;
    opt.HDF5_filename= '';
    opt.HDF5_data_counter = 0;
    opt.HDF5_compression = 0;
    opt.z = 0;
    opt.TwoD = 0;
    opt.from1 = 0;
    opt.topology = cell(divisions);
    opt.geometry = cell(divisions);
    opt.current_data = 't';
    opt.flipTimeSpace = 0;
    opt.rectilinear = 0;
    opt.path = '';
    opt.formatdouble = '%12.16f';
%%%%%%%%%%%%%%%%%%%%%
% Path calcul
    offset2 = strfind(filename,'/');
    if ~isempty(offset2)
        opt.path =filename(1:offset2);
    end
%%%%%%%%%%%%%%%%%%%%
% checking HDF5 and compression
%
    for k = 1:length(varargin)
        if strcmpi(varargin{k},'HDF5') 
            if exist('hdf5write')
                opt.HDF5 = 1;
                for j = 1:length(varargin)
                    if strcmpi(varargin{j},'gzip')
                        if  exist('H5ML.id')
                            avail = H5Z.filter_avail('H5Z_FILTER_DEFLATE');
                            if ~avail
                                disp('Warning: gzip filter not available. Using normal output.')
                            else
                                H5Z_FILTER_CONFIG_ENCODE_ENABLED = H5ML.get_constant_value('H5Z_FILTER_CONFIG_ENCODE_ENABLED');
                                H5Z_FILTER_CONFIG_DECODE_ENABLED = H5ML.get_constant_value('H5Z_FILTER_CONFIG_DECODE_ENABLED');
                                filter_info = H5Z.get_filter_info('H5Z_FILTER_DEFLATE');
                                if ( ~bitand(filter_info,H5Z_FILTER_CONFIG_ENCODE_ENABLED) || ~bitand(filter_info,H5Z_FILTER_CONFIG_DECODE_ENABLED) )
                                    disp('Warning: gzip filter not available for encoding and decoding. Using normal output.')
                                    else
                                        opt.HDF5_compression = 1;
                                        disp('Compression available')
                                    end
                                end
                            else
                                disp('Error: Compression not available. Using normal output (HDMF5).')
                            end
                        end
                    end
                else
                    disp('Error: HDF5 not available. Using normal output (ASCII).')
                end
            end
            if strcmpi(varargin{k}, 'from1')
                opt.from1=1;
            end
            if length(varargin{k}) >=2
                if strcmpi(varargin{k}(1:2),'z=')
                    opt.z = str2num(varargin{k}(3:end)) ;
                end 
            end
            if length(varargin{k}) >=length('formatdouble=')
                if strcmpi(varargin{k},'formatdouble=')
                    opt.formatdouble = varargin{k}(14:end) ;
                end     
            end
            if strcmpi(varargin{k}, 'flipTimeSpace')
                opt.flipTimeSpace=1;
            end
            if strcmpi(varargin{k}, 'rectilinear')
                opt.rectilinear = 1;
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forming the filename
%
        if size(filename,2) >= 5
            if strcmpi(filename(end-4:end),'.xdmf') == 0
                opt.HDF5_filename = [filename(length(opt.path)+1:end) '.h5'];
                filename = [filename '.xdmf' ];
            else
                opt.HDF5_filename = [filename(length(opt.path)+1:end-6) '.h5'];
            end
        else
            opt.HDF5_filename = [filename(length(opt.path)+1:end) '.h5'];
            filename = [filename '.xdmf' ];
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% check the input values

        if size(nodes,2) == 2;
            opt.TwoD = 1;
            nodes = [nodes opt.z*ones(size(nodes,1),1) ];
        end

        if opt.TwoD
            divisions = [divisions 1 ];
            elements(3) = 1;
        end

        posmin = min(nodes,[],1);
        posmax = max(nodes,[],1);

        if opt.rectilinear 
%Gen_nodes = zeros(prod(elements+1), 3);
            if size(elements,2) == 2 
                Gen_elements = zeros(prod(elements(1:2)),4);
            else
                Gen_elements = zeros(prod(elements),8);     
            end

            XX = linspace(posmin(1),posmax(1),elements(1)+1); 
            YY = linspace(posmin(2),posmax(2),elements(2)+1);
            ZZ = linspace(posmin(3),posmax(3),elements(3)+1);
            [y,x,z] = meshgrid(YY,XX,ZZ);
            Gen_nodes = [x(:) y(:) z(:)];

%E = 1:elements(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cpt= 1;
            if elements(3) == 0
                for j=0:elements(2)-1
                    for i=0:elements(1)-1
                        Gen_elements(cpt,:) = [1 2  elements(1)+3 elements(1)+2 ]+i+j*(elements(1)+1);
                        cpt =cpt +1;
                    end
                end
            else
                for k =0:elements(3)-1 
                    for j=0:elements(2)-1
                        for i=0:elements(1)-1
                            Gen_elements(cpt,:) = [[1 2  elements(1)+3 elements(1)+2]+(elements(1)+1)*(elements(2)+1)*(k) [1 2  elements(1)+3 elements(1)+2]+(elements(1)+1)*(elements(2)+1)*(k+1) ]+i+j*(elements(1)+1);
                            cpt =cpt +1;
                        end
                    end
                end
            end
            elements = Gen_elements-1;
            nodes = Gen_nodes;
        end
%nnodes = size(nodes,1);
        nelem= size(elements,1);

        if max(size(nodes_fields,2) ~= size(nodes_fields_names,2))
            display(['Error: the number of cells in nodes_fields is ' num2str(size(nodes_fields,2)) '  and the number of cells in nodes_fields_names is ' num2str(size(nodes_fields_names,2))])
            status = -1;
            return;
        end

        if max(size(elements_fields,2) ~= size(elements_fields_names,2))
            display(['Error: the number of cells in elements_fields is ' num2str(size(elements_fields,2)) '  and the number of cells in elements_fields_names is ' num2str(size(elements_fields_names,2))])
            status = -1;
            return;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if opt.from1 == 1
            elements = elements -1;
        end

        file = fopen(filename,'wt');

        if (file<0)
            disp(['Error opening file : ' filename] )
            status = -1;
            return
        end

%to create an empty file
        if opt.HDF5 
            if opt.HDF5_compression == 1
                opt.HDF5_file = H5F.create([opt.path opt.HDF5_filename],'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');
            else
                dset_details.Location = '/info';
                dset_details.Name = 'date';
                hdf5write([opt.path opt.HDF5_filename], dset_details,date );
            end
        end

        print_header(file,filename)

        posmax = posmax + norm(posmax-posmin)/100000;

        X= linspace(posmin(1),posmax(1),divisions(1)+1);
        Y= linspace(posmin(2),posmax(2),divisions(2)+1);
        Z= linspace(posmin(3),posmax(3),divisions(3)+1);

        delta = [X(2) Y(2) Z(2)]-[X(1) Y(1) Z(1)];
        Xzero = [X(1) Y(1) Z(1)];
        elements_to_write = zeros(nelem,1);
        for e=1:nelem
            pos = nodes(elements(e,1)+1,:)-Xzero;
            elements_to_write(e) = sum(floor(pos./delta).*[ 1 (divisions(1)) (divisions(1))*(divisions(2))]);
        end

        maxtime = max(max(size(nodes_fields,1),size(elements_fields,1)),1);

        if opt.flipTimeSpace
            cpt= 0;
            fprintf(file,'    <Grid Name="Domain Space x Time" GridType="Collection" >\n');
            for x=1:divisions(1)
                for y=1:divisions(2)
                    for z=1:divisions(3)
                        elem_to_write = (elements_to_write == cpt);
                        if sum(elem_to_write) >0
                            fprintf(file,'    <Grid Name="Domain %ix%ix%i" GridType="Collection" CollectionType="Temporal" >\n',x,y,z);
                            opt = write_time(file, opt, maxtime);
                            for t=1:maxtime
                                disp(['Writing Domain : ' num2str(x) 'x' num2str(y) 'x' num2str(z) '   TimeStep :' num2str(t)])
                                name = 'Time 1';
                                opt = write_grid(file,[x y z t],nodes,elements, elem_to_write,nodes_fields, elements_fields, nodes_fields_names, elements_fields_names, opt,name);
                            end
                            fprintf(file,'    </Grid>\n');
                        end
                        cpt=cpt+1;
                    end
                end
            end
            fprintf(file,'    </Grid>\n');
        else
            fprintf(file,'    <Grid Name="Domain Space x Time" GridType="Collection" CollectionType="Temporal" >\n');
            opt = write_time(file, opt, maxtime);
            for t=1:maxtime
                fprintf(file,'    <Grid Name="Time %i" GridType="Collection">\n',t);
                cpt= 0;
                for x=1:divisions(1)
                    for y=1:divisions(2)
                        for z=1:divisions(3)
                            elem_to_write = (elements_to_write == cpt);
                            if sum(elem_to_write) >0
                                disp(['Writing Domain : ' num2str(x) 'x' num2str(y) 'x' num2str(z) '   TimeStep :' num2str(t)])
                                name = ['Domain ' num2str(x) 'x' num2str(y) 'x' num2str(z) ];
                                opt = write_grid(file,[x y z t],nodes,elements, elem_to_write,nodes_fields, elements_fields, nodes_fields_names, elements_fields_names, opt,name);
                            end
                            cpt=cpt+1;
                        end
                    end
                end
                fprintf(file,'    </Grid>\n');
            end
            fprintf(file,'    </Grid>\n');
        end

        if opt.HDF5_compression == 1
            H5F.close(opt.HDF5_file);
        end


        print_footer(file)
        fclose(file);

%try
        status = 0;
    catch
        status = -1;
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = write_time(file, opt, maxtime)


fprintf(file,'    <Time TimeType="List">\n');
fprintf(file,'     <DataItem Format="XML" NumberType="Float" Dimensions="%i">\n',maxtime);
fprintf(file,' %f',(1:maxtime)-1);
fprintf(file,'\n     </DataItem>\n');
fprintf(file,'    </Time>\n');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = write_grid(file, partition ,nodes, elements, elements_to_write, nodes_fields, elements_fields, nodes_fields_names, elements_fields_names, opt,name) 

opt.TimeStep = partition(4);
opt.Time = opt.TimeStep;

indice_elements_to_write = find(elements_to_write);

%%%node filtering
nodes_to_write = zeros(size(nodes,1),1);
for e=indice_elements_to_write
for n=1:size(elements,2)
    nodes_to_write(elements(e,n)+1) = 1;
end
end

new_node_numbering =zeros(1,size(nodes,1));
cpt= 0;
for i = 1:size(nodes_to_write,1)
new_node_numbering(i) = cpt;
cpt = cpt +nodes_to_write(i);    
end


fprintf(file,'        <Grid Name="%s" >\n', name);


switch size(elements,2) 
case 8
element_type = 'Hexahedron';
case 6
element_type = 'Wedge';
case 3
element_type = 'Triangle';
case 4
element_type = 'Quadrilateral';    
otherwise
disp('Error Element type not implemented')
return;
end

fprintf(file,'            <Topology TopologyType="%s" NumberOfElements="%d" >\n' ,element_type,sum(elements_to_write));
if opt.TimeStep == 1
opt= print_DataItem(file, indice_elements_to_write,new_node_numbering(elements+1),opt,'%i');
opt.topology{partition(1),partition(2),partition(3)} =  opt.current_data;
else
if opt.HDF5
    fprintf(file,'%s', opt.topology{partition(1),partition(2),partition(3)});
else
    fprintf(file, '<DataItem Reference="XML" Dimensions="%i %i">\n', length(indice_elements_to_write), size(elements,2) );
    if opt.flipTimeSpace
        fprintf(file, '/Xdmf/Domain/Grid/Grid[@Name="Domain %ix%ix%i"]/Grid[@Name="%s"]/Topology/DataItem \n',partition(1),partition(2),partition(3),name);
    else
        fprintf(file, '/Xdmf/Domain/Grid/Grid[@Name="Time 1"]/Grid[@Name="%s"]/Topology/DataItem \n',name);
    end
    fprintf(file,'</DataItem>\n');
end
end

fprintf(file,'            </Topology>\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(file,'			<Geometry GeometryType="XYZ">\n' );
indice_nodes_to_write = find(nodes_to_write);
if opt.TimeStep == 1
opt= print_DataItem(file, indice_nodes_to_write,nodes,opt, opt.formatdouble);
opt.geometry{partition(1),partition(2),partition(3)} =  opt.current_data;
else
if opt.HDF5
    fprintf(file,'%s', opt.geometry{partition(1),partition(2),partition(3)});
else
    fprintf(file, '<DataItem Reference="XML" Dimensions="%i %i">\n', length(indice_elements_to_write), size(elements,2) );
    if opt.flipTimeSpace
        fprintf(file, '/Xdmf/Domain/Grid/Grid[@Name="Domain %ix%ix%i"]/Grid[@Name="%s"]/Geometry/DataItem \n',partition(1),partition(2),partition(3),name);
    else 
        fprintf(file, '/Xdmf/Domain/Grid/Grid[@Name="Time 1"]/Grid[@Name="%s"]/Geometry/DataItem \n',name);
    end
    fprintf(file,'</DataItem>\n');
end
end

fprintf(file,'			</Geometry>\n' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(nodes_fields,2)
if opt.TimeStep <= size(nodes_fields,1) 
    switch size(nodes_fields{opt.TimeStep,i},2)
    case 1
        fprintf(file,'        <Attribute Name="%s" Center="Node" ', nodes_fields_names{i} );
        fprintf(file,' AttributeType="Scalar" ');
        fprintf(file,' >\n');
        opt= print_DataItem(file, indice_nodes_to_write,nodes_fields{opt.TimeStep,i},opt,opt.formatdouble);
        fprintf(file,'        </Attribute>\n' );
    case 3 
        fprintf(file,'        <Attribute Name="%s" Center="Node" ', nodes_fields_names{i} );
        fprintf(file,' AttributeType="Vector" ');
        fprintf(file,' >\n');
        opt= print_DataItem(file, indice_nodes_to_write,nodes_fields{opt.TimeStep,i},opt,opt.formatdouble);
        fprintf(file,'        </Attribute>\n' );
    case 6
        fprintf(file,'        <Attribute Name="%s" Center="Node" ', nodes_fields_names{i} );
        fprintf(file,' AttributeType="Tensor6" ');
        fprintf(file,' >\n');
        opt= print_DataItem(file, indice_nodes_to_write,nodes_fields{opt.TimeStep,i},opt,opt.formatdouble);
        fprintf(file,'        </Attribute>\n' );
    case 9
        fprintf(file,'        <Attribute Name="%s" Center="Node" ', nodes_fields_names{i} );
        fprintf(file,' AttributeType="Tensor" ');
        fprintf(file,' >\n');
        opt= print_DataItem(file, indice_nodes_to_write,nodes_fields{opt.TimeStep,i},opt,opt.formatdouble);
        fprintf(file,'        </Attribute>\n' );
    otherwise
        for j=1:size(nodes_fields{opt.TimeStep,i},2)
            fprintf(file,'        <Attribute Name="%s" Center="Node" ', [nodes_fields_names{i} '_' num2str(j)] );
            fprintf(file,' AttributeType="Scalar" ');
            fprintf(file,' >\n');
            opt = print_DataItem(file, indice_nodes_to_write,nodes_fields{opt.TimeStep,i}(:,j),opt,opt.formatdouble);
            fprintf(file,'        </Attribute>\n' );
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(elements_fields,2)
if opt.TimeStep <= size(elements_fields,1) 
    switch size(elements_fields{opt.TimeStep,i},2)
    case 1
        fprintf(file,'        <Attribute Name="%s" Center="Cell" ', elements_fields_names{i} );
        fprintf(file,' AttributeType="Scalar" ');
        fprintf(file,' >\n');
        opt = print_DataItem(file, indice_elements_to_write,elements_fields{opt.TimeStep,i},opt,opt.formatdouble);
        fprintf(file,'        </Attribute>\n' );   
    case 3
        fprintf(file,'        <Attribute Name="%s" Center="Cell" ', elements_fields_names{i} );
        fprintf(file,' AttributeType="Vector" ');
        fprintf(file,' >\n');
        opt = print_DataItem(file, indice_elements_to_write,elements_fields{opt.TimeStep,i},opt,opt.formatdouble);
        fprintf(file,'        </Attribute>\n' );
    case 6
        fprintf(file,'        <Attribute Name="%s" Center="Cell" ', elements_fields_names{i} );
        fprintf(file,' AttributeType="Tensor6" ');
        fprintf(file,' >\n');
        opt = print_DataItem(file, indice_elements_to_write,elements_fields{opt.TimeStep,i},opt,opt.formatdouble);  
        fprintf(file,'        </Attribute>\n' );

    case 9 
        fprintf(file,'        <Attribute Name="%s" Center="Cell" ', elements_fields_names{i} );
        fprintf(file,' AttributeType="Tensor" ');
        fprintf(file,' >\n');
        opt = print_DataItem(file, indice_elements_to_write,elements_fields{opt.TimeStep,i},opt,opt.formatdouble);  
        fprintf(file,'        </Attribute>\n' );
    otherwise
        for j=1:size(elements_fields{opt.TimeStep,i},2)
            fprintf(file,'        <Attribute Name="%s" Center="Cell" ', [elements_fields_names{i} '_' num2str(j)] );
            fprintf(file,' AttributeType="Scalar" ');
            fprintf(file,' >\n');
            opt = print_DataItem(file, indice_elements_to_write,elements_fields{opt.TimeStep,i}(:,j),opt,opt.formatdouble);
            fprintf(file,'        </Attribute>\n' );
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(file,'        </Grid>\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = print_DataItem(file, indice_to_write, field, opt, format)
data = '';

data = [data sprintf('            <DataItem ')];
if opt.HDF5 == 1 && numel(field) > 100
data = [data sprintf(' Format="HDF"' )];
else
data = [data sprintf(' Format="XML"')];
end
if strcmpi(format,'%i')
data = [data sprintf(' NumberType="int"' )];
else 
data = [data sprintf(' NumberType="Float"' )];
end
data = [data sprintf(' Dimensions="%s">\n', [ num2str(length(indice_to_write)) ' ' num2str(size(field,2)) ] )];

if opt.HDF5 == 1 && numel(field) > 100
if  opt.HDF5_compression ==1
% Create dataspace.  Setting maximum size to [] sets the maximum
% size to be the current size.  Remember to flip the dimensions.
%
DIMS = size(field(indice_to_write,:));
space = H5S.create_simple (2, fliplr(DIMS), []);

%
% Create the dataset creation property list, add the gzip
% compression filter and set the chunk size.  Remember to flip
% the chunksize.
%
dcpl = H5P.create('H5P_DATASET_CREATE');

H5P.set_deflate (dcpl, 9);

%H5P.set_chunk (dcpl, fliplr([1 1]));
H5P.set_chunk (dcpl, fliplr(DIMS));

%DATASET =  ['/dataset' num2str(opt.HDF5_data_counter) '/' 'DataItem' num2str(opt.HDF5_data_counter) ];
DATASET=  ['DataItem' num2str(opt.HDF5_data_counter) ];
if strcmpi(format,'%i')
    H5P.set_chunk (dcpl, fliplr(DIMS));
%
% Create the dataset.
%
    dset = H5D.create(opt.HDF5_file,DATASET,'H5T_STD_I32LE',space,dcpl);
%
% Write the data to the dataset.
%
    H5D.write(dset,'H5T_NATIVE_INT','H5S_ALL','H5S_ALL','H5P_DEFAULT',int32(field(indice_to_write,:))');
else
%
% Create the dataset.
%
    dset= H5D.create(opt.HDF5_file,DATASET,'H5T_IEEE_F64LE',space,dcpl);
%
% Write the data to the dataset.
%
    H5D.write(dset,'H5T_NATIVE_DOUBLE','H5S_ALL','H5S_ALL','H5P_DEFAULT',field(indice_to_write,:)');
end
H5P.close(dcpl);
H5D.close(dset);
H5S.close(space); 
opt.HDF5_data_counter = opt.HDF5_data_counter +1;

data = [data sprintf('%s:%s' , opt.HDF5_filename ,DATASET)];
else
dset_details.Location = ['/dataset' num2str(opt.HDF5_data_counter)];
dset_details.Name = [ 'DataItem' num2str(opt.HDF5_data_counter)  ];

if strcmpi(format,'%i')
    hdf5write([opt.path opt.HDF5_filename], dset_details,int32(field(indice_to_write,:))','WriteMode','append');
else 
    hdf5write([opt.path opt.HDF5_filename], dset_details,field(indice_to_write,:)' ,'WriteMode','append');
end

opt.HDF5_data_counter = opt.HDF5_data_counter +1;

data = [data sprintf('%s:%s/%s' , opt.HDF5_filename ,dset_details.Location, dset_details.Name )];
end
data = [data sprintf('\n            </DataItem>\n' )];
fprintf(file,data);
else 
format_ = '';
for i=1:size(field,2)
format_ = [format_ ' ' format];
end
format_ = [format_ '\n'];

fprintf(file,data);
fprintf(file, format_ , field(indice_to_write,:)' );
fprintf(file, '            </DataItem>\n' );
end
opt.current_data = data;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function print_header(file, filename)

fprintf(file,'<?xml version="1.0" ?>\n');
fprintf(file,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
fprintf(file,'<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude" >\n');
fprintf(file,'    <Domain Name="%s">\n', filename);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function print_footer(file)
fprintf(file,'    </Domain>\n');
fprintf(file,'</Xdmf>\n');
end
