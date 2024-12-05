function status = output_mesh_PXDMF(filename, nodes, elements, names, nodes_fields, cell_fields, nodes_fields_names, cell_fields_names, varargin)
%
% Function to write a PXDMF file:
%
%    output_mesh_PXDMF(filename, nodes, elements, names ,nodes_fields, elements_fields, nodes_fields_names, elements_fields_names, <options>, ...)
%
%    filename              : file name Note : it will add .pxdmf if the filename doesn't have it
%    nodes                 : a cell with the matrix (number of nodes, 2 or 3) containing the positions of the nodes (only nodes used by the elements are printed out)
%                            If z coordinate is absent z =0 is assumed
%                            One cell for each PGD dimentions.
%    elements              : a cell with the connectivity matrix. size(number_of_elements, number_of_nodes_per_element)
%                           elements supported  1nodecell, 2 segment, 3 Triangle, 4 Quadrilateral,  6 Wedges ,  8 Hexahedrons and  mixed meshes (all type suported, but )
%                     NOTE : index start from 0 NOT 1 (see option 'from1')
%    names                 : is a  cell of cell containging the names of
%                            each dimention in each PGD dimention.
%    nodes_fields          : a cell with the nodes fields (PGD_dim, fields). field is a matrix with 1 row for each mode.
%                           Ex. nodes_fields{1,1} = [1 2 3 4; 1 2 5 4] , one field (with 4 nodes and 2 modes), ). Note : one line for each mode .
%    elements_fields       : a cell with the elements fields (PGD_dim, fields). field is a matrix  1 row for each mode.
%                           Ex. elements_fields{1,1} = [[1 0 0 ];[.2 .5 .6 ]] ,one field (3 elements, 2 modes).
%                            Note : one line for each mode.
%    nodes_fields_names    : a cell with the names of each node field. Ex. nodes_fields_names{1,1} = "Temp"
%    elements_fields_names : a cell with the names of each element field. Ex. elements_fields_names{1,1} = "Damage"
%    options : 'HDF5' for binary output. NOTE : a new file is created filename.h5  (default ASCII)
%              'gzip' option can be used with 'HDF5' compress the information. (default : not compressed)
%              'z=...' force the z coordinate to a value different of 0 (only for 2D meshes). (default z = 0)
%              'from1' the connectivity matrix numbering start from 1 (matlab type) (default : connectivity star form 0)
%              'precision=...' to change the precision ('single' or
%              'double')
%    status = 0 if the file was successfully written
%    status = -1 if error
%
%    Note : to used mixed a matrix a extra row must be added with the element type number (not yet implemented, sorry)
%
%
%    By Felipe Bordeu - GEM 2010
%    Felipe.Bordeu@ec-nantes.fr
%
%    Version date : 05/10/2010
%    Version 1.01


opt.HDF5 = 0;
opt.HDF5_filename= '';
opt.HDF5_data_counter = 0;
opt.HDF5_compression = 0;
opt.z = 0;
opt.TwoD = 0;
opt.from1 = 0;
opt.topology = cell(1);
opt.geometry = cell(1);
opt.current_data = 't';
opt.flipTimeSpace = 0;
opt.rectilinear = 0;
opt.path = '';
opt.formatdouble = 'double';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    if length(varargin{k}) >=2
        if strcmpi(varargin{k}(1:2),'z=')
            opt.z = str2num(varargin{k}(3:end)) ;
        end
    end
    if length(varargin{k}) >=length('precision=')
        if strcmpi(varargin{k}(1:10),'precision=')
            opt.formatdouble = varargin{k}(11:end) ;
        end
    end
    if strcmpi(varargin{k}, 'from1')
        opt.from1=1;
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
if size(filename,2) >= 6
    if strcmpi(filename(end-5:end),'.pxdmf') == 0
        opt.HDF5_filename = [filename(length(opt.path)+1:end) '.h5'];
        filename = [filename '.pxdmf' ];
    else
        opt.HDF5_filename = [filename(length(opt.path)+1:end-6) '.h5'];
    end
else
    opt.HDF5_filename = [filename(length(opt.path)+1:end) '.h5'];
    filename = [filename '.pxdmf' ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% check the input values
pgd_dims = size(nodes,1) ;

opt.dim = zeros(pgd_dims,1);

for i = 1:pgd_dims
    if size(nodes{i,1},2) == 2;
        opt.dim = 2;
        nodes{i,1} = [nodes{i,1} opt.z*ones(size(nodes{i,1},1),1) ];
    end
    if size(nodes{i,1},2) == 1;
        opt.dim = 1;
        nodes{i,1} = [nodes{i,1} opt.z*ones(size(nodes{i,1},1),1) opt.z*ones(size(nodes{i,1},1),1)];
    end
end
disp('******************************' )
disp(['File :' filename ])
if(opt.HDF5)
    disp(['Path :' opt.path ])
    disp(['File :' opt.HDF5_filename ])
    if(opt.HDF5_compression )
        disp('Compression : Active')
    end
    disp(['Number Of PGD Dimentions :' num2str(size(nodes,1)) ])
    
    for i = 1:size(nodes,1)
        disp([' PGD Dim :' num2str(i) '  Number of nodes :' num2str(size(nodes{i,1},1),'%5d') '  Number of Elements :' num2str(size(elements{i,1},1)) ]);
    end
    disp('Nodes Fields')
    for i = 1:size(nodes_fields,2)
        disp([ ' Name : ' nodes_fields_names{i} ' Modes : ' num2str(size(nodes_fields{i},1)) ]);
    end
    disp('Cell Fields')
    for i = 1:size(cell_fields,2)
        disp([ ' Name : ' cell_fields_names{i} ' Modes : ' num2str(size(cell_fields{i},1)) ]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % normalization of solution
% for j =1:size(nodes_fields,2)
%     NumberOfModes = size(nodes_fields{1,j},1);
%     alph = ones(NumberOfModes,1);
%     for i = 2:pgd_dims
%         alphN = sqrt(sum(nodes_fields{i,j}.^2,2));
%         alph = alph.*alphN;
%         nodes_fields{i,j} =  bsxfun(@times, nodes_fields{i,j}, 1./alphN);
%
%     end
%     nodes_fields{1,j} =  bsxfun(@times, nodes_fields{1,j}, alph);
% end

%if opt.TwoD
%    divisions = [divisions 1 ];
%    elements(3) = 1;
%end

posmin = cell(pgd_dims,1);
posmax = cell(pgd_dims,1);

nnodes = cell(pgd_dims,1);
nelem= cell(pgd_dims,1);


for i = 1:pgd_dims
    posmin{i} = min(nodes{i},[],1);
    posmax{i} = max(nodes{i},[],1);
    
    nnodes{i} = size(nodes{i},1);
    nelem{i}= size(elements{i},1);
end


%if max(size(nodes_fields,2) ~= size(nodes_fields_names,2))
%    display(['Error: the number of cells in nodes_fields is ' num2str(size(nodes_fields,2)) '  and the number of cells in nodes_fields_names is ' num2str(size(nodes_fields_names,2))])
%    status = -1;
%    return;
%end

%if max(size(elements_fields,2) ~= size(elements_fields_names,2))
%    display(['Error: the number of cells in elements_fields is ' num2str(size(elements_fields,2)) '  and the number of cells in elements_fields_names is ' num2str(size(elements_fields_names,2))])
%    status = -1;
%    return;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt.from1 == 1
    for i = 1:pgd_dims
        elements{i} = elements{i} -1;
    end
end



file = fopen(filename,'wt');

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

if (file<0)
    disp(['Error opening file : ' filename] )
    status = -1;
    return
end

print_header(file,filename)

% posmax = posmax + norm(posmax-posmin)/100000;
%
% X= linspace(posmin(1),posmax(1),divisions(1)+1);
% Y= linspace(posmin(2),posmax(2),divisions(2)+1);
% Z= linspace(posmin(3),posmax(3),divisions(3)+1);
%
% delta = [X(2) Y(2) Z(2)]-[X(1) Y(1) Z(1)];
% Xzero = [X(1) Y(1) Z(1)];
% elements_to_write = zeros(nelem,1);
% for e=1:nelem
%     pos = nodes(elements(e,1)+1,:)-Xzero;
%     elements_to_write(e) = sum(floor(pos./delta).*[ 1 (divisions(1)) (divisions(1))*(divisions(2))]);
% end

% maxtime = max(max(size(nodes_fields,1),size(elements_fields,1)),1);

pgd_dim_dim = zeros(pgd_dims,1);

for i= 1:pgd_dims
    pgd_dim_dim(i) =  length(names{i});
    fprintf(file,'    <Grid Name="PGD%i"  >\n',i);
    fprintf(file,'      <Information Name="Dims" Value="%i"  />\n',pgd_dim_dim(i));
    for j=1:pgd_dim_dim(i)
        fprintf(file,'      <Information Name="Dim%i" Value="%s"  />\n',j-1,names{i}{j});
    end
    opt = write_grid_pgd(file,nodes{i},elements{i},nodes_fields(i,:), cell_fields(i,:), nodes_fields_names, cell_fields_names, opt);
    
end

if opt.HDF5_compression == 1
    H5F.close(opt.HDF5_file);
end


print_footer(file)
fclose(file);


try
    status = 0;
catch
    status = -1;
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = write_grid_pgd(file,nodes,elements,nodes_fields, elements_fields, nodes_fields_names, elements_fields_names, opt)

switch size(elements,2)
    case 8
        element_type = 'Hexahedron';
    case 6
        element_type = 'Wedge';
    case 4
        element_type = 'Quadrilateral';
    case 3
        element_type = 'Triangle';
    case 2
        element_type = 'Polyline';
    case 1
        element_type = 'Polyvertex';
    otherwise
        disp('Error Element type not implemented')
        return;
end

fprintf(file,'            <Topology TopologyType="%s" NumberOfElements="%d" ' ,element_type,size(elements,1));
if strcmp(element_type, 'Polyline'); fprintf(file,' NodesPerElement="2" '); end
if strcmp(element_type, 'Polyvertex'); fprintf(file,' NodesPerElement="1" '); end
fprintf(file,' >\n');

opt= print_DataItem(file, 1:size(elements,1),elements,opt,'%i');


fprintf(file,'            </Topology>\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(file,'			<Geometry GeometryType="XYZ">\n' );

opt= print_DataItem(file, 1:size(nodes,1),nodes,opt, opt.formatdouble);


fprintf(file,'			</Geometry>\n' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:size(nodes_fields,2)
    for mode=1:size(nodes_fields{i},1)
        fprintf(file,'        <Attribute Name="%s_%i" Center="Node" ', nodes_fields_names{i},mode-1 );
        fprintf(file,' AttributeType="Scalar" ');
        fprintf(file,' >\n');
        opt= print_DataItem(file, 1:size(nodes_fields{i}(mode,:),1),nodes_fields{i}(mode,:),opt,opt.formatdouble);
        fprintf(file,'        </Attribute>\n' );
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(elements_fields,2)
    for mode=1:size(elements_fields{i},1)
        fprintf(file,'        <Attribute Name="%s_%i" Center="Cell" ', elements_fields_names{i},mode-1  );
        fprintf(file,' AttributeType="Scalar" ');
        fprintf(file,' >\n');
        opt = print_DataItem(file,  1:size(elements_fields{i}(mode,:),1),elements_fields{i}(mode,:),opt,opt.formatdouble);
        fprintf(file,'        </Attribute>\n' );
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
    data = [data sprintf(' Format="XML"' )];
end

if strcmpi(format,'%i')
    data = [data sprintf(' NumberType="int"' )];
else
    if(strcmpi(format,'single'))
        data = [data sprintf(' NumberType="float" Precision="4"' )];
    else
        data = [data sprintf(' NumberType="float" Precision="8"' )];
    end
end
data = [data sprintf(' Dimensions="%s">\n', [ num2str(length(indice_to_write)) ' ' num2str(size(field,2)) ] )];

if opt.HDF5 == 1 && numel(field) > 100
    if  opt.HDF5_compression ==1 && numel(field) > 1000
        % Create dataspace.  Setting maximum size to [] sets the maximum
        % size to be the current size.  Remember to flip the dimensions.
        %
        DIMS = size(field(indice_to_write,:));
        %space = H5S.create_simple (2, fliplr(DIMS), []);
        space = H5S.create_simple (2, DIMS, []);
        %
        % Create the dataset creation property list, add the gzip
        % compression filter and set the chunk size.  Remember to flip
        % the chunksize.
        %
        dcpl = H5P.create('H5P_DATASET_CREATE');
        
        H5P.set_deflate (dcpl, 9);
        
        %H5P.set_chunk (dcpl, fliplr([1 1]));
        %H5P.set_chunk (dcpl, fliplr(DIMS));
        H5P.set_chunk (dcpl, DIMS);
        
        %DATASET =  ['/dataset' num2str(opt.HDF5_data_counter) '/' 'DataItem' num2str(opt.HDF5_data_counter) ];
        DATASET=  ['DataItem' num2str(opt.HDF5_data_counter) ];
        if strcmpi(format,' %i')
            %H5P.set_chunk (dcpl, fliplr(DIMS));
            %
            % Create the dataset.
            %
            dset = H5D.create(opt.HDF5_file,DATASET,'H5T_STD_I32LE',space,dcpl);
            %
            % Write the data to the dataset.
            %
            H5D.write(dset,'H5T_NATIVE_INT','H5S_ALL','H5S_ALL','H5P_DEFAULT',int32(field(indice_to_write,:))');
        else
            if(strcmpi(format,'single'))
                % Create the dataset.
                dset= H5D.create(opt.HDF5_file,DATASET,'H5T_IEEE_F32LE',space,dcpl);
                % Write the data to the dataset.
                H5D.write(dset,'H5T_NATIVE_FLOAT','H5S_ALL','H5S_ALL','H5P_DEFAULT',single(field(indice_to_write,:)'));
            else
                % Create the dataset.
                dset= H5D.create(opt.HDF5_file,DATASET,'H5T_IEEE_F64LE',space,dcpl);
                % Write the data to the dataset.
                H5D.write(dset,'H5T_NATIVE_DOUBLE','H5S_ALL','H5S_ALL','H5P_DEFAULT',double(field(indice_to_write,:)'));
            end
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
            hdf5write([opt.path opt.HDF5_filename ], dset_details,int32(field(indice_to_write,:))','WriteMode','append');
        else
            if(strcmpi(format,'single'))
                hdf5write([opt.path opt.HDF5_filename ], dset_details, single(field(indice_to_write,:)'),'WriteMode','append');
            else
                hdf5write([opt.path opt.HDF5_filename ], dset_details, double(field(indice_to_write,:)'),'WriteMode','append');
            end
        end
        
        opt.HDF5_data_counter = opt.HDF5_data_counter +1;
        
        data = [data sprintf('%s:%s/%s' , opt.HDF5_filename ,dset_details.Location, dset_details.Name )];
    end
    data = [data sprintf('\n            </DataItem>\n' )];
    fprintf(file,data);
else
    if(strcmpi(format,'single'))
        format = '%12.8f';
    else
        format = '%12.16f';
    end
    
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
