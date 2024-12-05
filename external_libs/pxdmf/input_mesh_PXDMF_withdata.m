function data = input_mesh_PXDMF_withdata(filename,varargin)
%
% Function to write a XDMF file:
%
%    data = input_mesh_PXDMF_withdata(filename, <options>, ...)
%
%    filename : file name Note : it will add .pxdmf if the filename doesn't
%               have it
%    data.nodes : a cell with the matrix (number of nodes, 2 or 3) containing
%                 the positions of the nodes (only nodes used by the elements 
%                 are printed out). One cell for each PGD dimentions.
%    data.elements : a cell with the connectivity matrix. 
%                    size(number_of_elements, number_of_nodes_per_element) 
%             NOTE : index start from 0 NOT 1 (see option 'from1')
%    data.names : is a  cell of cell containging the names of each dimention 
%                 in the each PGD dimention.
%    data.nodes_fields : a cell with the nodes fields (PGD_dim, fields). 
%                        field is a matrix with 1 row for each mode.
%     Ex. nodes_fields{1,1} = [1 2 3 4; 1 2 5 4] , one field (with 4 nodes 
%          and 2 modes), ). Note : one line for each mode .
%    data.cell_fields : a cell with the elements fields (PGD_dim, fields). field is a matrix  1 row for each mode. 
%                           Ex. elements_fields{1,1} = [[1 0 0 ];[.2 .5 .6 ]] ,one field (3 elements, 2 modes).
%                            Note : one line for each mode. 
%    data. nodes_fields_names    : a cell with the names of each node field. Ex. nodes_fields_names{1,1} = "Temp"
%    elements_fields_names : a cell with the names of each element field. Ex. elements_fields_names{1,1} = "Damage"
%    options : 'from1' will add 1 to the connectivity matrix (matlab type) (default : nothing is added )
%    status = 0 if the file was successfully written
%    status = -1 if error
%
%    Note : to used mixed a matrix a extra row must be added with the element type number (not yet implemented, sorry)
%
%
%    By Felipe Bordeu - GEM 2010
%    Felipe.Bordeu@ec-nantes.fr
%    
%    Version date : 12/10/2010

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
opt.dtd_file_created = 0;
opt.dtd_file ='Xdmf.dtd';
opt.path = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path calcul
offset2 = strfind(filename,'/');
if ~isempty(offset2)
    opt.path =filename(1:offset2(end));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'from1')
        opt.from1=1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forming the filename
%
if size(filename,2) >= 6
    if strcmpi(filename(end-5:end),'.pxdmf') == 0
        filename = [filename '.pxdmf' ];
    end
else
    filename = [filename '.pxdmf' ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dtd file generation
offset2 = strfind(filename,'/');

if isempty(offset2)
    opt.dt_file ='Xdmf.dtd';
else
    opt.dt_file =[ filename(1:offset2(end)) 'Xdmf.dtd'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating the xdmf.dtd
disp(opt.dt_file)
dtdfile = exist( opt.dt_file,'file');
if dtdfile ~= 2
    opt.dtd_file_created = 1; 
    dtdfile = fopen(opt.dt_file,'w');
    fprintf(dtdfile,'<?xml version="1.0" encoding="UTF-8" ?> 					  \n');
    fprintf(dtdfile,'<!ELEMENT Xdmf (Domain) >									  \n');
    fprintf(dtdfile,'<!ELEMENT Domain (Grid+) >									  \n');
    fprintf(dtdfile,'<!ELEMENT Grid ((Topology,Geometry),Information+,Attribute?) >\n');
    fprintf(dtdfile,'<!ELEMENT Information ANY > 								  \n');
    fprintf(dtdfile,'<!ATTLIST Information Name CDATA #REQUIRED >				  \n');
    fprintf(dtdfile,'<!ATTLIST Information Value CDATA #REQUIRED >				  \n');
    fprintf(dtdfile,'<!ELEMENT Topology (DataItem) > 							  \n');
    fprintf(dtdfile,'<!ELEMENT Geometry (DataItem) > 							  \n');
    fprintf(dtdfile,'<!ELEMENT Attribute (DataItem) >							  \n');
    fprintf(dtdfile,'<!ELEMENT DataItem (#PCDATA) >								  \n');
    fclose(dtdfile);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reading the pxdmf file 
xDoc = xmlread(filename);
xRoot = xDoc.getDocumentElement;
allGrids = xRoot.getElementsByTagName('Grid');
nb_Grids = allGrids.getLength;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seting up the output structure
data.filename = filename;
data.nodes = cell(nb_Grids,1);
data.cells = cell(nb_Grids,1);
data.names = cell(nb_Grids,1);
data.nodes_fields = cell(nb_Grids,0);
data.cell_fields = cell(nb_Grids,0);
data.nodes_fields_names= cell(1,0);
data.cell_fields_names= cell(1,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% straction of name in each dimention
for i = 0:nb_Grids-1
    Grid = allGrids.item(i);
    Infos = Grid.getElementsByTagName('Information');
    nb_dims = GetNbOfDimsInGrid(Infos);
    dim_name = cell(nb_dims,1);
    for j = 0: nb_dims -1 
        for k = 0:Infos.getLength-1
            if strcmpi(Infos.item(k).getAttribute('Name') ,['Dim'  num2str(j)])
                dim_name{j+1} = char(Infos.item(k).getAttribute('Value'));
            end
        end
    end
    data.names{i+1} = dim_name;
    DataItemGeometry = Grid.getElementsByTagName('Geometry').item(0).getElementsByTagName('DataItem').item(0);
    data.nodes{i+1} = Read_data_Item(DataItemGeometry, opt);

    DataItemTopology = Grid.getElementsByTagName('Topology').item(0).getElementsByTagName('DataItem').item(0);
    data.cells{i+1} = Read_data_Item(DataItemTopology, opt) + opt.from1;
    Attributes = Grid.getElementsByTagName('Attribute');
    nb_Attributes = Attributes.getLength;
    for j = 0:nb_Attributes-1
        name = char(Attributes.item(j).getAttribute('Name'));

        type = Attributes.item(j).getAttribute('Center');
        offset1 = strfind(name,'_');
        field_name = strtrim(deblank(name(1:offset1(end)-1)));
        mode = str2double(strtrim(deblank(name(offset1(end)+1:end))));
        DataItemAttributes = Attributes.item(j).getElementsByTagName('DataItem').item(0);
        AttributeData = Read_data_Item(DataItemAttributes, opt);   

        if strcmpi(type,'Node')
            TF = strcmpi(field_name,data.nodes_fields_names);
            index = 0;
            for k =1:size(TF,2)
                if TF(1,k) == 1
                    index = k;
                end
            end
            if index==0
                data.nodes_fields_names{size(data.nodes_fields_names,1),size(data.nodes_fields_names,2)+1} = field_name;
                index = size(data.nodes_fields_names,2);
            end
            if(size(AttributeData,2)~=1) 
            AttributeData =   reshape(AttributeData',[],1);
        end
        data.nodes_fields{i+1,index}(mode+1,:) = AttributeData;

    else
        TF = strcmpi(field_name,data.cell_fields_names);
        index = 0;
        for k =1:size(TF,2)
            if TF(1,k) == 1
                index = k;
            end
        end
        if index==0
            data.cell_fields_names{size(data.cell_fields_names,1),size(data.cell_fields_names,2)+1} = field_name;
            index = size(data.cell_fields_names,2);
        end
        if(size(AttributeData,2)~=1) 
        AttributeData =   reshape(AttributeData',[],1);
    end
    data.cell_fields{i+1,index}(mode+1,:) = AttributeData;
end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elimination of dtd file 
if opt.dtd_file_created == 1; 
delete([opt.path opt.dtd_file])
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = Read_data_Item(DataItem, opt)
format = char(DataItem.getAttribute('Format'));
if strcmpi(format, 'XML')
data = Read_data_Item_XML(DataItem, opt);
else
data = Read_data_Item_H5(DataItem, opt);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = Read_data_Item_XML(DataItem, ~)
%Dimentions = str2num(DataItem.getAttribute('Dimensions'));
%NumberType = char(DataItem.getAttribute('NumberType'));
data = str2num(DataItem.getFirstChild.getData);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = Read_data_Item_H5(DataItem, opt)
Dimentions = str2num(DataItem.getAttribute('Dimensions'));
%NumberType = char(DataItem.getAttribute('NumberType'));
H5_file_and_path = char(DataItem.getFirstChild.getData);
offset1 = strfind(H5_file_and_path,':');
offset2 = strfind(H5_file_and_path,'/');
file = strtrim(deblank(H5_file_and_path(1:offset1-1)));
status = H5F.is_hdf5([opt.path file]);
if status ==0
disp(['File ' file ' is not a h5 file'])
return
else
if status == -1
disp(['File ' file ' corrupted. Sorry'])
return
end
end

if (isempty(offset2))
path ='';
offset2 = offset1;
else
path = strtrim(deblank(H5_file_and_path(offset1+1:offset2(end))));
end
dataset = strtrim(deblank(H5_file_and_path(offset2(end)+1:end)));

%Open the File
H5_file = H5F.open([opt.path file], 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

%Open the Dataset
datasetID = H5D.open(H5_file, [path dataset]);
%Read with no offset
data = H5D.read(datasetID, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT');

if (size(data,1) ~= Dimentions(1))
data =data';
end


H5F.close(H5_file);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nb_dims = GetNbOfDimsInGrid(Infos)
for i=0:Infos.getLength-1
if  strcmpi(Infos.item(i).getAttribute('Name'), 'Dims')
break;
end
end
nb_dims = str2num(Infos.item(i).getAttribute('Value'));

end
