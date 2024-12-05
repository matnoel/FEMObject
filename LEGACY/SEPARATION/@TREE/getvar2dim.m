function var2dim=getvar2dim(T,var,deep)
var2dim=T.var2dim;
if nargin==2
    var2dim=cell2mat(var2dim(var));
elseif nargin==3
    var2dim=var2dim(var);
    var2dim=cellfun(@(v) v(deep),var2dim);
end