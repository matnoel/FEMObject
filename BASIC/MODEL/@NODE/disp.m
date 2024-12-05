function disp(u)
% function disp(u)

disp(struct('nbnode',u.nbnode,'syscoord',class(getsyscoord(u))));
table = gettable(u);
if ~isempty(table)
    disp('Table des noeuds')
    disp(gettable(u));
end

if ~isempty(u.lsenrichnode)
    fprintf('[lsenrichnode,lsnumber,lsenrichtype,lsenrichnature] = \n')
    disp([u.lsenrichnode,u.lsnumber,u.lsenrichtype,u.repnature])
    fprintf('lsenrichnature : \n')
    disp([mat2cell([1:length(u.lsenrichnature)]',ones(length(u.lsenrichnature),1),1),u.lsenrichnature(:)]);
end