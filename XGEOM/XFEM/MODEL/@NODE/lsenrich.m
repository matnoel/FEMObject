function u = lsenrich(u,numnode,ls,choix,k)
% function u = lsenrich(u,numnode,ls,choix,k)
% numnode : liste de noeuds à enrichir
% ls : levelset (de nature material ou crack)
% attention : ne prevoit qu'un type d'enrichissement par noeud
% on encrase l'ancien enrichissement du noeud si doublon
% choix : 'tip' ou 'support' info complementaire pour 
% k : dans le cas ou enrichissement tip : precise le numero du tip
if ~isempty(numnode) & isenrich(ls)    

rep = ismember(u.lsenrichnode,numnode);
u.lsenrichnode(rep)=[];
u.repnature(rep)=[];
u.lsenrichtype(rep)=[];
u.lsnumber(rep)=[];

switch getnature(ls)
    case 'crack'     
        if nargin>=4 && strcmpi(choix,'tip')
 u.lsenrichnature = [u.lsenrichnature , {['tip' num2str(k)]}];     
 enrichtype = getenrichtypetip(ls,k);
        elseif nargin>=4 && strcmpi(choix,'support')
 u.lsenrichnature = [u.lsenrichnature , {'support'}];      
 enrichtype = getenrichtypesupport(ls);
        else
            error('preciser choix ''tip'' ou ''support''')
        end
    case 'material'
        if nargin<4
            choix='material';
        end
u.lsenrichnature = [u.lsenrichnature , {choix}]; 
enrichtype = getenrichtype(ls);
    otherwise
        error('type de levelset non definie pour l''enrichissement')
end

u.lsenrichnode = [u.lsenrichnode;numnode(:)];
u.lsenrichtype = [u.lsenrichtype;repmat(enrichtype,numel(numnode),1)];
u.lsnumber = [u.lsnumber;repmat(getnumber(ls),numel(numnode),1)];
u.repnature = [u.repnature;repmat(length(u.lsenrichnature),numel(numnode),1)];



end
