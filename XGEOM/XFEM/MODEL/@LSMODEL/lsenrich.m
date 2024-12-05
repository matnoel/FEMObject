function S = lsenrich(S,varargin)
% function S = lsenrich(S,varargin)

%S = deleteenrich(S);
S = calc_connec(S);

for i=1:getnblevelsets(S)
   ls = getlevelset(S,i);
   switch getnature(ls)
       case 'crack'
         S=lsenrich_crack(S,ls,varargin{:});  
       case 'material'
         S=lsenrich_material(S,ls,varargin{:});    
   end
end
