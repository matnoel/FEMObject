function x = PCCELL(varargin)
% function x = PCCELL(a,PC)
% a : PCARRAY ou PCRADIAL ou cell

if isa(varargin{1},'PCCELL')
    x = varargin{1};
elseif isa(varargin{1},'PCARRAY')
    s = size(varargin{1});
    PC = getPC(varargin{1});
    x.value = squeeze(num2cell(double(varargin{1}),1:length(s)-1));
    x=class(x,'PCCELL',PC);

else    
    PC = getclassinvarargin('POLYCHAOS',varargin);
    C = getclassinvarargin('cell',varargin);

    if getP(PC)+1~=length(C)
        error('le nombre de cellule ne correspond pas a la dimension du chaos')
    end
    x.value = C ;
    x=class(x,'PCCELL',getPC(PC));
end

