function c = FORMSSUM(a,b,varargin)

if nargin==0
    c.forms = {};
    c = class(c,'FORMSSUM');
    superiorto('VARFORM')
elseif nargin==1 && isa(a,'FORMSSUM')
    c = a;
elseif nargin==2 && getn(a)~=getn(b)
    error('on ne peut sommer que des formes du meme type')
elseif nargin==2
    if isa(a,'FORMSSUM') && ~isa(b,'FORMSSUM')
        c=a;
        c.forms = [a.forms,{b}];
    elseif isa(a,'FORMSSUM') && isa(b,'FORMSSUM')
        c=a;
        c.forms = [a.forms,b.forms];
    elseif ~isa(a,'FORMSSUM') && isa(b,'FORMSSUM')
        c=b;
        c.forms = [{a},b.forms];
    elseif ~isa(a,'FORMSSUM') && ~isa(b,'FORMSSUM')
        c.forms = {a,b};
        c = class(c,'FORMSSUM');
        superiorto('VARFORM')
    else
        error('pas les bons arguments')
    end

elseif nargin>2
    c = FORMSSUM(a,b);
    for i=1:nargin-2
        c = FORMSSUM(c,varargin{i});
    end
end

