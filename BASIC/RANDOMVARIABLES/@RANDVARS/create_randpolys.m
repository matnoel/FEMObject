function H = create_randpolys(a,varargin)
% function H = create_randpolys(a,'fedim',i,'femesh',n)
% FE au niveau stochastique pour les dimensions i
% n : nombre de subdivisions de [0,1]
% n : tableau de (length(i)) cellules. n{j} contient le decoupage de [0,1]
% pour la dimension i(j)
%
% function H = create_randpolys(a,'waveletdim',i,'waveletlevel',n)
% Ondelettes au niveau stochastique pour les dimensions i
% n : niveau de resolution (0, 1, ...)
% n : tableau de (length(i)) cellules. n{j} contient le niveau de
% resolution pour la dimension i(j)
%
% function H = create_randpolys(a,'pcgdim',i)
% PC generalise au niveau stochastique pour les dimensions i
% appel de RANDPOLY(a{i}) pour determiner la base polynomiale de la VA a{i}
%
% function H = create_randpolys(a,'pcg')
% PC generalise au niveau stochastique pour toutes les dimensions
% appel de RANDPOLYS(a)
%
% function H = create_randpolys(a,'manudim',i,H)
% RANDPOLYS H pour les dimensions stochastiques i

H = getclassin('RANDPOLY',varargin);
if ~isempty(H)
    H = RANDPOLYS(H);
else
    H = getclassin('RANDPOLYS',varargin);
end


if isempty(H) || length(H)~=a.M
    
    H=RANDPOLYS();
    H(1:a.M)=POLYHERMITE();
    manudim=myunique(getcharin('manudim',varargin));
    pcgdim=myunique(getcharin('pcgdim',varargin));
    fedim=myunique(getcharin('fedim',varargin));
    waveletdim=myunique(getcharin('waveletdim',varargin));
    if ~isempty(manudim)
        [temp,manudim] = ismember(manudim,a);
        H(manudim)=getclassin('RANDPOLYS',varargin);
    end
    if ischarin('pcg',varargin)
        H=RANDPOLYS(a);
    elseif ischarin('pcgdim',varargin)
        [temp,pcgdim] = ismember(pcgdim,a);
        H(pcgdim) = RANDPOLYS(a.RV(pcgdim));
    end
    if ~isempty(fedim)
        [temp,fedim] = ismember(fedim,a);
        n=getcharin('femesh',varargin);
        if ~isa(n,'cell') || length(n)~=length(fedim)
            error('femesh n''est pas correct')
        end
        p=getcharin('order',varargin);
        for k=1:length(fedim)
            H(fedim(k)) = POLYFE(n{k},p(k));
        end
    end
    if ~isempty(waveletdim)
        [temp,waveletdim] = ismember(waveletdim,a);
        n=getcharin('waveletlevel',varargin);
        if ~isa(n,'cell') || length(n)~=length(waveletdim)
            error('waveletlevel n''est pas correct')
        end
        p=getcharin('order',varargin);
        for k=1:length(waveletdim)
            H(waveletdim(k)) = POLYWAVELETS(n{k},p(k));
        end
    end
end

H=setnumber(H,getnumber(a));

function n=myunique(n)

if isa(n,'cell')
    n=[n{:}];
end
n=sort(unique(n));

return