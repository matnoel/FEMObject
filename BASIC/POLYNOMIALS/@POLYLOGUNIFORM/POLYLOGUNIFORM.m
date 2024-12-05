function h = POLYLOGUNIFORM(r)
% function h = POLYLOGUNIFORM(r)
% Orthonormal polynomials with respect to the
% measure of a a loguniform random variable
% r : RVLOGUNIFORM
if nargin==0
    r = RVLOGUNIFORM();
end
    h=struct();
    hp = POLYRANDVAR(r);
    h = class(h,'POLYLOGUNIFORM',hp);
