function pourcentage(m,n,step)
%function pourcentage(m,n,step)

if nargin==2
    step=10;
end
if (mod(m , floor(n/step)) == 0)
    fprintf(' %d%% /', ceil(100 * (m/n)));
    if ceil(100 * (m/n)) == 100
        fprintf('\n');
    end
end
