function plotParam(glob,patches,param_out,param_patch,param_in,varargin)
% function plotParam(glob,patches,param_out,param_patch,param_in,varargin)
% Display parameters
% glob: Global
% patches: Patches
% param_out: double or FunctionalTensor
% param_patch: 1-by-numel(patches) cell of double or FunctionalTensor
% param_in: 1-by-numel(patches) cell of double or FunctionalTensor

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

n = numel(patches);

figure('Name','Parameter')
% set(gcf,'Name','Parameter')
clf

subplot(3,n,1:n)
S_out = getfinalmodelpart(glob.S,0);
if israndom(param_out)
    t = mean(param_out);
    p_out = t.data;
else
    p_out = param_out;
end
plot(FENODEFIELD(p_out),removebc(S_out),varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)
for k=1:n
    subplot(3,n,n+k)
    patch = patches.patches{k};
    if israndom(param_patch{k})
        t = mean(param_patch{k});
        p_patch = t.data;
    else
        p_patch = param_patch{k};
    end
    plot(FENODEFIELD(p_patch),patch.S,varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    set(gca,'FontSize',p.Results.FontSize)
    
    subplot(3,n,2*n+k)
    S_in = getfinalmodelpart(glob.S,k);
    if israndom(param_in{k})
        t = mean(param_in{k});
        p_in = t.data;
    else
        p_in = param_in{k};
    end
    if size(p_in,1)==1
        plot(FENODEFIELD(p_in),S_in,varargin{:});
    else
        P_in = calc_P(glob.S,S_in);
        plot(FENODEFIELD(P_in*p_in),S_in,varargin{:});
    end
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    set(gca,'FontSize',p.Results.FontSize)
end

end
