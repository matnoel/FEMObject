function PC = adapt(PC,fun,x,algorithm,bulkparam)
% function PC = adapt(PC,fun,x,algorithm,bulkparam)
% Adapt POLYCHAOS PC by performing an iteration of algorithm
% PC : POLYCHAOS
% fun : function handle
% x : matrix sample
% algorithm : adaptive algorithm for the construction of a nested sequence of finite monotone/lower multi-index sets ('TP' or 'PD' or 'TD' or 'MS' or 'RMS'), 'RMS' by default
%             'TP' or 'PD': isotropic Tensor Product (or Partial Degree) polynomial space
%                    multidimensional space of polynomials of partial degree less or equal to p (in each variable)
%                    update the partial degree by 1 in each dimension at each iteration
%             'TD':  isotropic Total Degree polynomial space
%                    multidimensional space of polynomials of total degree less or equal to p
%                    update the total degree by 1 at each iteration
%             'MS':  Margin Search strategy
%                    add a smallest monotone/lower subset S_n of the margin M_n of a given monotone/lower set A_n
%                    for which energy(S_n)>=bulkparam*energy(M_n), where bulkparam is a bulk parameter
%             'RMS': Reduced Margin Search strategy
%                    add a smallest monotone/lower subset S_n of the reduced margin M_n of a given monotone/lower set A_n
%                    for which energy(S_n)>=bulkparam*energy(M_n), where bulkparam is a bulk parameter
% bulkparam : bulk parameter in (0,1) such that energy(S_n)>=bulkparam*energy(M_n), 0.5 by default
%             bulkparam = 1 selects all multi-indices in the (reduced) margin M_n
%             bulkparam = 0 selects the multi-index in the (reduced) margin M_n
%                           corresponding to the maximal norm of the expansion coefficients

if nargin<4 || isempty(algorithm)
    algorithm = 'RMS';
end
if nargin<5 || isempty(bulkparam)
    bulkparam = 0.5;
end

switch algorithm
    case {'TP','PD','TD'}
        h = RANDPOLYS(PC);
        p = getorder(PC);
        if strcmp(algorithm,'TP') || strcmp(algorithm,'PD')
            PC = POLYCHAOS(h,p+1,'typebase',2);
        elseif strcmp(algorithm,'TD')
            PC = POLYCHAOS(h,p+1,'typebase',1);
        end
    case {'MS','RMS'}
        if strcmp(algorithm,'MS')
            ind_marg = getindices(PC,'margin');
        elseif strcmp(algorithm,'RMS')
            ind_marg = getindices(PC,'reducedmargin');
        end
        PC_tot = addindices(PC,ind_marg,'update');
        A_tot = polyval(PC_tot,x);
        [u_tot,~] = fun(A_tot);
        u_tot = u_tot';
        ind_tot = getindices(PC_tot);
        [~,loc] = ismember(ind_marg,ind_tot,'rows'); % ind_marg = ind_tot(loc,:)
        u_marg = u_tot(:,loc);
        norm_u_marg = sqrt(sum(u_marg.^2,1));
        if strcmp(algorithm,'MS')
            PC_marg = setindices(PC,ind_marg,'update');
            env = calc_envelope(PC_marg,norm_u_marg);
            [~,I] = sort(env,'descend');
        elseif strcmp(algorithm,'RMS')
            [~,I] = sort(norm_u_marg,'descend');
        end
        energy = cumsum(norm_u_marg(I).^2);
        nb_ind = find(energy>=repmat(bulkparam.*energy(end),size(energy)),1);
        PC = addindices(PC,ind_marg(I(1:nb_ind),:),'update');
    otherwise
        error(['Algorithm ' algorithm ' not yet implemeneted'])
end
