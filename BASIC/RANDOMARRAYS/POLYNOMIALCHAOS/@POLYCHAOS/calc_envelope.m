function env=calc_envelope(PC,u)
% function env=calc_envelope(PC,u)
% Compute the montonone envelope (or monotone majorant) of bounded sequence u

ind=PC.indices;
M=PC.M;
P=length(PC);

if length(u)~=P
    error(['The length ' length(u) ' of the sequence is different from the number of unknown coefficients ' P])
end

env = u;
for i = 1:P
%     ind_sup = [];
%     for j=1:P
%         if all(ind(j,1:M)>=ind(i,1:M))
%             ind_sup = [ind_sup;j];
%         end
%     end
    ind_sup = all(ind(:,1:M)>=repmat(ind(i,1:M),[P 1]),2);
    env(i) = max(abs(u(ind_sup)));
end

end





