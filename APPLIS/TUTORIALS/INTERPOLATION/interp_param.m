clear all

n = 100;
X = linspace(0,2*pi,n)';
uf = @(x,mu) cos(x+mu).*sin(x*mu);

% Definition of the training set and computation of the snapshots
m_snap = 100;
Xi = 2*pi*rand(1,m_snap);
snapshots = zeros(n,m_snap);
for i = 1:m_snap
    snapshots(:,i) = uf(X,Xi(i));
end

figure(1)
clf
plot(snapshots)


%% Computation of the reference solution on a testing set
m_test = 1000;
Xi_test = 2*pi*rand(1,m_test);
ref_test = zeros(n,m_test);
for i = 1:m_test
    ref_test(:,i) = uf(X,Xi_test(i));
end


error_2 = @(x) norm(ref_test-x,'fro')/norm(ref_test,'fro');
error_inf = @(x) max(max(abs(ref_test-x)))/max(max(abs(ref_test)));

%% SHEPARD
Jk = mat2cell(snapshots,n,ones(1,m_snap));
p = 3.5; % smoothness coefficient in Shepard interpolation
shep = @(mu) shepard_interp(Jk,Xi,mu,p);

shepard_test = zeros(n,m_test);

for i = 1:m_test
    shepard_test(:,i) = shep(Xi_test(i));
end

shepard_2 = error_2(shepard_test);
shepard_inf = error_inf(shepard_test);

%% EIM
clear opts
opts.max_basis = 40;
opts.tol = 1e-14;

[approx,Q,indxi,indx] = eim(snapshots,opts);
norm(approx-snapshots)
Qred = Q(indx,:);

eim_test = zeros(numel(indx),m_test);
for i = 1:m_test
    eim_test(:,i) = uf(X(indx),Xi_test(i));
end
eim_test = Q*(Qred\eim_test);

eim_2 = error_2(eim_test);
eim_inf = error_inf(eim_test);

%% EIM with projection
opts.orth = true;
opts.dot = speye(n); % The approximation is now defined with a projection
[approx2,Q2,indxi2] = eim(snapshots,opts);

eim_proj_test = Q2*((Q2'*opts.dot*Q2) \ (Q2'*ref_test));

eim_proj_2 = error_2(eim_proj_test);
eim_proj_inf = error_inf(eim_proj_test);

%% POD of the snapshots + EIM approximation of the basis
tol = 1e-14;
[Qdeim,inddeim] = pod_deim(snapshots,tol);

deim_test = zeros(numel(inddeim),m_test);
for i = 1:m_test
    deim_test(:,i) = uf(X(inddeim),Xi_test(i));
end
Qdred = Qdeim(inddeim,:);
deim_test = Qdeim*(Qdred\deim_test);

deim_2 = error_2(deim_test);
deim_inf = error_inf(deim_test);

%%
fprintf('\nshepard_2 = %d\n',shepard_2)
fprintf('eim_2 = %d\n',eim_2)
fprintf('eim_proj_2 = %d\n',eim_proj_2)
fprintf('deim_2 = %d\n\n',deim_2)

fprintf('shepard_inf = %d\n',shepard_inf)
fprintf('eim_inf = %d\n',eim_inf)
fprintf('eim_proj_inf = %d\n',eim_proj_inf)
fprintf('deim_inf = %d\n',deim_inf)
