function SM=init_SEPMODEL(SM)

F={};
for d=1:SM.dim
    model=SM.F{d};
    
    if isa(model,'MODEL')
        % Modele deterministe
        Fnew=struct();
        Fnew.model  = model;
        Fnew.type   = 'SPACE';
        Fnew.masse  = getmasse(model);
        Fnew.metric = getmetric(model);
        Fnew.nbddl  = getnbddl(model);
        
        F=[F {Fnew}];
        
    elseif isa(model,'TIMEMODEL')
        % Modele temporel
    	Fnew=struct();
        Fnew.model  = model;
        Fnew.type   = 'TIME';
        Fnew.masse  = getMtrimatrix(model);
        Fnew.metric = getMmatrix(model);
        Fnew.nbddl  = getnbtimedof(model);
        
        F=[F {Fnew}];
        
    elseif isa(model,'POLYCHAOS')
        % Modele stochastique
    	Fnew=struct();
        
        if ~iscalcmasse(model),  model=calc_masse(model);  end
        if ~iscalcmetric(model), model=calc_metric(model); end
        Fnew.model  = model;
        Fnew.type   = 'STOCH';
        Fnew.masse  = getmasse(model);
%         Fnew.metric = getmetric(model);
        Fnew.metric = getmetric(calc_metric(model));
        Fnew.nbddl  = getP(model)+1;
        
        F=[F {Fnew}];
        
    elseif isa(model,'POLYCHAOSTP')
        % Modele stochastique sous format
        % PCTPMODEL (equivalent de SEPMODEL
        % pour le stochastique)
        
        subdim=getnbdim(model);
        Fnew=cell(1,subdim);
        for sd=1:subdim
            smodel=getpcgroup(model,sd);
            if ~iscalcmasse(smodel),  smodel=calc_masse(smodel);  end
            if ~iscalcmetric(smodel), smodel=calc_metric(smodel); end
            Fnew{sd}.model  = smodel;
            Fnew{sd}.type   = 'STOCH';
            Fnew{sd}.masse  = getmasse(smodel);
            Fnew{sd}.metric = getmetric(smodel);
            Fnew{sd}.nbddl  = getn(smodel);
        end
        
        F=[F Fnew];
        
    else
        error('la dimension d n''est pas encore prise en compte...')
    end
end

SM.F=F;

if isempty(SM.mapping)
    SM.masse  = cellfun(@(f) f.masse ,F,'UniformOutput',0);
    SM.masse  = SEP(SM.masse);
    SM.metric = cellfun(@(f) f.metric,F,'UniformOutput',0);
    SM.metric = SEP(SM.metric);
else
    SM.masse = calc_trimatrix(SM,TRILINFORM(0,0,0),[],[],[]);
    SM.metric= calc_matrix(SM,BILINFORM(0,0),[],[]);
end





% 
% 
% for d=1:SM.dim
%     if isa(SM.F{d},'double') && length(SM.F{d})==1
%         % Proceder par colocation, point par point
%         % -> DECONSEILLE...
%         
%         nddl = SM.F{d};
%         SM.F{d}={};
%         SM.F{d}.model=nddl;
%         J=[1:nddl]';
%         SM.F{d}.masse     = SPARSETENSOR(repmat(J,1,3),ones(nddl,1));
%         SM.F{d}.metric    = sparse(J,J,ones(nddl,1));
%         SM.F{d}.nbddl     = nddl;
%         SM.F{d}.isrand    = 0; % ATTENTION : ambiguite...
%         SM.F{d}.isspatial = 0; % ATTENTION : ambiguite...
%     else
%         % On donne des vrais modeles :
%         
%         model=SM.F{d};
%         SM.F{d}={};
%         SM.F{d}.model=model;
%         % METRIQUE :
%         try    SM.F{d}.metric= getmetric(model);
%         catch, SM.F{d}.metric= getmetric(calc_metric(model));
%         end
%         % MASSE :
%         try    SM.F{d}.masse = getmasse(model);
%         catch, SM.F{d}.masse = getmasse(calc_masse(model));
%         end
%         
%         % Quelques infos utiles sur les dimensions
%         if isa(SM.F{d}.model,'MODEL')
%             SM.F{d}.nbddl=getnbddl(SM.F{d}.model);
%             SM.F{d}.isrand=0;
%             SM.F{d}.isspatial=1;
%         else
%             SM.F{d}.nbddl=getn(SM.F{d}.model);
%             SM.F{d}.isrand=1;
%             SM.F{d}.isspatial=0;
%         end
%     end
% end
% 
% 
% SM.masse = calc_trimatrix(SM,TRILINFORM(0,0,0),[],[],[]);
% SM.metric= calc_matrix(SM,BILINFORM(0,0),[],[]);
% 
% 
% 
% 
% 


