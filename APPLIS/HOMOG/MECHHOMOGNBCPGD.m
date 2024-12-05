classdef MECHHOMOGNBCPGD
    % MECHHOMOGNBCPGD class
    %   Homogenization of elasticity tensor using NATURAL boundary
    %   conditions
    %   A XYZ separation is used here
    
    properties
        size_domain;
        nb_elem;
        dim;
        nb_homog_mode;
        domain;
        model;
        X;Y;Z;
        ls;
        A_uf;
        A_uf_pen;
        epsilon1;
        epsilon2;
        L;
        mean_op;
        deriv_op;
        lin
        BC;
        u_uf;
        strain;
        stress;
        vonmises;
        strainenergydensity;
        result;
        result_svd;
        C_ps; % Elasticity tensor for the outer (positive) side of the LS
        C_ns; % Elasticity tensor for the inner (negative) side of the LS
        C_s;  % Separate elasticity tensor
        C;    % Final elasticity tensor
    end
    
    methods
        function obj=MECHHOMOGNBCPGD(size_domain,nb_elem,ls,C_ns,C_ps,delta,epsilon1,epsilon2,varargin)
            obj.epsilon1=epsilon1;
            obj.epsilon2=epsilon2;
            obj.dim=length(size_domain);
            if obj.dim~=2 && obj.dim~=3
                error('Wrong number of dimensions in size_domain')
            end
            obj.size_domain=size_domain;
            if length(nb_elem)~=obj.dim
                error('Wrong number of dimensions in size_domain')
            end
            if obj.dim==2
                obj.nb_homog_mode=3;
            elseif obj.dim==3
                obj.nb_homog_mode=6;
            end
            obj.nb_elem=nb_elem;
            
            obj.domain=cell(obj.dim,1);
            obj.model=cell(obj.dim,1);
            
            obj.C_ns=C_ns;
            obj.C_ps=C_ps;
            obj.lin=cell(1,obj.nb_homog_mode);
            obj.BC=cell(1,obj.nb_homog_mode);
            obj.u_uf=cell(obj.nb_homog_mode,1);
            obj.strain=cell(obj.nb_homog_mode,1);
            obj.stress=cell(obj.nb_homog_mode,1);
            obj.vonmises=cell(obj.nb_homog_mode,1);
            obj.strainenergydensity=cell(obj.nb_homog_mode,1);
            obj.result=cell(obj.nb_homog_mode,1);
            obj.C=zeros(obj.nb_homog_mode);
            
            obj=obj.buildmodels();
            obj=obj.computelin();
            obj=obj.definebcs();
            obj=obj.processls(ls,delta,varargin{:});
            obj=obj.buildoperators();
        end
        function obj=solve(obj,varargin)
            for i=1:obj.nb_homog_mode
                obj=obj.solve_mode(i,varargin{:});
            end
        end
        
        function obj=solve_mode(obj,i,varargin)
            A=obj.A_uf_pen;
            for j=1:obj.dim
                A=freematrix(A,j,obj.model{j});
            end
            solver = SEPSOLVER(getdim(A),varargin{:});
            [obj.u_uf{i},obj.result{i}] = solve(A,obj.BC{i},solver);
            for j=1:obj.dim
                obj.u_uf{i}=freevector(obj.u_uf{i},j,obj.model{j});
            end
        end
        
        function obj=findeffectiveproperties(obj,varargin)
            volume=1;
            for i=1:obj.dim
                volume=volume*getvolume(obj.domain{i});
            end
            h=obj.size_domain./obj.nb_elem;
            integ_cell_field_op=SEPMATRIX(obj.dim+1,1);
            for i=1:obj.dim
                integ_cell_field_op.F{i}=h(i)*ones(1,obj.nb_elem(i));
            end
            integ_cell_field_op.F{obj.dim+1}=eye(obj.nb_homog_mode);
            for i=1:obj.nb_homog_mode
                obj.C(:,i)=expand(integ_cell_field_op*obj.strain{i});
            end
            obj.C=obj.C/volume;
            obj.C=inv(obj.C);
        end
        
        function obj=computefields(obj,varargin)
            for i=1:obj.nb_homog_mode
                obj=obj.computefields_mode(i,varargin{:});
            end
        end
        
        function obj=computefields_mode(obj,i,varargin)
            % Compute strain, stress and vonmises
            % Compute gradx u, grady u and gradz u
            % Derivative operator
            obj.mean_op=cell(obj.dim,1);
            obj.deriv_op=cell(obj.dim,1);
            for j=1:obj.dim
                r=obj.nb_elem(j);
                obj.mean_op{j}=diag(1/2*ones(r+1,1))+diag(1/2*ones(r,1),1);
                obj.mean_op{j}(end,:)=[];
                h=obj.size_domain(j)/r;
                obj.deriv_op{j}=diag(-1/h*ones(r+1,1))+diag(1/h*ones(r,1),1);
                obj.deriv_op{j}(end,:)=[];
            end
            
            % STRAIN
            if obj.dim==2
                strainop=SEPMATRIX({obj.deriv_op{1},obj.mean_op{2},obj.L{1};...
                    obj.mean_op{1},obj.deriv_op{2},obj.L{2}});
            elseif obj.dim==3
                strainop=SEPMATRIX({obj.deriv_op{1},obj.mean_op{2},obj.mean_op{3},obj.L{1};...
                    obj.mean_op{1},obj.deriv_op{2},obj.mean_op{3},obj.L{2};...
                    obj.mean_op{1},obj.mean_op{2},obj.deriv_op{3},obj.L{3}...
                    });
            end
            obj.strain{i}=strainop*obj.u_uf{i};
            
            strains=multisvd(obj.strain{i},'maxorder',obj.strain{i}.m,varargin{:});
            if strains.m<obj.strain{i}.m
                obj.strain{i}=strains;
            end
            clear strains;
            
            % STRESS
            centerop=SEPMATRIX(obj.dim+1,1);
            for j=1:obj.dim
                centerop.F{j}=obj.mean_op{j};
            end
            centerop.F{obj.dim+1}=eye(obj.nb_homog_mode);
            
            Cg=diag(centerop*obj.C_s,1:obj.dim);
            obj.stress{i}=Cg*obj.strain{i};
            
            sigx=expand(obj.stress{i});
            stresss=multisvd(sigx,'reference',sigx,'maxorder',obj.stress{i}.m,varargin{:});
            
            if stresss.m<obj.stress{i}.m
                obj.stress{i}=stresss;
            end
            clear stresss;
            
            % VONMISES
            
            vm=zeros(obj.nb_elem);
            if obj.dim==2
                for k=1:obj.nb_elem(1)
                    for l=1:obj.nb_elem(2)
                        vm(k,l)=abs(sqrt(3/2*(sigx(k,l,1)^2+sigx(k,l,2)^2+2*sigx(k,l,3)^2)...
                            -1/2*(sigx(k,l,1)+sigx(k,l,2))^2));
                    end
                end
            else
                for k=1:obj.nb_elem(1)
                    for l=1:obj.nb_elem(2)
                        for m=1:obj.nb_elem(3)
                            % use of 'abs' to remove complex values
                            vm(k,l,m)=abs(sqrt(3/2*(sigx(k,l,m,1)^2+sigx(k,l,m,2)^2+sigx(k,l,m,3)^2+...
                                2*(sigx(k,l,m,4)^2+sigx(k,l,m,5)^2+sigx(k,l,m,6)^2))...
                                -1/2*(sigx(k,l,m,1)+sigx(k,l,m,2)+sigx(k,l,m,2))^2));
                        end
                    end
                end
            end
            obj.vonmises{i}=multisvd(vm,'reference',vm,...
                'maxorder',obj.stress{i}.m,varargin{:});
            clear sigx;
            clear vm;
            % STRAIN ENERGY DENSITY
            sig=transpose(diag(obj.stress{i},1:obj.dim),obj.stress{i}.dim);
            obj.strainenergydensity{i}=multisvd(expand(sig*obj.strain{i}),...
                'maxorder',obj.stress{i}.m*obj.strain{i}.m,varargin{:});
            clear sig;
        end
        
        function []=export(obj,filename)
            fid = fopen(filename, 'w');
            if fid == -1
                error('Cannot open file for writing.');
            end
            nl = sprintf('\n');
            fwrite(fid, ['<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">' nl]);
            if obj.dim==2
                uu=expand(obj.ls);
                [N M] = size(uu); % The expanded LS is in the wrong order
                spacing=(obj.size_domain)./[M N];
                fwrite(fid, ['<ImageData WholeExtent="1 ' num2str(M) ' 1 ' num2str(N) ' 1 ' num2str(1) '"' nl ...
                    'Origin="0 0 0" Spacing="'  num2str(spacing(1)) ' ' num2str(spacing(2)) ' ' num2str(1) '">' nl]);
                fwrite(fid, ['<Piece Extent="1 ' num2str(M) ' 1 ' num2str(N) ' 1 1">' nl]);
                fwrite(fid, ['<PointData>' nl]);
                fwrite(fid, ['<DataArray type="Float32" Name="LevelSet" format="ascii">' nl]);
                for j=1:N
                    for i=1:M
                        fwrite(fid, [num2str(uu(j,i)) ' ']);
                    end
                end
                fwrite(fid,nl);
                fwrite(fid, ['</DataArray>' nl]);
                for d=1:obj.nb_homog_mode
                    uu=expand(obj.u_uf{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="DisplacementMode' num2str(d) '" format="ascii" NumberOfComponents="3">' nl]);
                    for j=1:N
                        for i=1:M
                            fwrite(fid, [num2str(uu(i,j,1)) ' ' num2str(uu(i,j,2)) ' ' num2str(0) ' ']);
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                fwrite(fid, ['</PointData>' nl]);
                fwrite(fid, ['<CellData>' nl]);
                for d=1:obj.nb_homog_mode
                    uu=expand(obj.strain{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="Strain' num2str(d) '" format="ascii" NumberOfComponents="9">' nl]);
                    for j=1:(N-1)
                        for i=1:(M-1)
                            fwrite(fid, [num2str(uu(i,j,1)) ' ' num2str(uu(i,j,3)/2) ' ' num2str(0) ' '...
                                num2str(uu(i,j,3)/2) ' ' num2str(uu(i,j,2)) ' ' num2str(0) ' '...
                                num2str(0) ' ' num2str(0) ' ' num2str(0) ' ']);
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                for d=1:obj.nb_homog_mode
                    uu=expand(obj.stress{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="Stress' num2str(d) '" format="ascii" NumberOfComponents="9">' nl]);
                    for j=1:(N-1)
                        for i=1:(M-1)
                            fwrite(fid, [num2str(uu(i,j,1)) ' ' num2str(uu(i,j,3)) ' ' num2str(0) ' '...
                                num2str(uu(i,j,3)) ' ' num2str(uu(i,j,2)) ' ' num2str(0) ' '...
                                num2str(0) ' ' num2str(0) ' ' num2str(0) ' ']);
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                for d=1:obj.nb_homog_mode
                    uu=expand(obj.vonmises{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="VonMises' num2str(d) '" format="ascii" NumberOfComponents="1">' nl]);
                    for j=1:(N-1)
                        for i=1:(M-1)
                            fwrite(fid, [num2str(uu(i,j)) ' ']);
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                for d=1:obj.nb_homog_mode
                    uu=expand(obj.strainenergydensity{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="StrainEnergyDensity' num2str(d) '" format="ascii" NumberOfComponents="1">' nl]);
                    for j=1:(N-1)
                        for i=1:(M-1)
                            fwrite(fid, [num2str(uu(i,j)) ' ']);
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                fwrite(fid, ['</CellData>' nl]);
            else
                uu=expand(obj.ls);
                [N M P] = size(uu); % The expanded LS is in the wrong order
                spacing=(obj.size_domain)./[M N P];
                fwrite(fid, ['<ImageData WholeExtent="1 ' num2str(M) ' 1 ' num2str(N) ' 1 ' num2str(P) '"' nl ...
                    'Origin="0 0 0" Spacing="'  num2str(spacing(1)) ' ' num2str(spacing(2)) ' ' num2str(spacing(3)) '">' nl]);
                fwrite(fid, ['<Piece Extent="1 ' num2str(M) ' 1 ' num2str(N) ' 1 ' num2str(P) '">' nl]);
                fwrite(fid, ['<PointData>' nl]);
                fwrite(fid, ['<DataArray type="Float32" Name="LevelSet" format="ascii">' nl]);
                for k=1:P
                    for j=1:N
                        for i=1:M
                            fwrite(fid, [num2str(uu(j,i,k)) ' ']);
                        end
                    end
                end
                fwrite(fid,nl);
                fwrite(fid, ['</DataArray>' nl]);
                for d=1:obj.nb_homog_mode
                    uu=expand(obj.u_uf{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="DisplacementMode' num2str(d) '" format="ascii" NumberOfComponents="3">' nl]);
                    for k=1:P
                        for j=1:N
                            for i=1:M
                                fwrite(fid, [num2str(uu(i,j,k,1)) ' ' num2str(uu(i,j,k,2)) ' ' num2str(uu(i,j,k,3)) ' ']);
                            end
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                fwrite(fid, ['</PointData>' nl]);
                fwrite(fid, ['<CellData>' nl]);
                for d=1:obj.nb_homog_mode
                    uu=expand(obj.strain{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="Strain' num2str(d) '" format="ascii" NumberOfComponents="9">' nl]);
                    for k=1:(P-1)
                        for j=1:(N-1)
                            for i=1:(M-1)
                                fwrite(fid, [num2str(uu(i,j,k,1)) ' ' num2str(uu(i,j,k,6)/2) ' ' num2str(uu(i,j,k,5)/2) ' '...
                                    num2str(uu(i,j,k,6)/2) ' ' num2str(uu(i,j,k,2)) ' ' num2str(uu(i,j,k,4)/2) ' '...
                                    num2str(uu(i,j,k,5)/2) ' ' num2str(uu(i,j,k,4)/2) ' ' num2str(uu(i,j,k,3)) ' ']);
                            end
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                for d=1:obj.nb_homog_mode
                    uu=expand(obj.stress{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="Stress' num2str(d) '" format="ascii" NumberOfComponents="9">' nl]);
                    for k=1:(P-1)
                        for j=1:(N-1)
                            for i=1:(M-1)
                                fwrite(fid, [num2str(uu(i,j,k,1)) ' ' num2str(uu(i,j,k,6)) ' ' num2str(uu(i,j,k,5)) ' '...
                                    num2str(uu(i,j,k,6)) ' ' num2str(uu(i,j,k,2)) ' ' num2str(uu(i,j,k,4)) ' '...
                                    num2str(uu(i,j,k,5)) ' ' num2str(uu(i,j,k,4)) ' ' num2str(uu(i,j,k,3)) ' ']);
                            end
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                for d=1:obj.nb_homog_mode
                    uu=expand(obj.vonmises{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="VonMises' num2str(d) '" format="ascii" NumberOfComponents="1">' nl]);
                    for k=1:(P-1)
                        for j=1:(N-1)
                            for i=1:(M-1)
                                fwrite(fid, [num2str(uu(i,j,k)) ' ']);
                            end
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                for d=1:obj.nb_homog_mode
                    uu=expand(obj.strainenergydensity{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="StrainEnergyDensity' num2str(d) '" format="ascii" NumberOfComponents="1">' nl]);
                    for k=1:(P-1)
                        for j=1:(N-1)
                            for i=1:(M-1)
                                fwrite(fid, [num2str(uu(i,j,k)) ' ']);
                            end
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                fwrite(fid, ['</CellData>' nl]);
            end
            fwrite(fid, ['</Piece>' nl]);
            fwrite(fid, ['</ImageData>' nl]);
            fwrite(fid, ['</VTKFile>' nl]);
            fclose(fid);
        end
    end
    
    
    methods (Access = private)
        function obj=buildmodels(obj)
            % build 1D models
            for i=1:obj.dim
                obj.domain{i}=DOMAIN(1,0,obj.size_domain(i));
                obj.model{i}=mesh(obj.domain{i},obj.nb_elem(i));
                obj.model{i}=createddlnode(obj.model{i},DDL('u'));
            end
        end
        
        function obj=computelin(obj)
            l = setfree(LINFORM(0,1),0);
            for i=1:obj.dim
                obj.lin{i}=calc_vector(l,obj.model{i});
            end
        end
        
        function obj=definebcs(obj)
            if obj.dim==2
                xm=zeros(obj.nb_elem(1)+1,1);
                xm(1)=-1;
                xp=zeros(obj.nb_elem(1)+1,1);
                xp(obj.nb_elem(1)+1)=1;
                
                ym=zeros(obj.nb_elem(2)+1,1);
                ym(1)=-1;
                yp=zeros(obj.nb_elem(2)+1,1);
                yp(obj.nb_elem(2)+1)=1;
                
                obj.BC{1} = SEPMATRIX({xm,obj.lin{2},[1;0];xp,obj.lin{2},[1;0]});
                obj.BC{2} = SEPMATRIX({obj.lin{1},ym,[0;1];obj.lin{1},yp,[0;1]});
                obj.BC{3} = SEPMATRIX({xm,obj.lin{2},[0;1];xp,obj.lin{2},[0;1];...
                    obj.lin{1},ym,[1;0];obj.lin{1},yp,[1;0]});
            elseif obj.dim==3
                xm=zeros(obj.nb_elem(1)+1,1);
                xm(1)=-1;
                xp=zeros(obj.nb_elem(1)+1,1);
                xp(obj.nb_elem(1)+1)=1;
                
                ym=zeros(obj.nb_elem(2)+1,1);
                ym(1)=-1;
                yp=zeros(obj.nb_elem(2)+1,1);
                yp(obj.nb_elem(2)+1)=1;
                
                zm=zeros(obj.nb_elem(3)+1,1);
                zm(1)=-1;
                zp=zeros(obj.nb_elem(3)+1,1);
                zp(obj.nb_elem(3)+1)=1;
                
                obj.BC{1} = SEPMATRIX({xm,obj.lin{2},obj.lin{3},[1;0;0];...
                    xp,obj.lin{2},obj.lin{3},[1;0;0]});
                obj.BC{2} = SEPMATRIX({obj.lin{1},ym,obj.lin{3},[0;1;0];...
                    obj.lin{1},yp,obj.lin{3},[0;1;0]});
                obj.BC{3} = SEPMATRIX({obj.lin{1},obj.lin{2},zm,[0;0;1];...
                    obj.lin{1},obj.lin{2},zp,[0;0;1]});
                obj.BC{4} = SEPMATRIX({obj.lin{1},ym,obj.lin{3},[0;0;1];...
                    obj.lin{1},yp,obj.lin{3},[0;0;1];...
                    obj.lin{1},obj.lin{2},zm,[0;1;0];...
                    obj.lin{1},obj.lin{2},zp,[0;1;0]});
                obj.BC{5} = SEPMATRIX({xm,obj.lin{2},obj.lin{3},[0;0;1];...
                    xp,obj.lin{2},obj.lin{3},[0;0;1];...
                    obj.lin{1},obj.lin{2},zm,[1;0;0];...
                    obj.lin{1},obj.lin{2},zp,[1;0;0]});
                obj.BC{6} = SEPMATRIX({xm,obj.lin{2},obj.lin{3},[0;1;0];...
                    xp,obj.lin{2},obj.lin{3},[0;1;0];...
                    obj.lin{1},ym,obj.lin{3},[1;0;0];...
                    obj.lin{1},yp,obj.lin{3},[1;0;0]});
            end
        end
        
        function obj=processls(obj,ls,delta,varargin)
            % ls=neg_side+(pos_side-neg_side)*((1+tanh(ls/delta))/2);
            ls=1-(1+tanh(ls/delta))/2;
            % SVD
            [obj.ls,obj.result_svd]=multisvd(ls,varargin{:});
        end
        
        function obj=buildoperators(obj)
            % Build C_s
            Cdiff=obj.C_ns-obj.C_ps;
            Calpha=ones(1,1+obj.ls.m);
            Calpha(2:end)=obj.ls.alpha;
            
            obj.C_s=SEPMATRIX(obj.dim+1,Calpha);
            for j=1:obj.dim
                obj.C_s.F{1,j}=ones(obj.nb_elem(j)+1,1);
            end
            obj.C_s.F{1,obj.dim+1}=obj.C_ps;
            for i=1:obj.ls.m
                obj.C_s.F{1+i,1}=obj.ls.F{i,2}; % SVD of the LS in the wrong order
                obj.C_s.F{1+i,2}=obj.ls.F{i,1};
                obj.C_s.F{1+i,obj.dim+1}=Cdiff;
            end
            if obj.dim==3
                for i=1:obj.ls.m
                    obj.C_s.F{1+i,3}=obj.ls.F{i,3};
                end
            end
            % Build A_uf
            if obj.dim==2
                obj.L{1}=zeros(3,2);
                obj.L{1}(1,1)=1;
                obj.L{1}(3,2)=1;
                
                obj.L{2}=zeros(3,2);
                obj.L{2}(2,2)=1;
                obj.L{2}(3,1)=1;
                
                obj.A_uf=SEPMATRIX(3);
                
                for i=1:obj.C_s.m
                    Kx=BILINFORM(1,1,obj.C_s.F{i,1},0);
                    Kx=setfree(Kx,0);
                    Kx=Kx{obj.model{1}}(:,:);
                    Mx=BILINFORM(0,0,obj.C_s.F{i,1},0);
                    Mx=setfree(Mx,0);
                    Mx=Mx{obj.model{1}}(:,:);
                    Dlx=BILINFORM(1,0,obj.C_s.F{i,1},0);
                    Dlx=setfree(Dlx,0);
                    Dlx=Dlx{obj.model{1}}(:,:);
                    Drx=Dlx';
                    
                    Ky=BILINFORM(1,1,obj.C_s.F{i,2},0);
                    Ky=setfree(Ky,0);
                    Ky=Ky{obj.model{2}}(:,:);
                    My=BILINFORM(0,0,obj.C_s.F{i,2},0);
                    My=setfree(My,0);
                    My=My{obj.model{2}}(:,:);
                    Dly=BILINFORM(1,0,obj.C_s.F{i,2},0);
                    Dly=setfree(Dly,0);
                    Dly=Dly{obj.model{2}}(:,:);
                    Dry=Dly';
                    
                    d11=obj.L{1}'*obj.C_s.F{i,3}*obj.L{1};
                    d12=obj.L{1}'*obj.C_s.F{i,3}*obj.L{2};
                    d21=obj.L{2}'*obj.C_s.F{i,3}*obj.L{1};
                    d22=obj.L{2}'*obj.C_s.F{i,3}*obj.L{2};
                    
                    obj.A_uf=obj.A_uf+SEPMATRIX({Kx,My,d11;Dlx,Dry,d12;Drx,Dly,d21;Mx,Ky,d22},obj.C_s.alpha(i)*ones(1,4));
                end
                % Penalization
                lx = setfree( LINFORM(0,1),0);
                Lx = calc_vector(lx,obj.model{1});
                Ly = calc_vector(lx,obj.model{2});
                obj.A_uf_pen = obj.A_uf + SEPMATRIX({Lx*Lx',Ly*Ly',eye(2)},obj.epsilon1);
                
                x=getcoord(getnode(obj.model{1}));
                y=getcoord(getnode(obj.model{2}));
                
                mx = setfree( BILINFORM(0,0,1),0);
                Mx = calc_matrix(mx,obj.model{1});
                My = calc_matrix(mx,obj.model{2});
                
                d11 = [1 0;0 0];
                d12 = [0 1;0 0];
                d21 = [0 0;1 0];
                d22 = [0 0;0 1];
                
                obj.A_uf_pen = obj.A_uf +...
                    SEPMATRIX({Lx*Lx', My*(y*y')*My',d11;...
                    Lx*x'*Mx',My*y*Ly',d12;...
                    Mx*x*Lx',Ly*y'*My',d21;...
                    Mx*(x*x')*Mx',Ly*Ly',d22},obj.epsilon2*[1 -1 -1 1]);
                
            elseif obj.dim==3
                obj.L{1}=zeros(6,3);
                obj.L{1}(1,1)=1;
                obj.L{1}(5,3)=1;
                obj.L{1}(6,2)=1;
                
                obj.L{2}=zeros(6,3);
                obj.L{2}(2,2)=1;
                obj.L{2}(4,3)=1;
                obj.L{2}(6,1)=1;
                
                obj.L{3}=zeros(6,3);
                obj.L{3}(3,3)=1;
                obj.L{3}(4,2)=1;
                obj.L{3}(5,1)=1;
                
                obj.A_uf=SEPMATRIX(4);
                
                for i=1:obj.C_s.m
                    Kx=BILINFORM(1,1,obj.C_s.F{i,1},0);
                    Kx=setfree(Kx,0);
                    Kx=Kx{obj.model{1}}(:,:);
                    Mx=BILINFORM(0,0,obj.C_s.F{i,1},0);
                    Mx=setfree(Mx,0);
                    Mx=Mx{obj.model{1}}(:,:);
                    Dlx=BILINFORM(1,0,obj.C_s.F{i,1},0);
                    Dlx=setfree(Dlx,0);
                    Dlx=Dlx{obj.model{1}}(:,:);
                    Drx=Dlx';
                    
                    Ky=BILINFORM(1,1,obj.C_s.F{i,2},0);
                    Ky=setfree(Ky,0);
                    Ky=Ky{obj.model{2}}(:,:);
                    My=BILINFORM(0,0,obj.C_s.F{i,2},0);
                    My=setfree(My,0);
                    My=My{obj.model{2}}(:,:);
                    Dly=BILINFORM(1,0,obj.C_s.F{i,2},0);
                    Dly=setfree(Dly,0);
                    Dly=Dly{obj.model{2}}(:,:);
                    Dry=Dly';
                    
                    Kz=BILINFORM(1,1,obj.C_s.F{i,3},0);
                    Kz=setfree(Kz,0);
                    Kz=Kz{obj.model{3}}(:,:);
                    Mz=BILINFORM(0,0,obj.C_s.F{i,3},0);
                    Mz=setfree(Mz,0);
                    Mz=Mz{obj.model{3}}(:,:);
                    Dlz=BILINFORM(1,0,obj.C_s.F{i,3},0);
                    Dlz=setfree(Dlz,0);
                    Dlz=Dlz{obj.model{3}}(:,:);
                    Drz=Dlz';
                    
                    d11=obj.L{1}'*obj.C_s.F{i,4}*obj.L{1};
                    d12=obj.L{1}'*obj.C_s.F{i,4}*obj.L{2};
                    d13=obj.L{1}'*obj.C_s.F{i,4}*obj.L{3};
                    d21=obj.L{2}'*obj.C_s.F{i,4}*obj.L{1};
                    d22=obj.L{2}'*obj.C_s.F{i,4}*obj.L{2};
                    d23=obj.L{2}'*obj.C_s.F{i,4}*obj.L{3};
                    d31=obj.L{3}'*obj.C_s.F{i,4}*obj.L{1};
                    d32=obj.L{3}'*obj.C_s.F{i,4}*obj.L{2};
                    d33=obj.L{3}'*obj.C_s.F{i,4}*obj.L{3};
                    
                    obj.A_uf=obj.A_uf+...
                        SEPMATRIX({...
                        Kx,My,Mz,d11;Dlx,Dry,Mz,d12;Dlx,My,Drz,d13; ...
                        Drx,Dly,Mz,d21;Mx,Ky,Mz,d22;Mx,Dly,Drz,d23; ...
                        Drx,My,Dlz,d31;Mx,Dry,Dlz,d32;Mx,My,Kz,d33},...
                        obj.C_s.alpha(i)*ones(1,9));
                end
                % Penalization
                lx = setfree( LINFORM(0,1),0);
                Lx = calc_vector(lx,obj.model{1});
                Ly = calc_vector(lx,obj.model{2});
                Lz = calc_vector(lx,obj.model{3});
                obj.A_uf_pen = obj.A_uf + SEPMATRIX({Lx*Lx',Ly*Ly',Lz*Lz',eye(3)},obj.epsilon1);
                
                x=getcoord(getnode(obj.model{1}));
                y=getcoord(getnode(obj.model{2}));
                z=getcoord(getnode(obj.model{3}));
                
                mx = setfree( BILINFORM(0,0,1),0);
                Mx = calc_matrix(mx,obj.model{1});
                My = calc_matrix(mx,obj.model{2});
                Mz = calc_matrix(mx,obj.model{3});
                
                d11 = [1 0 0;0 0 0;0 0 0];
                d12 = [0 1 0;0 0 0;0 0 0];
                d13 = [0 0 1;0 0 0;0 0 0];
                d21 = d12';
                d22 = [0 0 0;0 1 0;0 0 0];
                d23 = [0 0 0;0 0 1;0 0 0];
                d31 = d13';
                d32 = d32';
                d33 = [0 0 0;0 0 0;0 0 1];
                
                obj.A_uf_pen = obj.A_uf +...
                    SEPMATRIX({Lx*Lx',Ly*Ly',Mz*(z*z')*Mz',d22;...
                    Lx*Lx',Ly*y'*My,Mz*z*Lz',d23;...
                    Lx*Lx',My*y*Ly',Lz*z'*Mz',d32;...
                    Lx*Lx',My*(y*y')*My',Lz*Lz',d33;...
                    Lx*Lx',Ly*Ly',Mz*(z*z')*Mz',d11;...
                    Lx*x'*Mx,Ly*Ly',Mz*z*Lz',d13;...
                    Mx*x*Lx',Ly*Ly',Lz*z'*Mz',d31;...
                    Mx*(x*x')*Mx',Ly*Ly',Lz*Lz',d33;...
                    Lx*Lx',My*(y*y')*My',Lz*Lz',d11;...
                    Lx*x'*Mx',My*y*Ly',Lz*Lz',d12;...
                    Mx*x*Lx',Ly*y'*My',Lz*Lz',d21;...
                    Mx*(x*x')*Mx',Ly*Ly',Lz*Lz',d22;...
                    },obj.epsilon2*[1 -1 -1 1 1 -1 -1 1 1 -1 -1 1]);
                
            end
        end
    end
end
