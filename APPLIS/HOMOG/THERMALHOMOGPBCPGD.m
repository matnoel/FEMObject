classdef THERMALHOMOGPBCPGD
    % THERMALHOMOGPBCPGD class
    %   Homogenization of thermal conductivity using PERIODIC boundary
    %   conditions
    %   A XYZ separation is used here
    properties
        size_domain;
        nb_elem;
        dim
        domain;
        model;
        X;Y;Z;
        ls;
        epsi;
        A_unfree;
        A_unfree_pen;
        BC; %modes macro
        u_unfree_tot; %modes totaux
        u_free_tot;
        u_unfree; %modes micro
        u_free;
        q;
        q_tot;
        result;
        result_svd;
        C;
    end
    methods
        function obj=THERMALHOMOGPBCPGD(size_domain,nb_elem,ls,neg_side,pos_side,delta,epsi,varargin)
            obj.dim=length(size_domain);
            if obj.dim~=2 && obj.dim~=3
                error('Wrong number of dimensions in size_domain')
            end
            obj.size_domain=size_domain;
            if length(nb_elem)~=obj.dim
                error('Wrong number of dimensions in size_domain')
            end
            obj.nb_elem=nb_elem;
            obj.epsi=epsi;
            
            obj.domain=cell(obj.dim,1);
            obj.model=cell(obj.dim,1);
            
            obj.BC=cell(1,obj.dim);
            obj.u_unfree_tot=cell(obj.dim,1);
            obj.u_free_tot=cell(obj.dim,1);
            obj.u_unfree=cell(obj.dim,1);
            obj.u_free=cell(obj.dim,1);
            obj.q=cell(obj.dim,1);
            obj.q_tot=cell(obj.dim,1);
            
            obj.result=cell(obj.dim,1);
            obj.C=zeros(obj.dim);
            
            obj=obj.buildmodels();
            obj=obj.definebcs();
            obj=obj.processls(ls,neg_side,pos_side,delta,varargin{:});
            obj=obj.buildoperators();
        end
        function obj=solve(obj,varargin)
            for i=1:obj.dim
                obj=obj.solve_mode(i,varargin{:});
            end
        end
        
        function obj=solve_mode(obj,i,varargin)
            b=obj.A_unfree*obj.BC{i};
            A=obj.A_unfree_pen;
            
            BC_free=obj.BC{i};
            for j=1:obj.dim
                b=freevector(b,j,obj.model{j});
                A=freematrix(A,j,obj.model{j});
                BC_free=freevector(BC_free,j,obj.model{j});
            end
            solver = SEPSOLVER(getdim(A),varargin{:});
            [obj.u_free{i},obj.result{i}] = solve(A,-b,solver);
            
            obj.u_free_tot{i}=BC_free+obj.u_free{i};
            
            obj.u_unfree{i}=obj.u_free{i};
            for j=1:obj.dim
                obj.u_unfree{i}=unfreevector(obj.u_unfree{i},j,obj.model{j});
            end
            obj.u_unfree_tot{i}=obj.BC{i}+obj.u_unfree{i};
        end
        
        function obj=findeffectiveconductivity(obj,varargin)
            volume=1;
            for i=1:obj.dim
                volume=volume*getvolume(obj.domain{i});
            end
            for i=1:obj.dim
                for j=1:obj.dim
                    Au=multisvd(expand(obj.A_unfree*obj.u_unfree_tot{j}),varargin{:});
                    obj.C(i,j)=prodscal(obj.u_unfree_tot{i},Au);
                end
            end
            obj.C=obj.C/volume;
        end
        
        function obj=computemicroflux(obj,varargin)
            for i=1:obj.dim
                obj=obj.computemicroflux_mode(i,varargin{:});
            end
        end
        
        function obj=computemicroflux_mode(obj,i,varargin)
            if obj.dim==2
                ex=[1;0];
                ey=[0;1];
                hx=obj.size_domain(1)/obj.nb_elem(1);
                hy=obj.size_domain(2)/obj.nb_elem(2);
                %Derivative operator
                Dx=diag(-1/hx*ones(obj.nb_elem(1)+1,1))+diag(1/hx*ones(obj.nb_elem(1),1),1);
                Dx(end,:)=[];
                Dy=diag(-1/hy*ones(obj.nb_elem(2)+1,1))+diag(1/hy*ones(obj.nb_elem(2),1),1);
                Dy(end,:)=[];
                %Mean operator
                Mx=diag(1/2*ones(obj.nb_elem(1)+1,1))+diag(1/2*ones(obj.nb_elem(1),1),1);
                Mx(end,:)=[];
                My=diag(1/2*ones(obj.nb_elem(2)+1,1))+diag(1/2*ones(obj.nb_elem(2),1),1);
                My(end,:)=[];
                
                dx=SEPMATRIX(2);
                dx.alpha=obj.u_unfree{i}.alpha;
                dx.F=obj.u_unfree{i}.F;
                dy=SEPMATRIX(2);
                dy.alpha=obj.u_unfree{i}.alpha;
                dy.F=obj.u_unfree{i}.F;
                
                %Mean LS
                lsc=obj.ls;
                for j=1:lsc.m;
                    lsc.F{j,1}=Mx*obj.ls.F{j,2};
                    lsc.F{j,2}=My*obj.ls.F{j,1};
                end
                for j=1:size(dx.alpha,2)
                    dx.F{j,1}=Dx*dx.F{j,1};
                    dx.F{j,2}=My*dx.F{j,2};
                    dy.F{j,1}=Mx*dy.F{j,1};
                    dy.F{j,2}=Dy*dy.F{j,2};
                end
                qx=SEPMATRIX(3);
                qx.alpha=zeros(1,lsc.m*dx.m);
                l=1;
                for j=1:lsc.m
                    for k=1:dx.m
                        qx.alpha(l)=-lsc.alpha(j)*dx.alpha(k);
                        l=l+1;
                    end
                end
                l=1;
                for j=1:lsc.m
                    for k=1:size(dx.alpha,2)
                        qx.F{l,1}=dx.F{k,1}.*lsc.F{j,1};
                        qx.F{l,2}=dx.F{k,2}.*lsc.F{j,2};
                        qx.F{l,3}=ex;
                        l=l+1;
                    end
                end
                qy=SEPMATRIX(3);
                qy.alpha=zeros(1,lsc.m*dy.m);
                l=1;
                for j=1:lsc.m
                    for k=1:size(dy.alpha,2)
                        qy.alpha(l)=-lsc.alpha(j)*dy.alpha(k);
                        l=l+1;
                    end
                end
                l=1;
                for j=1:lsc.m
                    for k=1:size(dy.alpha,2)
                        qy.F{l,1}=dy.F{k,1}.*lsc.F{j,1};
                        qy.F{l,2}=dy.F{k,2}.*lsc.F{j,2};
                        qy.F{l,3}=ey;
                        l=l+1;
                    end
                end
                obj.q{i}=qx+qy;
                obj.q{i}=multisvd(obj.q{i},'maxorder',obj.q{i}.m,varargin{:});
            elseif obj.dim==3
                ex=[1;0;0];
                ey=[0;1;0];
                ez=[0;0;1];
                
                hx=obj.size_domain(1)/obj.nb_elem(1);
                hy=obj.size_domain(2)/obj.nb_elem(2);
                hz=obj.size_domain(3)/obj.nb_elem(3);
                
                %Derivative operator
                Dx=diag(-1/hx*ones(obj.nb_elem(1)+1,1))+diag(1/hx*ones(obj.nb_elem(1),1),1);
                Dx(end,:)=[];
                Dy=diag(-1/hy*ones(obj.nb_elem(2)+1,1))+diag(1/hy*ones(obj.nb_elem(2),1),1);
                Dy(end,:)=[];
                Dz=diag(-1/hz*ones(obj.nb_elem(3)+1,1))+diag(1/hz*ones(obj.nb_elem(3),1),1);
                Dz(end,:)=[];
                %Mean operator
                Mx=diag(1/2*ones(obj.nb_elem(1)+1,1))+diag(1/2*ones(obj.nb_elem(1),1),1);
                Mx(end,:)=[];
                My=diag(1/2*ones(obj.nb_elem(2)+1,1))+diag(1/2*ones(obj.nb_elem(2),1),1);
                My(end,:)=[];
                Mz=diag(1/2*ones(obj.nb_elem(3)+1,1))+diag(1/2*ones(obj.nb_elem(3),1),1);
                Mz(end,:)=[];
                
                
                dx=SEPMATRIX(3);
                dx.alpha=obj.u_unfree{i}.alpha;
                dx.F=obj.u_unfree{i}.F;
                dy=SEPMATRIX(3);
                dy.alpha=obj.u_unfree{i}.alpha;
                dy.F=obj.u_unfree{i}.F;
                dz=SEPMATRIX(3);
                dz.alpha=obj.u_unfree{i}.alpha;
                dz.F=obj.u_unfree{i}.F;
                
                %Mean LS
                lsc=obj.ls;
                for j=1:lsc.m;
                    lsc.F{j,1}=Mx*obj.ls.F{j,2};
                    lsc.F{j,2}=My*obj.ls.F{j,1};
                    lsc.F{j,3}=Mz*obj.ls.F{j,3};
                end
                for j=1:size(dx.alpha,2)
                    dx.F{j,1}=Dx*dx.F{j,1};
                    dx.F{j,2}=My*dx.F{j,2};
                    dx.F{j,3}=Mz*dx.F{j,3};
                    dy.F{j,1}=Mx*dy.F{j,1};
                    dy.F{j,2}=Dy*dy.F{j,2};
                    dy.F{j,3}=Mz*dy.F{j,3};
                    dz.F{j,1}=Mx*dz.F{j,1};
                    dz.F{j,2}=My*dz.F{j,2};
                    dz.F{j,3}=Dz*dz.F{j,3};
                end
                qx=SEPMATRIX(4);
                qx.alpha=zeros(1,lsc.m*dx.m);
                l=1;
                for j=1:lsc.m
                    for k=1:dx.m
                        qx.alpha(l)=-lsc.alpha(j)*dx.alpha(k);
                        l=l+1;
                    end
                end
                l=1;
                for j=1:lsc.m
                    for k=1:size(dx.alpha,2)
                        qx.F{l,1}=dx.F{k,1}.*lsc.F{j,1};
                        qx.F{l,2}=dx.F{k,2}.*lsc.F{j,2};
                        qx.F{l,3}=dx.F{k,3}.*lsc.F{j,3};
                        qx.F{l,4}=ex;
                        l=l+1;
                    end
                end
                qy=SEPMATRIX(3);
                qy.alpha=zeros(1,lsc.m*dy.m);
                l=1;
                for j=1:lsc.m
                    for k=1:size(dy.alpha,2)
                        qy.alpha(l)=-lsc.alpha(j)*dy.alpha(k);
                        l=l+1;
                    end
                end
                l=1;
                for j=1:lsc.m
                    for k=1:size(dy.alpha,2)
                        qy.F{l,1}=dy.F{k,1}.*lsc.F{j,1};
                        qy.F{l,2}=dy.F{k,2}.*lsc.F{j,2};
                        qy.F{l,3}=dy.F{k,3}.*lsc.F{j,3};
                        qy.F{l,4}=ey;
                        l=l+1;
                    end
                end
                qz=SEPMATRIX(3);
                qz.alpha=zeros(1,lsc.m*dz.m);
                l=1;
                for j=1:lsc.m
                    for k=1:size(dz.alpha,2)
                        qz.alpha(l)=-lsc.alpha(j)*dz.alpha(k);
                        l=l+1;
                    end
                end
                l=1;
                for j=1:lsc.m
                    for k=1:size(dy.alpha,2)
                        qz.F{l,1}=dz.F{k,1}.*lsc.F{j,1};
                        qz.F{l,2}=dz.F{k,2}.*lsc.F{j,2};
                        qz.F{l,3}=dz.F{k,3}.*lsc.F{j,3};
                        qz.F{l,4}=ez;
                        l=l+1;
                    end
                end
                obj.q{i}=qx+qy+qz;
                obj.q{i}=multisvd(obj.q{i},'maxorder',obj.q{i}.m,varargin{:});
            end
        end
        
        function obj=computeflux(obj,varargin)
            for i=1:obj.dim
                obj=obj.computeflux_mode(i,varargin{:});
            end
        end
        
        function obj=computeflux_mode(obj,i,varargin)
            if obj.dim==2
                ex=[1;0];
                ey=[0;1];
                hx=obj.size_domain(1)/obj.nb_elem(1);
                hy=obj.size_domain(2)/obj.nb_elem(2);
                %Derivative operator
                Dx=diag(-1/hx*ones(obj.nb_elem(1)+1,1))+diag(1/hx*ones(obj.nb_elem(1),1),1);
                Dx(end,:)=[];
                Dy=diag(-1/hy*ones(obj.nb_elem(2)+1,1))+diag(1/hy*ones(obj.nb_elem(2),1),1);
                Dy(end,:)=[];
                %Mean operator
                Mx=diag(1/2*ones(obj.nb_elem(1)+1,1))+diag(1/2*ones(obj.nb_elem(1),1),1);
                Mx(end,:)=[];
                My=diag(1/2*ones(obj.nb_elem(2)+1,1))+diag(1/2*ones(obj.nb_elem(2),1),1);
                My(end,:)=[];
                
                dx=SEPMATRIX(2);
                dx.alpha=obj.u_unfree_tot{i}.alpha;
                dx.F=obj.u_unfree_tot{i}.F;
                dy=SEPMATRIX(2);
                dy.alpha=obj.u_unfree_tot{i}.alpha;
                dy.F=obj.u_unfree_tot{i}.F;
                
                %Mean LS
                lsc=obj.ls;
                for j=1:lsc.m;
                    lsc.F{j,1}=Mx*obj.ls.F{j,2};
                    lsc.F{j,2}=My*obj.ls.F{j,1};
                end
                for j=1:size(dx.alpha,2)
                    dx.F{j,1}=Dx*dx.F{j,1};
                    dx.F{j,2}=My*dx.F{j,2};
                    dy.F{j,1}=Mx*dy.F{j,1};
                    dy.F{j,2}=Dy*dy.F{j,2};
                end
                qx=SEPMATRIX(3);
                qx.alpha=zeros(1,lsc.m*dx.m);
                l=1;
                for j=1:lsc.m
                    for k=1:dx.m
                        qx.alpha(l)=-lsc.alpha(j)*dx.alpha(k);
                        l=l+1;
                    end
                end
                l=1;
                for j=1:lsc.m
                    for k=1:size(dx.alpha,2)
                        qx.F{l,1}=dx.F{k,1}.*lsc.F{j,1};
                        qx.F{l,2}=dx.F{k,2}.*lsc.F{j,2};
                        qx.F{l,3}=ex;
                        l=l+1;
                    end
                end
                qy=SEPMATRIX(3);
                qy.alpha=zeros(1,lsc.m*dy.m);
                l=1;
                for j=1:lsc.m
                    for k=1:size(dy.alpha,2)
                        qy.alpha(l)=-lsc.alpha(j)*dy.alpha(k);
                        l=l+1;
                    end
                end
                l=1;
                for j=1:lsc.m
                    for k=1:size(dy.alpha,2)
                        qy.F{l,1}=dy.F{k,1}.*lsc.F{j,1};
                        qy.F{l,2}=dy.F{k,2}.*lsc.F{j,2};
                        qy.F{l,3}=ey;
                        l=l+1;
                    end
                end
                obj.q_tot{i}=qx+qy;
                obj.q_tot{i}=multisvd(obj.q_tot{i},'maxorder',obj.q_tot{i}.m,varargin{:});
            elseif obj.dim==3
                ex=[1;0;0];
                ey=[0;1;0];
                ez=[0;0;1];
                
                hx=obj.size_domain(1)/obj.nb_elem(1);
                hy=obj.size_domain(2)/obj.nb_elem(2);
                hz=obj.size_domain(3)/obj.nb_elem(3);
                
                %Derivative operator
                Dx=diag(-1/hx*ones(obj.nb_elem(1)+1,1))+diag(1/hx*ones(obj.nb_elem(1),1),1);
                Dx(end,:)=[];
                Dy=diag(-1/hy*ones(obj.nb_elem(2)+1,1))+diag(1/hy*ones(obj.nb_elem(2),1),1);
                Dy(end,:)=[];
                Dz=diag(-1/hz*ones(obj.nb_elem(3)+1,1))+diag(1/hz*ones(obj.nb_elem(3),1),1);
                Dz(end,:)=[];
                %Mean operator
                Mx=diag(1/2*ones(obj.nb_elem(1)+1,1))+diag(1/2*ones(obj.nb_elem(1),1),1);
                Mx(end,:)=[];
                My=diag(1/2*ones(obj.nb_elem(2)+1,1))+diag(1/2*ones(obj.nb_elem(2),1),1);
                My(end,:)=[];
                Mz=diag(1/2*ones(obj.nb_elem(3)+1,1))+diag(1/2*ones(obj.nb_elem(3),1),1);
                Mz(end,:)=[];
                
                
                dx=SEPMATRIX(3);
                dx.alpha=obj.u_unfree_tot{i}.alpha;
                dx.F=obj.u_unfree_tot{i}.F;
                dy=SEPMATRIX(3);
                dy.alpha=obj.u_unfree_tot{i}.alpha;
                dy.F=obj.u_unfree_tot{i}.F;
                dz=SEPMATRIX(3);
                dz.alpha=obj.u_unfree_tot{i}.alpha;
                dz.F=obj.u_unfree_tot{i}.F;
                
                %Mean LS
                lsc=obj.ls;
                for j=1:lsc.m;
                    lsc.F{j,1}=Mx*obj.ls.F{j,2};
                    lsc.F{j,2}=My*obj.ls.F{j,1};
                    lsc.F{j,3}=Mz*obj.ls.F{j,3};
                end
                for j=1:size(dx.alpha,2)
                    dx.F{j,1}=Dx*dx.F{j,1};
                    dx.F{j,2}=My*dx.F{j,2};
                    dx.F{j,3}=Mz*dx.F{j,3};
                    dy.F{j,1}=Mx*dy.F{j,1};
                    dy.F{j,2}=Dy*dy.F{j,2};
                    dy.F{j,3}=Mz*dy.F{j,3};
                    dz.F{j,1}=Mx*dz.F{j,1};
                    dz.F{j,2}=My*dz.F{j,2};
                    dz.F{j,3}=Dz*dz.F{j,3};
                end
                qx=SEPMATRIX(4);
                qx.alpha=zeros(1,lsc.m*dx.m);
                l=1;
                for j=1:lsc.m
                    for k=1:dx.m
                        qx.alpha(l)=-lsc.alpha(j)*dx.alpha(k);
                        l=l+1;
                    end
                end
                l=1;
                for j=1:lsc.m
                    for k=1:size(dx.alpha,2)
                        qx.F{l,1}=dx.F{k,1}.*lsc.F{j,1};
                        qx.F{l,2}=dx.F{k,2}.*lsc.F{j,2};
                        qx.F{l,3}=dx.F{k,3}.*lsc.F{j,3};
                        qx.F{l,4}=ex;
                        l=l+1;
                    end
                end
                qy=SEPMATRIX(3);
                qy.alpha=zeros(1,lsc.m*dy.m);
                l=1;
                for j=1:lsc.m
                    for k=1:size(dy.alpha,2)
                        qy.alpha(l)=-lsc.alpha(j)*dy.alpha(k);
                        l=l+1;
                    end
                end
                l=1;
                for j=1:lsc.m
                    for k=1:size(dy.alpha,2)
                        qy.F{l,1}=dy.F{k,1}.*lsc.F{j,1};
                        qy.F{l,2}=dy.F{k,2}.*lsc.F{j,2};
                        qy.F{l,3}=dy.F{k,3}.*lsc.F{j,3};
                        qy.F{l,4}=ey;
                        l=l+1;
                    end
                end
                qz=SEPMATRIX(3);
                qz.alpha=zeros(1,lsc.m*dz.m);
                l=1;
                for j=1:lsc.m
                    for k=1:size(dz.alpha,2)
                        qz.alpha(l)=-lsc.alpha(j)*dz.alpha(k);
                        l=l+1;
                    end
                end
                l=1;
                for j=1:lsc.m
                    for k=1:size(dy.alpha,2)
                        qz.F{l,1}=dz.F{k,1}.*lsc.F{j,1};
                        qz.F{l,2}=dz.F{k,2}.*lsc.F{j,2};
                        qz.F{l,3}=dz.F{k,3}.*lsc.F{j,3};
                        qz.F{l,4}=ez;
                        l=l+1;
                    end
                end
                obj.q_tot{i}=qx+qy+qz;
                obj.q_tot{i}=multisvd(obj.q_tot{i},'maxorder',obj.q_tot{i}.m,varargin{:});
            end
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
                [N M] = size(uu); %The expanded LS is in the wrong order
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
                for d=1:obj.dim
                    uu=expand(obj.u_unfree{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="ModeMicro' num2str(d) '" format="ascii" NumberOfComponents="1">' nl]);
                    for j=1:N
                        for i=1:M
                            fwrite(fid, [num2str(uu(i,j)) ' ']);
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                for d=1:obj.dim
                    uu=expand(obj.u_unfree_tot{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="ModeTot' num2str(d) '" format="ascii" NumberOfComponents="1">' nl]);
                    for j=1:N
                        for i=1:M
                            fwrite(fid, [num2str(uu(i,j)) ' ']);
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                fwrite(fid, ['</PointData>' nl]);
                fwrite(fid, ['<CellData>' nl]);
                for d=1:obj.dim
                    uu=expand(obj.q{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="FluxMicro' num2str(d) '" format="ascii" NumberOfComponents="3">' nl]);
                    for j=1:(N-1)
                        for i=1:(M-1)
                            fwrite(fid, [num2str(uu(i,j,1)) ' ' num2str(uu(i,j,2)) ' ' num2str(0) ' ']);
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                for d=1:obj.dim
                    uu=expand(obj.q_tot{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="FluxTot' num2str(d) '" format="ascii" NumberOfComponents="3">' nl]);
                    for j=1:(N-1)
                        for i=1:(M-1)
                            fwrite(fid, [num2str(uu(i,j,1)) ' ' num2str(uu(i,j,2)) ' ' num2str(0) ' ']);
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                fwrite(fid, ['</CellData>' nl]);
            else
                uu=expand(obj.ls);
                [N M P] = size(uu); %The expanded LS is in the wrong order
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
                for d=1:obj.dim
                    uu=expand(obj.u_unfree{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="ModeMicro' num2str(d) '" format="ascii" NumberOfComponents="1">' nl]);
                    for k=1:P
                        for j=1:N
                            for i=1:M
                                fwrite(fid, [num2str(uu(i,j,k)) ' ']);
                            end
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                for d=1:obj.dim
                    uu=expand(obj.u_unfree_tot{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="ModeTot' num2str(d) '" format="ascii" NumberOfComponents="1">' nl]);
                    for k=1:P
                        for j=1:N
                            for i=1:M
                                fwrite(fid, [num2str(uu(i,j,k)) ' ']);
                            end
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                fwrite(fid, ['</PointData>' nl]);
                fwrite(fid, ['<CellData>' nl]);
                for d=1:obj.dim
                    uu=expand(obj.q{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="FluxMicro' num2str(d) '" format="ascii" NumberOfComponents="3">' nl]);
                    for k=1:(P-1)
                        for j=1:(N-1)
                            for i=1:(M-1)
                                fwrite(fid, [num2str(uu(i,j,k,1)) ' ' num2str(uu(i,j,k,2)) ' ' num2str(uu(i,j,k,3)) ' ']);
                            end
                        end
                    end
                    fwrite(fid,nl);
                    fwrite(fid, ['</DataArray>' nl]);
                end
                for d=1:obj.dim
                    uu=expand(obj.q_tot{d});
                    fwrite(fid, ['<DataArray type="Float32" Name="FluxTot' num2str(d) '" format="ascii" NumberOfComponents="3">' nl]);
                    for k=1:(P-1)
                        for j=1:(N-1)
                            for i=1:(M-1)
                                fwrite(fid, [num2str(uu(i,j,k,1)) ' ' num2str(uu(i,j,k,2)) ' ' num2str(uu(i,j,k,3)) ' ']);
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
            %build 1D models
            for i=1:obj.dim
                obj.domain{i}=DOMAIN(1,0,obj.size_domain(i));
                obj.model{i}=mesh(obj.domain{i},obj.nb_elem(i));
                obj.model{i}=createddlnode(obj.model{i},DDL('T'));
            end
        end
        
        function obj=definebcs(obj)
            %build 1D bcs
            if obj.dim==2
                x = getcoord(getnode(obj.model{1}));
                y = getcoord(getnode(obj.model{2}));
                
                obj.BC{1} = SEPMATRIX({x(:),ones(length(y),1)});
                obj.BC{2} = SEPMATRIX({ones(length(x),1),y(:)});
            elseif obj.dim==3
                x = getcoord(getnode(obj.model{1}));
                y = getcoord(getnode(obj.model{2}));
                z = getcoord(getnode(obj.model{3}));
                
                obj.BC{1} = SEPMATRIX({x(:),ones(length(y),1),ones(length(z),1)});
                obj.BC{2} = SEPMATRIX({ones(length(x),1),y(:),ones(length(z),1)});
                obj.BC{3} = SEPMATRIX({ones(length(x),1),ones(length(y),1),z(:)});
            end
            for i=1:obj.dim
                obj.model{i}=addclperiodic(obj.model{i},POINT(0),POINT(obj.size_domain(i)),'T');
            end
        end
        
        function obj=processls(obj,ls,neg_side,pos_side,delta,varargin)
            ls=neg_side+(pos_side-neg_side)*((1+tanh(ls/delta))/2);
            %SVD
            [obj.ls,obj.result_svd]=multisvd(ls,varargin{:});
        end
        
        function obj=buildoperators(obj)
            if obj.dim==2
                obj.A_unfree=SEPMATRIX(2);
                for i=1:obj.ls.m
                    Ax = BILINFORM(1,1,obj.ls.F{i,2},0);
                    Ax = setfree(Ax,0);
                    Ax = Ax{obj.model{1}}(:,:);
                    Ay = BILINFORM(1,1,obj.ls.F{i,1},0);
                    Ay = setfree(Ay,0);
                    Ay = Ay{obj.model{2}}(:,:);
                    Bx = BILINFORM(0,0,obj.ls.F{i,2},0);
                    Bx = setfree(Bx,0);
                    Bx = Bx{obj.model{1}}(:,:);
                    By = BILINFORM(0,0,obj.ls.F{i,1},0);
                    By = setfree(By,0);
                    By = By{obj.model{2}}(:,:);
                    obj.A_unfree = obj.A_unfree + SEPMATRIX({Ax,By;Bx,Ay},[obj.ls.alpha(i),obj.ls.alpha(i)]);
                end
                % Penalisation
                lx = setfree( LINFORM(0,1),0);
                Lx = calc_vector(lx,obj.model{1});
                Ly = calc_vector(lx,obj.model{2});
                
                obj.A_unfree_pen = obj.A_unfree + SEPMATRIX({Lx*Lx',Ly*Ly'},obj.epsi);
                
            elseif obj.dim==3
                obj.A_unfree=SEPMATRIX(3);
                for i=1:obj.ls.m
                    Ax = BILINFORM(1,1,obj.ls.F{i,2},0);
                    Ax = setfree(Ax,0);
                    Ax = Ax{obj.model{1}}(:,:);
                    Ay = BILINFORM(1,1,obj.ls.F{i,1},0);
                    Ay = setfree(Ay,0);
                    Ay = Ay{obj.model{2}}(:,:);
                    Az = BILINFORM(1,1,obj.ls.F{i,3},0);
                    Az = setfree(Az,0);
                    Az = Az{obj.model{3}}(:,:);
                    Bx = BILINFORM(0,0,obj.ls.F{i,2},0);
                    Bx = setfree(Bx,0);
                    Bx = Bx{obj.model{1}}(:,:);
                    By = BILINFORM(0,0,obj.ls.F{i,1},0);
                    By = setfree(By,0);
                    By = By{obj.model{2}}(:,:);
                    Bz = BILINFORM(0,0,obj.ls.F{i,3},0);
                    Bz = setfree(Bz,0);
                    Bz = Bz{obj.model{3}}(:,:);
                    obj.A_unfree = obj.A_unfree + SEPMATRIX({Ax,By,Bz;Bx,Ay,Bz;Bx,By,Az},[obj.ls.alpha(i),obj.ls.alpha(i),obj.ls.alpha(i)]);
                end
                lx = setfree( LINFORM(0,1),0);
                Lx = calc_vector(lx,obj.model{1});
                Ly = calc_vector(lx,obj.model{2});
                Lz = calc_vector(lx,obj.model{3});
                
                obj.A_unfree_pen = obj.A_unfree + SEPMATRIX({Lx*Lx',Ly*Ly',Lz*Lz'},obj.epsi);
            end
        end
    end
end

