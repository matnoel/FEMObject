classdef THERMALHOMOGEBC
    % THERMALHOMOGEBC class
    %   Homogenization of thermal conductivity using ESSENTIAL boundary
    %   conditions
    
    properties
        size_domain;
        nb_elem;
        dim;
        domain;
        model;
        ls;
        A_unfree;
        A_free;
        u_unfree;
        u_free
        C;
    end
    
    methods
        function obj=THERMALHOMOGEBC(size_domain,nb_elem,ls,neg_side,pos_side)
            obj.dim=length(size_domain);
            obj.size_domain=size_domain;
            if (obj.dim ~= 2) && (obj.dim ~=3)
                error('Wrong number of dimensions in size_domain')
            end
            if (length(size_domain)~=length(nb_elem)) && (length(nb_elem)~=1)
                error('Error in the number of dimensions in nb_elem')
            end
            if length(nb_elem)==1
                obj.nb_elem=nb_elem*eyes(obj.dim,1);
            else
                obj.nb_elem=nb_elem;
            end
            
            obj=obj.buildmodels();
            obj=obj.definebcs();
            
            obj=obj.processls(ls,neg_side,pos_side);
            
            obj=obj.buildoperators();
            
            obj.u_unfree=cell([obj.dim 1]);
            obj.u_free=cell([obj.dim 1]);
            obj.C=zeros(obj.dim);
        end
        
        function obj=solve(obj)
            for i=1:obj.dim
                obj=obj.solve_mode(i);
            end
        end
        
        function obj=solve_mode(obj,dim)
            f0=calc_nonhomogeneous_vector(obj.model{dim},obj.A_unfree);
            obj.u_free{dim}=solvesingular(obj.A_free,-f0);
            obj.u_unfree{dim}=unfreevector(obj.model{dim},obj.u_free{dim});
        end
        
        function obj=findeffectiveconductivity(obj)
            for i=1:obj.dim
                for j=1:obj.dim
                    obj.C(i,j)=obj.u_unfree{i}'*obj.A_unfree*obj.u_unfree{j}/getvolume(obj.domain);
                end
            end
        end
        
        function []=plot_mode(obj,dim,varargin)
            dim_str=['x' 'y' 'z'];
            figure('Name',['Mode ' dim_str(dim)]);
            plot(obj.u_unfree{dim},obj.model{dim},varargin{:});
            colorbar;
        end
        
        function []=plot_ls(obj,varargin)
            figure('Name','Level-set');
            plot(obj.ls,obj.model{1},varargin{:});
            colorbar;
        end
    end
    
    methods (Access = private)
        function obj=buildmodels(obj)
            obj.domain=DOMAIN(obj.dim, zeros(1,obj.dim), obj.size_domain);
            obj.model=cell(obj.dim,1);
            if obj.dim == 2
                S=mesh(obj.domain,obj.nb_elem(1),obj.nb_elem(2));
            else
                S=mesh(obj.domain,obj.nb_elem(1),obj.nb_elem(2),obj.nb_elem(3));
            end
            for i=1:obj.dim
                obj.model{i}=S;
                obj.model{i}=createddlnode(obj.model{i},DDL('T'));
            end
        end
        
        function obj=definebcs(obj)
            if obj.dim==2
                for i=1:obj.dim
                    obj.model{i} = addcl(obj.model{i},getedge(obj.domain,1),'T',@(x) x(:,i));
                    obj.model{i} = addcl(obj.model{i},getedge(obj.domain,3),'T',@(x) x(:,i));
                    obj.model{i} = addcl(obj.model{i},getedge(obj.domain,2),'T',@(x) x(:,i));
                    obj.model{i} = addcl(obj.model{i},getedge(obj.domain,4),'T',@(x) x(:,i));
                end
            elseif obj.dim==3
                for i=1:obj.dim
                    obj.model{i} = addcl(obj.model{i},getface(obj.domain,1),'T',@(x) x(:,i));
                    obj.model{i} = addcl(obj.model{i},getface(obj.domain,3),'T',@(x) x(:,i));
                    obj.model{i} = addcl(obj.model{i},getface(obj.domain,2),'T',@(x) x(:,i));
                    obj.model{i} = addcl(obj.model{i},getface(obj.domain,4),'T',@(x) x(:,i));
                    obj.model{i} = addcl(obj.model{i},getface(obj.domain,5),'T',@(x) x(:,i));
                    obj.model{i} = addcl(obj.model{i},getface(obj.domain,6),'T',@(x) x(:,i));
                end
            end
        end
        
        function obj=processls(obj,ls,neg_side,pos_side)
            ls=lseval(ls,obj.model{1});
            ls=getvalue(ls);
            obj.ls=neg_side+(pos_side-neg_side)*double(ls>0);
        end
        
        function obj=buildoperators(obj)
            a=BILINFORM(1,1,obj.ls(:),0);
            a=setfree(a,0);
            obj.A_unfree=a{obj.model{1}}(:,:);
            obj.A_free=freematrix(obj.model{1},obj.A_unfree);
        end
    end
end

