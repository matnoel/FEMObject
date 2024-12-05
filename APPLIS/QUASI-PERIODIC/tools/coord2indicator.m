function carac = coord2indicator(SIZE,coord,tol)
% carac = coord2carac(SIZE,coord,tol)
% Convert Lambda or IxJ coordinates into characteristics operators
% tol is optional and asks for SVD method, which is slower but provide
% lower rank in some cases ; IxJ coordinates only.

% TODO : implement cell input format for coord, to handle multiple phases.

switch numel(SIZE)
    
    case 1
        
        if nargin > 2
            disp('coord2indicator: SVD method will not be applied in this format')
        end
        
        % For safety
        if isempty(coord)
            carac = {zeros(SIZE,1)} ;
            return
        end
        
        % Juste one column with 1 on lines whose number is in coord, 0
        % everywhere else. Sparsity may be applied.
        carac = sparse(coord,1,1,SIZE,1) ;
        if (numel(coord)/SIZE)>=0.5
            carac = full(carac) ;
        end
        carac = {carac} ;

    case 2
        
        % For safety
        if isempty(coord)
            carac = {zeros(SIZE(1),1) zeros(SIZE(2),1)} ;
            return
        end
        
        % Get I and J coordinates from coord 
        switch size(coord,2)
            case 1 % coord is in Lambda format : convert.
                [I,J] = ind2sub(SIZE,coord) ;
            case 2 % coord is in IxJ format
                I = coord(:,1) ;
                J = coord(:,2) ;
            otherwise
                error('Second argument has incorrect number of columns')
        end
        
        % Convert coordinate into boolean columns and store them in one
        % matrix for each dimension
        I_val  = unique(I) ;
        MI = zeros(SIZE(1),length(I_val)) ;
        MJ = zeros(SIZE(2),length(I_val)) ;
        for k = 1:length(I_val)
            MI(:,k) = sparse(I_val(k),1,1,SIZE(1),1) ;
            % For each coordinate value in I, sum boolean columns for every
            % matching value along J into one column :
            MJ(:,k) = sparse(J(I==I_val(k)),1,1,SIZE(2),1) ;
            % Thus there are no duplicate columns in MI
            % There could be in MJ though, hence the following
        end
        
        % Get indexes Original2Unique and vice versa for MJ rows
        [~,o2u,u2o] = unique(MJ','rows','stable') ;
        % UniqueMJ = MJ(o2u,:) and MJ = UniqueMJ(u2o,:) 
        
        %TODO : Is 'stable' necessary ? Doesn't seem like it.
        
        carac = cell(length(o2u),2) ;
        for k = 1:length(o2u)
            % Store MJ's unique columns
            carac{k,2} = MJ(:,o2u(k)) ;
            % Store sum of every corresponding column 
            carac{k,1} = sum( MI(:,u2o==k) ,2) ;
        end
        % Consequently, there are no duplicates in either carac{:,1} nor 
        % carac{:,2}. However, this methods doesn't cope with the case 
        % where a column is linearly tied to the others (i.e. det=0) 
        % In these cases, SVD method provides lower rank.
        
        % SVD method
        exactMethodIsEfficient = size(carac,1)<=prod(SIZE)*1e-3 ; % Factor 1e-3 is arbitrary
        if nargin > 2 && tol>0 && ~exactMethodIsEfficient
            % TODO : obsolete, use TENSALG's Truncator
            carac = SEPMATRIX(2) ;
            for k = 1:length(I)
                MI = sparse(I(k),1,1,SIZE(1),1) ;
                MJ = sparse(J(k),1,1,SIZE(2),1) ;
                carac = carac + SEPMATRIX({MI MJ}) ;
            end
            carac = svd(carac,tol) ;
            carac = getF(carac) ;
        elseif nargin > 2 && tol>0 && exactMethodIsEfficient
            % Ignore if exact method is "efficient enough"
            disp('WARNING coord2indicator: SVD method overriden')
        end
        
    case 3
        error('Not implemented')
        % It would be easy enough to paste the previous code here and add
        % an MK for the third coordinates. Use multisvd instead of svd
        % though.
        % It would be better to write a unique generic code designed to
        % handle any number of dimensions.
        
    otherwise
        error('First argument has incorrect number of elements')

end

end