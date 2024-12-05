function p = calc_gaussorder(elem,choice,varargin)
% function p = calc_gaussorder(elem,choice,varargin)

switch class(choice)
    case 'char'
        switch choice
            case {'mass','load','damp'}
                p = orderN(elem);
                p = 2*p;
            case 'rigi'
                p = orderB(elem);
                p = 2*p;
            case 'adv'
                p = orderN(elem)+orderB(elem);
            case 'rigitang'
                p = 4;
            case 'rigils'
                p = 2 * (orderN(elem)+orderB(elem));
            case {'massls','dampls'}
                p = 2 * 2*orderN(elem);
            case 'lsload'
                p = 2*orderN(elem) ;
            otherwise
                error('Undefined choice')
        end
    case 'double'
        p = choice;
    otherwise
        error('Wrong type of choice. Choice must be a char or a double')
        
end
