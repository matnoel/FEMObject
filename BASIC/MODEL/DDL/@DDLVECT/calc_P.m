function P = calc_P(ddl,sys)


switch ddl.type
    case 'TRANS'
        P = calc_P(sys);
    case 'ROTA'
        if getindim(sys)==2
            P = 1;
        elseif getindim(sys)==3
            P = calc_P(sys);
        else
            error(' ')
        end
    otherwise
        error(' ')
end
