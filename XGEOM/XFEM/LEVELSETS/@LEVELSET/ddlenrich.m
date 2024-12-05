function ddl = ddlenrich(ls,ddl,nature)

switch getnature(ls)

    case 'material'
        ddl=enrich(ddl);
        
end