function ddl = ddlenrich(ls,ddl,nature)

switch nature(1:3)
    case 'tip' 
      k = str2num(nature(4));
      switch getenrichtypetip(ls,k)
          case 1
      ddl = enrich(ddl);
          otherwise
              error('enrichtype non prevu pour l''enrichissement tip')            
      end
        
    case 'sup'
      switch getenrichtypesupport(ls)
          case 1
        ddl = enrich(ddl);  
          otherwise
              error('enrichtype non prevu pour l''enrichissement support')            
      end  
        
end