function conn = mesoConnectivity(model,direction,isPeriodic)
% conn = mesoConnectivity(model,direction,isPeriodic)
% Shortcut to call tool function mesoCOnnectivity with
% correct parameters.

    conn = mesoConnectivity(getOrder(model),getCellNum(model), ...
                            direction,isPeriodic) ;

end