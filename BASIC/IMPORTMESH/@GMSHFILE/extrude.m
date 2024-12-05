function u = extrude(u,vect,name,number,varargin)
% function u = extrude(u,vect,name,number)

u = stepcounter(u);
vect=vect(:)';
if length(vect)==2
    vect=[vect,0];
end
if ischarin('recombine',varargin)
    s = ['Extrude' valuesintobraces(vect) ...
        ' { ' name '{' num2str(number) '} ; Recombine; }\n' ];
else
    s = ['Extrude' valuesintobraces(vect) ...
        ' { ' name '{' num2str(number) '} ; }\n' ];
end
u = addstring(u,s);
