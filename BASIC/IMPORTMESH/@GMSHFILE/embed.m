function u = embed(u,name,number,nameparent,numberparent)
% function u = embed(u,name,number,nameparent,numberparent)

u = stepcounter(u);
s = [name '{' num2str(number) '} In ' nameparent '{' num2str(numberparent) '} ;\n' ];
u = addstring(u,s);
