function u = createboxfield(u,VIn,VOut,XMin,XMax,YMin,YMax,ZMin,ZMax,Thickness,number)
% function u = createboxfield(u,VIn,VOut,XMin,XMax,YMin,YMax,ZMin,ZMax,Thickness,number)

if nargin<2 || isempty(VIn)
    VIn = 1e+22;
end
if nargin<3 || isempty(VOut)
    VOut = 1e+22;
end
if nargin<4 || isempty(XMin)
    XMin = 0;
end
if nargin<5 || isempty(XMax)
    XMax = 0;
end
if nargin<6 || isempty(YMin)
    YMin = 0;
end
if nargin<7 || isempty(YMax)
    YMax = 0;
end
if nargin<8 || isempty(ZMin)
    ZMin = 0;
end
if nargin<9 || isempty(ZMax)
    ZMax = 0;
end
if nargin<10 || isempty(Thickness)
    Thickness = 0;
end
if nargin<11 || isempty(number)
    number = 1;
end

u = createfield(u,'Box',number);
u = setfieldattribute(u,'VIn',VIn,number);
u = setfieldattribute(u,'VOut',VOut,number);
u = setfieldattribute(u,'XMin',XMin,number);
u = setfieldattribute(u,'XMax',XMax,number);
u = setfieldattribute(u,'YMin',YMin,number);
u = setfieldattribute(u,'YMax',YMax,number);
u = setfieldattribute(u,'ZMin',ZMin,number);
u = setfieldattribute(u,'ZMax',ZMax,number);
u = setfieldattribute(u,'Thickness',Thickness,number);
