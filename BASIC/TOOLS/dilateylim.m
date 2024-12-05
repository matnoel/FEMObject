function varargout = dilateylim(fact)

a = ylim;
am = (a(2)+a(1))/2;
da = (a(2)-a(1));
a =  [am-da/2*fact,am+da/2*fact];

if nargout==1
    varargout{1}=a;
else
    ylim(a);
end

