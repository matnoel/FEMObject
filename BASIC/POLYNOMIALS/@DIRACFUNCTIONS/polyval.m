function hval = polyval(h,liste,x)

switch min(size(x))
  case 0
    hval = sparse(0,length(liste));
  case 1
    Q = numel(liste);
    x = repmat(x,1,Q);
    hval = repmat(1:Q,size(x,1),1);
    hval = double(sparse(x == hval));
  otherwise
    hval = polyval(h,liste,x(:));
    hval = reshape(full(hval),[size(x,1),size(x,2),length(liste)]);
end

end
