function [s,I] = isin(a,b)

[s,I] = ismember(a,b);
s=all(s);
