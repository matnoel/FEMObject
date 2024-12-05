function [H,s]=mat2vect(H)
s=size(H);
for r=1:H.m
    for d=1:H.dim
        H.F{r,d}=mat2vect(H.F{r,d});
    end
end