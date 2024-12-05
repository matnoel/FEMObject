function [H]=vect2mat(H)
for r=1:H.m
    for d=1:H.dim
        H.F{r,d}=vect2mat(H.F{r,d});
    end
end