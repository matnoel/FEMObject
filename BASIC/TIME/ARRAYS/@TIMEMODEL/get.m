function t = get(T,param)
% function t = get(T,param)

t = getfield(struct(T),param);