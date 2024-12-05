function [x,flag] = solve(A,b,varargin)
% function [x,flag] = solve(A,b,varargin)

x = A\b;
flag = 0;
