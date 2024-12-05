function type = vtkelemtype(elem)

if isa(elem,'ELEMPOINT')
    type = 1;
elseif isa(elem,'SEG2')
    type = 3;
elseif isa(elem,'TRI3')
    type = 5;
elseif isa(elem,'QUA4')
    type = 9;
elseif isa(elem,'TET4')
    type = 10;
elseif isa(elem,'CUB8')
    type = 12;
elseif isa(elem,'PRI6')
    type = 13;
elseif isa(elem,'PYR5')
    type = 14;
elseif isa(elem,'SEG3')
    type = 21;
elseif isa(elem,'TRI6')
    type = 22;
elseif isa(elem,'QUA8') || isa(elem,'QUA9')
    type = 23;
elseif isa(elem,'TET10')
    type = 24;
elseif isa(elem,'CUB20') || isa(elem,'CUB27')
    type = 25;
elseif isa(elem,'PRI15')
    type = 26;
else
    type = 0;
end

end

