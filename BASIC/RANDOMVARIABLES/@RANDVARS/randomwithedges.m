function A = randomwithedges(u,varargin)

A = random(u,varargin{:});
B = getedges(u,'cell');

for k=1:u.M
    A{k} = [A{k};B{k}];
end

    