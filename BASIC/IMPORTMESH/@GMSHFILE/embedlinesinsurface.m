function u = embedlinesinsurface(u,numberlines,numbersurface)
% function u = embedlinesinsurface(u,numberlines,numbersurface)

for k=1:length(numberlines)
    u = embedlineinsurface(u,numberlines(k),numbersurface);
end
