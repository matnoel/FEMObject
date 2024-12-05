function h = precalc_intxn(h,p)

for k=1:h.M
h.h{k} = precalc_intxn(h.h{k},p(k));
end
