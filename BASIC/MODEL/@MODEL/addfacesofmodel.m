function M=addfacesofmodel(M,A)

if isa(M,'MODEL') && isa(A,'MODEL')
    M.facets=[M.facets A.facets];
    M.ridges=[M.ridges A.ridges];
    M.peaks=[M.peaks A.peaks];
end

    