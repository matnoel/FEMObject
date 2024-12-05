modT = QPModelTest() ;
cAT = QPConductivityAssemblerTest() ;
dAT = QPDiffusionAssemblerTest() ;
psT = QPPatchesTest() ;
pbT = QPDiffusionProblemTest() ;

modRes = run(modT) ;
cARes = run(cAT) ;
dARes = run(dAT) ;
psRes = run(psT) ;
pbRes = run(pbT) ;

%% Display

modIncomplete = find([modRes.Incomplete]) ;
modFailed = setdiff(find([modRes.Failed]),modIncomplete) ;
cAIncomplete = find([cARes.Incomplete]) ;
cAFailed = setdiff(find([cARes.Failed]),cAIncomplete) ;
dAIncomplete = find([dARes.Incomplete]) ;
dAFailed = setdiff(find([dARes.Failed]),dAIncomplete) ;
psIncomplete = find([psRes.Incomplete]) ;
psFailed = setdiff(find([psRes.Failed]),psIncomplete) ;
pbIncomplete = find([pbRes.Incomplete]) ;
pbFailed = setdiff(find([pbRes.Failed]),pbIncomplete) ;

if any([modIncomplete cAIncomplete dAIncomplete psIncomplete pbIncomplete])
    disp('_Incomplete tests:')
    for i=modIncomplete; disp(modRes(i).Name) ; end
    for i=cAIncomplete; disp(cARes(i).Name) ; end
    for i=dAIncomplete; disp(dARes(i).Name) ; end
    for i=psIncomplete; disp(psRes(i).Name) ; end
    for i=pbIncomplete; disp(pbRes(i).Name) ; end
else
    disp('All tests complete.')
end

if any([modFailed cAFailed dAFailed psFailed pbFailed])
    disp('_Failed (yet complete) tests:')
    for i=modFailed; disp(modRes(i).Name) ; end
    for i=cAFailed; disp(cARes(i).Name) ; end
    for i=dAFailed; disp(dARes(i).Name) ; end
    for i=psFailed; disp(psRes(i).Name) ; end
    for i=pbFailed; disp(pbRes(i).Name) ; end
else
    disp('Every complete test passed.')
end