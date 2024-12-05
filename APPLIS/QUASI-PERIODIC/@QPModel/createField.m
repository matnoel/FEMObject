function t = createField(model,generator)
% t = createField(model,generator)

t = TuckerLikeTensor.create(generator,tensorSize(model)) ;
end