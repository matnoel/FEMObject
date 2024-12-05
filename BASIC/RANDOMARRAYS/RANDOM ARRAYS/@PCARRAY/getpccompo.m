function u=getpccompo(upc,p)

u=double(upc);
rep=cell(1,ndims(upc));
rep(:)={':'};
rep{upc.pcdim}=p;
u=u(rep{:});
