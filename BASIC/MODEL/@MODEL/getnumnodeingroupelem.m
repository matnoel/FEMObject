function num=getnumnodeingroupelem(M,i)
% function num=getnumnodeingroupelem(M,i)

 num=getconnec(getgroupelem(M,i));
 num = unique(num);
