function [fx1,gx1,fx2,gx2] = getEMT(data_cell,data_gene,elist,mlist)

for k=1:numel(data_gene)
    commongene(k)= sum(strcmp(data_gene(k),elist)); % Find overlap with E genes
end
for k=1:numel(data_gene)
    commongene2(k)= sum(strcmp(data_gene(k),mlist)); % Find overlap with M genes
end
clear EMTscoresum
for k=1:numel(data_cell(1,:))
 [fx1,gx1]=ecdf(data_cell(commongene==1,k)); %calculate ECDF of all E genes
[fx2,gx2]=ecdf(data_cell(commongene2==1,k)); %calculate ECDF of all M genes
if max(gx1)>max(gx2)
    fx2(end+1)=1;gx2(end+1)=max(gx1);
elseif max(gx1)<max(gx2)
    fx1(end+1)=1;gx1(end+1)=max(gx2);
end

end

