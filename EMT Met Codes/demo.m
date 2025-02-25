%% Read in your dataset

meth = readtable("methylation_data.csv") %read in just your methylation values for each sample (ranging from 0 to 1) - exclude gene column
meth(:,1) = []
meth = table2array(meth)


genelist=readtable("genelist.csv") %read in the gene column from your methylation dataset
genelist(:,1) = []
genelist = table2array(genelist)
genelist = string(genelist)

mlist=readtable("mlist.csv") %read in the Mesenchymal gene signatures
mlist = table2array(mlist)
mlist = string(mlist)

elist=readtable("elist.csv") %read in the Epithelial gene signatures
elist = table2array(elist)
elist = string(elist)

%% Calculate EMT score 
for k=1:76 %calculate for all your samples - change k accordingly
[x1,y1,x2,y2]=getEMT(meth(:,k),genelist,elist,mlist);
%This will be your AUC curve output:x1,x2,y1,y2
auc0_1=trapz(y1,x1);auc0_2=trapz(y2,x2);
EMTscoresum(k)=auc0_2-auc0_1; %final score 
end
