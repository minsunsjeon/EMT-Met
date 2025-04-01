%% %% Demo: EMT Score Calculation

% Read in methylation matrix (values only - ranging from 0 to 1, excluding gene names)
meth = readtable("methylation_data.csv")
meth(:,1) = []
meth = table2array(meth)% Convert to numeric matrix

% Read in gene names (corresponding to meth row order)
genelist=readtable("genelist.csv") 
genelist(:,1) = []
genelist = table2array(genelist)
genelist = string(genelist)

%Read in the Mesenchymal gene signature list
mlist=readtable("mlist.csv") 
mlist = table2array(mlist)
mlist = string(mlist)

%Read in the Epithelial gene signature list
elist=readtable("elist.csv") 
elist = table2array(elist)
elist = string(elist)

%% Calculate EMT score 
n = size(meth, 2); % number of samples
EMTscoresum = zeros(1, n); % preallocate

for k=1:n 
[x1,y1,x2,y2]=getEMT(meth(:,k),genelist,elist,mlist);
%This will be your AUC curve output:x1,x2,y1,y2
auc0_1=trapz(y1,x1);auc0_2=trapz(y2,x2);
EMTscoresum(k)=auc0_2-auc0_1; %final score 
end
