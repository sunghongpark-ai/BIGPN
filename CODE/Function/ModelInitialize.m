function model = ModelInitialize(cohort,network,pathway,param)

model = cohort;

%% PPI Normalization

Wnet = network.Wnet; Wnet = Wnet.*(Wnet>=0.15);
Wnet(Wnet>0) = 1./(1+exp((-1)*(Wnet(Wnet>0)-mean(Wnet(Wnet>0)))./std(Wnet(Wnet>0))));
diagW = 1./sqrt(sum(Wnet,2)); diagW(isinf(diagW)==1) = 0;
Dnet = diag(diagW); Lnet = eye(model.NumGene)-(Dnet*Wnet*Dnet);
model.Wnet = Wnet; model.Lnet = Lnet;

%% Pathway Structure

model.IdxStruct = pathway.IdxStruct; model.SetStruct = pathway.SetStruct;
model.SetPath = pathway.SetPath; model.NumPath = pathway.NumPath;
model.SetPath0 = pathway.SetPath0; model.NumPath0 = pathway.NumPath0;
model.SetPath1 = pathway.SetPath1; model.NumPath1 = pathway.NumPath1;
model.SetPath2 = pathway.SetPath2; model.NumPath2 = pathway.NumPath2;
model.SetPath3 = pathway.SetPath3; model.NumPath3 = pathway.NumPath3;
model.SetPath4 = pathway.SetPath4; model.NumPath4 = pathway.NumPath4;

model.IdxPath0 = zeros(model.NumPath0,1); for idx = 1:model.NumPath0  model.IdxPath0(idx) = find(model.SetPath0(idx)==model.SetPath);   end
model.IdxPath1 = zeros(model.NumPath1,1); for idx = 1:model.NumPath1  model.IdxPath1(idx) = find(model.SetPath1(idx)==model.SetPath);   end
model.IdxPath2 = zeros(model.NumPath2,1); for idx = 1:model.NumPath2  model.IdxPath2(idx) = find(model.SetPath2(idx)==model.SetPath);   end
model.IdxPath3 = zeros(model.NumPath3,1); for idx = 1:model.NumPath3  model.IdxPath3(idx) = find(model.SetPath3(idx)==model.SetPath);   end
model.IdxPath4 = zeros(model.NumPath4,1); for idx = 1:model.NumPath4  model.IdxPath4(idx) = find(model.SetPath4(idx)==model.SetPath);   end

LeafPath1 = cell(model.NumPath1,1); for idx = 1:model.NumPath1   LeafPath1{idx} = find(ismember(model.SetPath,model.SetStruct(find(model.SetPath1(idx)==model.SetStruct(:,2)),1))==1);  end
LeafPath2 = cell(model.NumPath2,1); for idx = 1:model.NumPath2   LeafPath2{idx} = find(ismember(model.SetPath,model.SetStruct(find(model.SetPath2(idx)==model.SetStruct(:,2)),1))==1);  end
LeafPath3 = cell(model.NumPath3,1); for idx = 1:model.NumPath3   LeafPath3{idx} = find(ismember(model.SetPath,model.SetStruct(find(model.SetPath3(idx)==model.SetStruct(:,2)),1))==1);  end
LeafPath4 = cell(model.NumPath4,1); for idx = 1:model.NumPath4   LeafPath4{idx} = find(ismember(model.SetPath,model.SetStruct(find(model.SetPath4(idx)==model.SetStruct(:,2)),1))==1);  end

model.Mpath1 = zeros(model.NumPath1,model.NumPath); for idx = 1:model.NumPath1  model.Mpath1(idx,LeafPath1{idx}) = 1; end
model.Mpath2 = zeros(model.NumPath2,model.NumPath); for idx = 1:model.NumPath2  model.Mpath2(idx,LeafPath2{idx}) = 1; end
model.Mpath3 = zeros(model.NumPath3,model.NumPath); for idx = 1:model.NumPath3  model.Mpath3(idx,LeafPath3{idx}) = 1; end
model.Mpath4 = zeros(model.NumPath4,model.NumPath); for idx = 1:model.NumPath4  model.Mpath4(idx,LeafPath4{idx}) = 1; end
model.Mpath = [model.Mpath1;model.Mpath2;model.Mpath3;model.Mpath4];

%% CVindexSetting

model.NumIter = param.NumIter; model.NumFold = param.NumFold;

CVdata = zeros(model.NumIter,model.NumSubj);

check = sum([model.Yabt;model.Ygfa;model.Ynfl;model.Ytau]);

index0 = find(check==0);

index1abt = find((check==1)&(model.Yabt==1));
index1gfa = find((check==1)&(model.Ygfa==1));
index1nfl = find((check==1)&(model.Ynfl==1));
index1tau = find((check==1)&(model.Ytau==1));
index1 = unique([index1abt,index1gfa,index1nfl,index1tau]);

index2abtgfa = find((check==2)&(model.Yabt==1)&(model.Ygfa==1));
index2abtnfl = find((check==2)&(model.Yabt==1)&(model.Ynfl==1));
index2abttau = find((check==2)&(model.Yabt==1)&(model.Ytau==1));
index2gfanfl = find((check==2)&(model.Ygfa==1)&(model.Ynfl==1));
index2gfatau = find((check==2)&(model.Ygfa==1)&(model.Ytau==1));
index2nfltau = find((check==2)&(model.Ynfl==1)&(model.Ytau==1));
index2 = unique([index2abtgfa,index2abtnfl,index2abttau,index2gfanfl,index2gfatau,index2nfltau]);

index3abt = find((check==3)&(model.Yabt==0));
index3gfa = find((check==3)&(model.Ygfa==0));
index3nfl = find((check==3)&(model.Ynfl==0));
index3tau = find((check==3)&(model.Ytau==0));
index3 = unique([index3abt,index3gfa,index3nfl,index3tau]);

index4 = find(check==4);

index = unique([index0,index1,index2,index3,index4]);

if (model.NumSubj==length(index))

    for iter = 1:model.NumIter
        
        RandStream.setGlobalStream(RandStream('mt19937ar','Seed',iter))

        CVdata(iter,index0) = mod(randperm(length(index0)),model.NumFold)+1;
        CVdata(iter,index1abt) = mod(randperm(length(index1abt)),model.NumFold)+1;
        CVdata(iter,index1gfa) = mod(randperm(length(index1gfa)),model.NumFold)+1;
        CVdata(iter,index1nfl) = mod(randperm(length(index1nfl)),model.NumFold)+1;
        CVdata(iter,index1tau) = mod(randperm(length(index1tau)),model.NumFold)+1;
        CVdata(iter,index2abtgfa) = mod(randperm(length(index2abtgfa)),model.NumFold)+1;
        CVdata(iter,index2abtnfl) = mod(randperm(length(index2abtnfl)),model.NumFold)+1;
        CVdata(iter,index2abttau) = mod(randperm(length(index2abttau)),model.NumFold)+1;
        CVdata(iter,index2gfanfl) = mod(randperm(length(index2gfanfl)),model.NumFold)+1;
        CVdata(iter,index2gfatau) = mod(randperm(length(index2gfatau)),model.NumFold)+1;
        CVdata(iter,index2nfltau) = mod(randperm(length(index2nfltau)),model.NumFold)+1;
        CVdata(iter,index3abt) = mod(randperm(length(index3abt)),model.NumFold)+1;
        CVdata(iter,index3gfa) = mod(randperm(length(index3gfa)),model.NumFold)+1;
        CVdata(iter,index3nfl) = mod(randperm(length(index3nfl)),model.NumFold)+1;
        CVdata(iter,index3tau) = mod(randperm(length(index3tau)),model.NumFold)+1;
        CVdata(iter,index4) = mod(randperm(length(index4)),model.NumFold)+1;

    end

else
    disp("Error Occured in CVindexSetting")
end

CVfold = []; for fold = 1:model.NumFold  list = [1:model.NumFold]'; list(fold) = []; CVfold = [CVfold;[repmat(fold,size(list)),list]];  end
CVlist = []; for iter = 1:model.NumIter  CVlist = [CVlist;[repmat(iter,[size(CVfold,1) 1]),CVfold]];  end

model.CVdata = CVdata; model.CVlist = CVlist;

%% Training Parameters

model.MaxEpoch = param.MaxEpoch; model.LearnRate = param.LearnRate; model.RegGamma = param.RegGamma;