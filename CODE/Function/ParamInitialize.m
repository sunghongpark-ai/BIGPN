function model = ParamInitialize(model)

%model = ModelTrain;

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',model.IdxIter))

sizeUnet = [model.NumGene,1]; model.Unet = ones(sizeUnet);

sizeWpath1 = size(model.Mpath1); model.Wpath1 = (2*rand(sizeWpath1)-1)*sqrt(6/sum(sizeWpath1));
sizeWpath2 = size(model.Mpath2); model.Wpath2 = (2*rand(sizeWpath2)-1)*sqrt(6/sum(sizeWpath2));
sizeWpath3 = size(model.Mpath3); model.Wpath3 = (2*rand(sizeWpath3)-1)*sqrt(6/sum(sizeWpath3));
sizeWpath4 = size(model.Mpath4); model.Wpath4 = (2*rand(sizeWpath4)-1)*sqrt(6/sum(sizeWpath4));
model.Wpath = [model.Wpath1;model.Wpath2;model.Wpath3;model.Wpath4];

sizeBabt = [model.NumPath,1]; model.Babt = (2*rand(sizeBabt)-1)*sqrt(6/sum(sizeBabt));
sizeBgfa = [model.NumPath,1]; model.Bgfa = (2*rand(sizeBgfa)-1)*sqrt(6/sum(sizeBgfa));
sizeBnfl = [model.NumPath,1]; model.Bnfl = (2*rand(sizeBnfl)-1)*sqrt(6/sum(sizeBnfl));
sizeBtau = [model.NumPath,1]; model.Btau = (2*rand(sizeBtau)-1)*sqrt(6/sum(sizeBtau));

model.sizeParam = [sizeUnet;sizeWpath1;sizeWpath2;sizeWpath3;sizeWpath4;sizeBabt;sizeBgfa;sizeBnfl;sizeBtau];
model.NumParam = prod(sizeUnet)+prod(sizeWpath1)+prod(sizeWpath2)+prod(sizeWpath3)+prod(sizeWpath4)+prod(sizeBabt)+prod(sizeBgfa)+prod(sizeBnfl)+prod(sizeBtau);
model.WeightParam = [model.Unet(:);model.Wpath1(:);model.Wpath2(:);model.Wpath3(:);model.Wpath4(:);model.Babt(:);model.Bgfa(:);model.Bnfl(:);model.Btau(:)];
model.WeightEpoch = cell(model.MaxEpoch,1); model.AdamParam = AdamInitialize(model.NumParam,model.LearnRate);

model.LossTrain = zeros(model.MaxEpoch,1); model.LossValid = zeros(model.MaxEpoch,1);