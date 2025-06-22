clc;clear; addpath("./Function")

load('Dataset.mat')

%% Parameter Setting
Parameter.NumIter = 100; Parameter.NumFold = 5;
% User-specific Parameters for Model Training
%Parameter.MaxEpoch = ;
%Parameter.LearnRate = ;
%Parameter.RegGamma = ;

%% Model Initialization
ModelInit = ModelInitialize(Cohort,Network,Pathway,Parameter);

RiskPabt = cell(size(ModelInit.CVlist,1),1);
RiskPgfa = cell(size(ModelInit.CVlist,1),1);
RiskPnfl = cell(size(ModelInit.CVlist,1),1);
RiskPtau = cell(size(ModelInit.CVlist,1),1);

%% Model Parameter Training
parpool(100)
parfor IdxModel = 1:size(ModelInit.CVlist,1)
    ModelTrain = ModelInit;
    ModelTrain.IdxIter = ModelTrain.CVlist(IdxModel,1);
    ModelTrain.FoldTest = ModelTrain.CVlist(IdxModel,2);
    ModelTrain.FoldValid = ModelTrain.CVlist(IdxModel,3);
    ModelTrain = DataIndexing(ModelTrain);
    ModelTrain = ParamInitialize(ModelTrain);
    ModelTrain = ParamTraining(ModelTrain);
    RiskPabt{IdxModel,1} = ModelTrain.Pabt;
    RiskPgfa{IdxModel,1} = ModelTrain.Pgfa;
    RiskPnfl{IdxModel,1} = ModelTrain.Pnfl;
    RiskPtau{IdxModel,1} = ModelTrain.Ptau;
end
poolobj = gcp('nocreate'); delete(poolobj);