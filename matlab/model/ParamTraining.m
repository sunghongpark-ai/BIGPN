function Fit = ParamTraining(Model, Split, Weight)

Adam = AdamInitialize(Model.NumParam, Model.LearnRate);
LabelTrain = Model.Y(:, Split.IdxTrain);
LabelValid = Model.Y(:, Split.IdxValid);

Fit.LossTrain = nan(Model.MaxEpoch, 1);
Fit.LossValid = nan(Model.MaxEpoch, 1);
Fit.BestEpoch = 0;
Fit.Weight = Weight;
BestLoss = Inf;

for Epoch = 1:Model.MaxEpoch
    Param = ParamReshape(Model, Weight);
    CacheTrain = ForwardPropagate(Model, Param, Split.IdxTrain);
    CacheValid = ForwardPropagate(Model, Param, Split.IdxValid);
    Fit.LossTrain(Epoch) = LossCalculation(CacheTrain.Logit, LabelTrain);
    Fit.LossValid(Epoch) = LossCalculation(CacheValid.Logit, LabelValid);
    if Fit.LossValid(Epoch) < BestLoss
        BestLoss = Fit.LossValid(Epoch);
        Fit.BestEpoch = Epoch;
        Fit.Weight = Weight;
    end
    if Epoch == Model.MaxEpoch
        break
    end
    Gradient = BackwardPropagate(Model, Param, CacheTrain, LabelTrain, Weight);
    [Weight, Adam] = ParameterUpdate(Weight, Gradient, Adam);
    if ~all(isfinite(Weight))
        break
    end
end

if Fit.BestEpoch == 0
    error('BIGPN:TrainingFailed', 'Validation loss was not finite in any epoch of model %d.', Split.IdxModel);
end
Cache = ForwardPropagate(Model, ParamReshape(Model, Fit.Weight), 1:Model.NumSubj);
Fit.Probability = Cache.Probability;

end
