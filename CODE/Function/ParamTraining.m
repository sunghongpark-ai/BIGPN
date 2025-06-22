function model = ParamTraining(model)

for IdxEpoch = 1:model.MaxEpoch
    model.IdxEpoch = IdxEpoch;
    model.WeightEpoch{model.IdxEpoch,1} = model.WeightParam;
    model = ParamReshape(model);
    model = ForwardPropagate(model);
    model = LossCalculation(model);
    model = BackwardPropagate(model);
    model = ParameterUpdate(model);
end

model.BestEpoch = find(model.LossValid==min(model.LossValid),1);
model.WeightParam = model.WeightEpoch{model.BestEpoch,1};
model = ParamReshape(model); model = ForwardPropagate(model);