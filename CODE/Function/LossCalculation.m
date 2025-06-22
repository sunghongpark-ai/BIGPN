function model = LossCalculation(model)

LabtTrain = (-1/model.NumTrain)*(model.Yabt(model.IdxTrain)*log(model.Pabt(model.IdxTrain)')+(1-model.Yabt(model.IdxTrain))*log(1-model.Pabt(model.IdxTrain)'));
LgfaTrain = (-1/model.NumTrain)*(model.Ygfa(model.IdxTrain)*log(model.Pgfa(model.IdxTrain)')+(1-model.Ygfa(model.IdxTrain))*log(1-model.Pgfa(model.IdxTrain)'));
LnflTrain = (-1/model.NumTrain)*(model.Ynfl(model.IdxTrain)*log(model.Pnfl(model.IdxTrain)')+(1-model.Ynfl(model.IdxTrain))*log(1-model.Pnfl(model.IdxTrain)'));
LtauTrain = (-1/model.NumTrain)*(model.Ytau(model.IdxTrain)*log(model.Ptau(model.IdxTrain)')+(1-model.Ytau(model.IdxTrain))*log(1-model.Ptau(model.IdxTrain)'));

model.LossTrain(model.IdxEpoch) = LabtTrain+LgfaTrain+LnflTrain+LtauTrain;

LabtValid = (-1/model.NumValid)*(model.Yabt(model.IdxValid)*log(model.Pabt(model.IdxValid)')+(1-model.Yabt(model.IdxValid))*log(1-model.Pabt(model.IdxValid)'));
LgfaValid = (-1/model.NumValid)*(model.Ygfa(model.IdxValid)*log(model.Pgfa(model.IdxValid)')+(1-model.Ygfa(model.IdxValid))*log(1-model.Pgfa(model.IdxValid)'));
LnflValid = (-1/model.NumValid)*(model.Ynfl(model.IdxValid)*log(model.Pnfl(model.IdxValid)')+(1-model.Ynfl(model.IdxValid))*log(1-model.Pnfl(model.IdxValid)'));
LtauValid = (-1/model.NumValid)*(model.Ytau(model.IdxValid)*log(model.Ptau(model.IdxValid)')+(1-model.Ytau(model.IdxValid))*log(1-model.Ptau(model.IdxValid)'));

model.LossValid(model.IdxEpoch) = LabtValid+LgfaValid+LnflValid+LtauValid;