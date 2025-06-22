function model = DataIndexing(model)

model.FoldTrain = find(contains(string([1:model.NumFold]'),string([model.FoldTest,model.FoldValid]))==0);
model.IdxTrain = find(contains(string(model.CVdata(model.IdxIter,:)),string(model.FoldTrain))==1);
model.IdxValid = find(model.CVdata(model.IdxIter,:)==model.FoldValid);
model.IdxTest = find(model.CVdata(model.IdxIter,:)==model.FoldTest);

model.NumTrain = length(model.IdxTrain);
model.NumValid = length(model.IdxValid);
model.NumTest = length(model.IdxTest);