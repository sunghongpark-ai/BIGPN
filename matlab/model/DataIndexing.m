function Split = DataIndexing(Model, IdxModel)

Split.IdxModel = IdxModel;
Split.IdxIter = Model.CVlist(IdxModel, 1);
Split.FoldTest = Model.CVlist(IdxModel, 2);
Split.FoldValid = Model.CVlist(IdxModel, 3);
Split.FoldTrain = setdiff(1:Model.NumFold, [Split.FoldTest, Split.FoldValid]);

Fold = Model.CVdata(Split.IdxIter, :);
Split.IdxTrain = find(ismember(Fold, Split.FoldTrain));
Split.IdxValid = find(Fold == Split.FoldValid);
Split.IdxTest = find(Fold == Split.FoldTest);

end
