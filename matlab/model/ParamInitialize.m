function Weight = ParamInitialize(Model, IdxIter)

Stream = RandStream('mt19937ar', 'Seed', Model.Seed(IdxIter));
Weight = zeros(Model.NumParam, 1);
Weight(Model.ParamIndex.U) = 1;
for Depth = 1:Model.NumLevel
    Level = Model.Level(Depth);
    Draw = rand(Stream, Level.NumRow, Model.NumPath);
    Weight(Model.ParamIndex.W{Depth}) = (2 * Draw(Level.Linear) - 1) * sqrt(6 / (Level.NumRow + Model.NumPath));
end
Draw = rand(Stream, Model.NumPath, Model.NumTarget);
Weight(Model.ParamIndex.B) = (2 * Draw(:) - 1) * sqrt(6 / (Model.NumPath + 1));

end
