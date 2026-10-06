function Gradient = BackwardPropagate(Model, Param, Cache, Label, Weight)

Delta = (Cache.Probability - Label) / size(Label, 2);
Gradient = 2 * Model.RegGamma * Weight;
Legacy = strcmp(Model.Gradient, 'legacy');

GradB = zeros(Model.NumPath, Model.NumTarget);
GradB(Model.GeneIndex, :) = Cache.Gene * Delta';
ErrorGene = Param.B(Model.GeneIndex, :) * Delta;
ErrorPath = cell(Model.NumLevel, 1);
for Depth = 1:Model.NumLevel
    Row = Model.Level(Depth).Index;
    GradB(Row, :) = Cache.Path{Depth} * Delta';
    ErrorPath{Depth} = Param.B(Row, :) * Delta;
end

for Depth = Model.NumLevel:-1:1
    Index = Model.ParamIndex.W{Depth};
    Outer = ErrorPath{Depth} * Cache.Source{Depth}';
    Gradient(Index) = Gradient(Index) + Outer(Model.Level(Depth).Compact) .* (Param.S{Depth} .* (1 - Param.S{Depth}));
    Back = Param.A{Depth}' * ErrorPath{Depth};
    if Depth == 1 || ~Legacy
        ErrorGene = ErrorGene + Back(1:Model.NumGene, :);
    end
    if Depth > 1
        ErrorPath{Depth - 1} = ErrorPath{Depth - 1} + Back(Model.NumGene + 1:end, :);
    end
end

if Legacy
    Adjoint = TransposedSolve(Param, Cache.Input * ErrorGene');
    GradU = diag(Adjoint) - diag(TransposedSolve(Param, Adjoint)) .* Param.U;
else
    GradU = sum(TransposedSolve(Param, ErrorGene) .* (Cache.Input - Cache.Gene), 2);
end
Gradient(Model.ParamIndex.U) = Gradient(Model.ParamIndex.U) + GradU;
Gradient(Model.ParamIndex.B) = Gradient(Model.ParamIndex.B) + GradB(:);

end

function Solution = TransposedSolve(Param, Right)

Solution = zeros(size(Right));
Solution(Param.Order, :) = linsolve(Param.Lower, linsolve(Param.Upper, Right, struct('UT', true, 'TRANSA', true)), struct('LT', true, 'TRANSA', true));

end
