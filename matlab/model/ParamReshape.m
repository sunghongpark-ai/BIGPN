function Param = ParamReshape(Model, Weight)

Param.U = Weight(Model.ParamIndex.U);
Param.B = reshape(Weight(Model.ParamIndex.B), Model.NumPath, Model.NumTarget);
Param.W = cell(Model.NumLevel, 1);
Param.S = cell(Model.NumLevel, 1);
Param.A = cell(Model.NumLevel, 1);
for Depth = 1:Model.NumLevel
    Level = Model.Level(Depth);
    Param.W{Depth} = Weight(Model.ParamIndex.W{Depth});
    Param.S{Depth} = 1 ./ (1 + exp(-Param.W{Depth}));
    Param.A{Depth} = sparse(Level.Row, Level.Col, Param.S{Depth}, Level.NumRow, Level.NumCol);
end
Propagation = Model.Laplacian + diag(Param.U);
if ~(rcond(Propagation) >= eps)
    error('BIGPN:SingularPropagation', 'The propagation matrix diag(U) + L is singular to working precision; lower LearnRate.');
end
[Param.Lower, Param.Upper, Param.Order] = lu(Propagation, 'vector');

end
