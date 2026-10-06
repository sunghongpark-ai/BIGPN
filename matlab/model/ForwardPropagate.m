function Cache = ForwardPropagate(Model, Param, Index)

Cache.Input = Model.X(:, Index);
Scaled = Param.U .* Cache.Input;
Cache.Gene = linsolve(Param.Upper, linsolve(Param.Lower, Scaled(Param.Order, :), struct('LT', true)), struct('UT', true));
Cache.Source = cell(Model.NumLevel, 1);
Cache.Path = cell(Model.NumLevel, 1);

Logit = Param.B(Model.GeneIndex, :)' * Cache.Gene;
Previous = zeros(0, size(Cache.Gene, 2));
for Depth = 1:Model.NumLevel
    Cache.Source{Depth} = [Cache.Gene; Previous];
    Cache.Path{Depth} = Param.A{Depth} * Cache.Source{Depth};
    Logit = Logit + Param.B(Model.Level(Depth).Index, :)' * Cache.Path{Depth};
    Previous = Cache.Path{Depth};
end

Cache.Logit = Logit;
Cache.Probability = 1 ./ (1 + exp(-Logit));

end
