function Model = ModelInitialize(Data, Options)

OptionCheck(Options);

Model.Target = Data.Target;
Model.NumTarget = numel(Model.Target);
Model.NumIter = Options.NumIter;
Model.NumFold = Options.NumFold;
Model.MaxEpoch = Options.MaxEpoch;
Model.LearnRate = Options.LearnRate;
Model.RegGamma = Options.RegGamma;
Model.EdgeThreshold = Options.EdgeThreshold;
Model.Gradient = char(Options.Gradient);
Model.Seed = Options.Seed + (0:Options.NumIter - 1)';

Model.Subject = Data.Subject;
Model.Protein = Data.Protein;
Model.X = Data.X;
Model.Y = Data.Y;
[Model.NumGene, Model.NumSubj] = size(Model.X);
Model.Laplacian = NetworkLaplacian(Data.Network, Model.EdgeThreshold);
Model = PathwayStructure(Model, Data.Node, Data.Depth, Data.Link);
Model = ParamLayout(Model);
[Model.CVdata, Model.CVlist] = CrossValidation(Model);

end

function OptionCheck(Options)

if ~isstruct(Options) || ~isscalar(Options)
    error('BIGPN:InvalidOption', 'Options must be a scalar struct.');
end
Required = {'NumIter', 'NumFold', 'MaxEpoch', 'LearnRate', 'RegGamma', 'Seed', 'EdgeThreshold', 'Gradient'};
Missing = Required(~isfield(Options, Required));
if ~isempty(Missing)
    error('BIGPN:MissingOption', 'Specify %s.', strjoin(Missing, ', '));
end
IntegerCheck(Options.NumIter, 'NumIter', 1);
IntegerCheck(Options.NumFold, 'NumFold', 3);
IntegerCheck(Options.MaxEpoch, 'MaxEpoch', 1);
IntegerCheck(Options.Seed, 'Seed', 0);
if Options.Seed + Options.NumIter - 1 > 2^32 - 1
    error('BIGPN:InvalidOption', 'Seed + NumIter - 1 must not exceed 2^32 - 1.');
end
RealCheck(Options.LearnRate, 'LearnRate', false);
RealCheck(Options.RegGamma, 'RegGamma', true);
RealCheck(Options.EdgeThreshold, 'EdgeThreshold', true);
if ~((ischar(Options.Gradient) && isrow(Options.Gradient)) || (isstring(Options.Gradient) && isscalar(Options.Gradient))) || ~any(strcmp(Options.Gradient, {'exact', 'legacy'}))
    error('BIGPN:InvalidOption', 'Gradient must be ''exact'' or ''legacy''.');
end

end

function IntegerCheck(Value, Name, Lower)

Valid = isnumeric(Value) && isscalar(Value) && isreal(Value) && isfinite(Value) && Value == round(Value) && Value >= Lower;
if ~Valid
    error('BIGPN:InvalidOption', '%s must be an integer of at least %d.', Name, Lower);
end

end

function RealCheck(Value, Name, AllowZero)

Valid = isnumeric(Value) && isscalar(Value) && isreal(Value) && isfinite(Value) && (Value > 0 || (AllowZero && Value == 0));
if ~Valid
    Bound = 'positive';
    if AllowZero
        Bound = 'nonnegative';
    end
    error('BIGPN:InvalidOption', '%s must be a finite %s scalar.', Name, Bound);
end

end

function Laplacian = NetworkLaplacian(Weight, Threshold)

NumGene = size(Weight, 1);
Weight(Weight < Threshold) = 0;
Edge = Weight > 0;
if any(Edge(:))
    Value = Weight(Edge);
    Score = zeros(size(Value));
    Spread = std(Value);
    if Spread > 0
        Score = (Value - mean(Value)) / Spread;
    end
    Weight(Edge) = 1 ./ (1 + exp(-Score));
end
Degree = sum(Weight, 2);
Scale = zeros(NumGene, 1);
Scale(Degree > 0) = 1 ./ sqrt(Degree(Degree > 0));
Laplacian = eye(NumGene) - (Scale .* Weight) .* Scale';

end

function Model = PathwayStructure(Model, Node, NodeDepth, Link)

Model.SetPath = Node;
Model.NumPath = numel(Node);
Model.NumLevel = max(NodeDepth);
NodeRow = zeros(Model.NumPath, 1);
Member = cell(Model.NumLevel + 1, 1);
for Depth = 0:Model.NumLevel
    Index = find(NodeDepth == Depth);
    NodeRow(Index) = 1:numel(Index);
    Member{Depth + 1} = Index;
end
Model.GeneIndex = Member{1};

Pair = unique(Link(:, [2, 1]), 'rows');
ParentDepth = NodeDepth(Pair(:, 1));
Level = cell(Model.NumLevel, 1);
for Depth = 1:Model.NumLevel
    NumRow = numel(Member{Depth + 1});
    NumCol = Model.NumGene;
    if Depth > 1
        NumCol = NumCol + numel(Member{Depth});
    end
    Parent = Pair(ParentDepth == Depth, 1);
    Child = Pair(ParentDepth == Depth, 2);
    [Linear, Order] = sort(NodeRow(Parent) + (Child - 1) * NumRow);
    Parent = Parent(Order);
    Child = Child(Order);
    Row = NodeRow(Parent);
    Col = NodeRow(Child) + Model.NumGene * (NodeDepth(Child) > 0);
    Level{Depth} = struct('Index', Member{Depth + 1}, 'NumRow', NumRow, 'NumCol', NumCol, 'Parent', Parent, 'Child', Child, 'Row', Row, 'Col', Col, 'Linear', Linear, 'Compact', Row + (Col - 1) * NumRow);
end
Model.Level = vertcat(Level{:});

end

function Model = ParamLayout(Model)

Count = [Model.NumGene; arrayfun(@(Level) numel(Level.Row), Model.Level); Model.NumPath * Model.NumTarget];
Edge = cumsum([0; Count]);
Block = arrayfun(@(IdxBlock) (Edge(IdxBlock) + 1:Edge(IdxBlock + 1))', (1:numel(Count))', 'UniformOutput', false);
Model.ParamIndex.U = Block{1};
Model.ParamIndex.W = Block(2:end - 1);
Model.ParamIndex.B = Block{end};
Model.NumParam = Edge(end);

end

function [CVdata, CVlist] = CrossValidation(Model)

Pattern = LabelPattern(Model.NumTarget);
[~, Group] = ismember(Model.Y', double(Pattern), 'rows');
Member = arrayfun(@(IdxGroup) find(Group == IdxGroup), (1:size(Pattern, 1))', 'UniformOutput', false);

CVdata = zeros(Model.NumIter, Model.NumSubj);
for IdxIter = 1:Model.NumIter
    Stream = RandStream('mt19937ar', 'Seed', Model.Seed(IdxIter));
    for IdxGroup = 1:numel(Member)
        CVdata(IdxIter, Member{IdxGroup}) = mod(randperm(Stream, numel(Member{IdxGroup})), Model.NumFold) + 1;
    end
end

for Fold = 1:Model.NumFold
    Empty = find(~any(CVdata == Fold, 2), 1);
    if ~isempty(Empty)
        error('BIGPN:InvalidDataset', 'Fold %d of iteration %d is empty; reduce NumFold.', Fold, Empty);
    end
end

[FoldValid, FoldTest, IdxIter] = ndgrid(1:Model.NumFold, 1:Model.NumFold, 1:Model.NumIter);
Keep = FoldValid ~= FoldTest;
CVlist = [IdxIter(Keep), FoldTest(Keep), FoldValid(Keep)];

end

function Pattern = LabelPattern(NumTarget)

Block = cell(NumTarget + 1, 1);
for NumPositive = 0:NumTarget
    PositiveMinority = 2 * NumPositive <= NumTarget;
    NumMinority = min(NumPositive, NumTarget - NumPositive);
    Minority = zeros(1, 0);
    if NumMinority > 0
        Minority = nchoosek(1:NumTarget, NumMinority);
    end
    NumCombination = size(Minority, 1);
    Combination = repmat(~PositiveMinority, NumCombination, NumTarget);
    Combination(sub2ind(size(Combination), repmat((1:NumCombination)', 1, NumMinority), Minority)) = PositiveMinority;
    Block{NumPositive + 1} = Combination;
end
Pattern = vertcat(Block{:});

end
