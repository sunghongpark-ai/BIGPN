function Data = DatasetRead(CohortFile, NetworkFile, PathwayFile, Targets)

Data.Target = TargetCheck(Targets);
[Data.Node, Data.Depth, Data.Link] = PathwayRead(char(PathwayFile));
Data.Protein = Data.Node(Data.Depth == 0);
Data.Network = NetworkRead(char(NetworkFile), Data.Protein);
[Data.Subject, Data.X, Data.Y] = CohortRead(char(CohortFile), Data.Protein, Data.Target);

end

function Target = TargetCheck(Targets)

if isempty(Targets) || ~(ischar(Targets) || iscellstr(Targets) || isstring(Targets))
    error('BIGPN:InvalidOption', 'Targets must be nonempty text.');
end
Target = reshape(cellstr(Targets), [], 1);
if ~all(cellfun(@isvarname, Target)) || numel(unique(Target)) ~= numel(Target)
    error('BIGPN:InvalidOption', 'Targets must be unique valid MATLAB identifiers.');
end

end

function [Node, Depth, Link] = PathwayRead(File)

[Header, Body] = CsvRead(File);
Column = ColumnIndex(File, Header, {'Node'; 'Level'; 'Parent'});
Name = Body(:, Column(1));
Parent = Body(:, Column(3));
Level = NumberParse(File, Body(:, Column(2)), 'Level');
if isempty(Name)
    error('BIGPN:InvalidDataset', '%s has no rows.', File);
end
if any(cellfun(@isempty, Name))
    error('BIGPN:InvalidDataset', '%s: Node must not be empty (line %d).', File, find(cellfun(@isempty, Name), 1) + 1);
end
if any(Level ~= round(Level) | Level < 0)
    error('BIGPN:InvalidDataset', '%s: Level must be a nonnegative integer (line %d).', File, find(Level ~= round(Level) | Level < 0, 1) + 1);
end

[Unique, ~, Group] = unique(Name);
Lowest = accumarray(Group, Level, [], @min);
Highest = accumarray(Group, Level, [], @max);
if any(Lowest ~= Highest)
    error('BIGPN:InvalidDataset', '%s: nodes declared at more than one level: %s.', File, NameList(Unique(Lowest ~= Highest)));
end
NumLevel = max(Lowest);
if NumLevel < 1 || ~all(ismember(0:NumLevel, Lowest))
    error('BIGPN:InvalidDataset', '%s: levels must run from 0 (proteins) to the top pathway level without gaps.', File);
end

[~, Order] = sortrows([Lowest, (1:numel(Lowest))']);
Node = Unique(Order);
Depth = Lowest(Order);
Position = zeros(numel(Unique), 1);
Position(Order) = 1:numel(Order);

[~, ~, ParentGroup] = unique(Parent);
if size(unique([Group, ParentGroup], 'rows'), 1) < numel(Group)
    error('BIGPN:InvalidDataset', '%s has duplicate Node-Parent rows.', File);
end
HasParent = ~cellfun(@isempty, Parent);
[Known, ParentIndex] = ismember(Parent(HasParent), Unique);
if ~all(Known)
    Missing = Parent(HasParent);
    error('BIGPN:InvalidDataset', '%s: parents not declared as nodes: %s.', File, NameList(unique(Missing(~Known))));
end
Link = [Position(Group(HasParent)), Position(ParentIndex)];
ParentDepth = Depth(Link(:, 2));
ChildDepth = Depth(Link(:, 1));
if any(ParentDepth < 1)
    error('BIGPN:InvalidDataset', '%s: proteins (level 0) cannot be parents: %s.', File, NameList(unique(Node(Link(ParentDepth < 1, 2)))));
end
Skip = ChildDepth ~= 0 & ChildDepth ~= ParentDepth - 1;
if any(Skip)
    error('BIGPN:InvalidDataset', '%s has %d links whose child is neither a protein nor a pathway one level below its parent.', File, nnz(Skip));
end
Childless = setdiff(find(Depth >= 1), Link(:, 2));
if ~isempty(Childless)
    error('BIGPN:InvalidDataset', '%s: pathways without children: %s.', File, NameList(Node(Childless)));
end
Link = sortrows(Link, [2, 1]);

end

function Weight = NetworkRead(File, Protein)

[Header, Body] = CsvRead(File);
Column = ColumnIndex(File, Header, {'Protein1'; 'Protein2'; 'Score'});
[KnownFirst, First] = ismember(Body(:, Column(1)), Protein);
[KnownSecond, Second] = ismember(Body(:, Column(2)), Protein);
if ~all(KnownFirst & KnownSecond)
    error('BIGPN:InvalidDataset', '%s: proteins absent from the pathway file: %s.', File, NameList(unique([Body(~KnownFirst, Column(1)); Body(~KnownSecond, Column(2))])));
end
Score = NumberParse(File, Body(:, Column(3)), 'Score');
if any(Score < 0)
    error('BIGPN:InvalidDataset', '%s: Score must be nonnegative (line %d).', File, find(Score < 0, 1) + 1);
end
if any(First == Second)
    error('BIGPN:InvalidDataset', '%s: self-interactions are not allowed (line %d).', File, find(First == Second, 1) + 1);
end
Pair = sort([First, Second], 2);
if size(unique(Pair, 'rows'), 1) < size(Pair, 1)
    error('BIGPN:InvalidDataset', '%s lists a protein pair more than once.', File);
end
NumGene = numel(Protein);
Weight = zeros(NumGene);
Weight(sub2ind([NumGene, NumGene], Pair(:, 1), Pair(:, 2))) = Score;
Weight = Weight + Weight';

end

function [Subject, X, Y] = CohortRead(File, Protein, Target)

[Header, Body] = CsvRead(File);
Label = strcat('Y', Target);
Expected = [{'ID'}; Label; Protein];
if numel(unique(Expected)) < numel(Expected)
    error('BIGPN:InvalidOption', 'The column names ID, %s and the protein names must be distinct.', strjoin(Label', ', '));
end
Column = ColumnIndex(File, Header, Expected);
Subject = Body(:, Column(1));
if isempty(Subject)
    error('BIGPN:InvalidDataset', '%s has no subjects.', File);
end
if any(cellfun(@isempty, Subject))
    error('BIGPN:InvalidDataset', '%s: ID must not be empty (line %d).', File, find(cellfun(@isempty, Subject), 1) + 1);
end
[UniqueSubject, First] = unique(Subject);
if numel(UniqueSubject) < numel(Subject)
    Repeat = Subject(setdiff((1:numel(Subject))', First));
    error('BIGPN:InvalidDataset', '%s: duplicate IDs: %s.', File, NameList(unique(Repeat)));
end
NumTarget = numel(Target);
Y = NumberParse(File, Body(:, Column(2:NumTarget + 1)), 'label')';
if any(Y(:) ~= 0 & Y(:) ~= 1)
    [Row, Line] = find(Y ~= 0 & Y ~= 1, 1);
    error('BIGPN:InvalidDataset', '%s: %s must be 0 or 1 (line %d).', File, Label{Row}, Line + 1);
end
X = NumberParse(File, Body(:, Column(NumTarget + 2:end)), 'protein')';

end

function Column = ColumnIndex(File, Header, Expected)

[UniqueHeader, First] = unique(Header);
if numel(UniqueHeader) < numel(Header)
    error('BIGPN:InvalidDataset', '%s: duplicate column names: %s.', File, NameList(unique(Header(setdiff(1:numel(Header), First)))));
end
[Found, Column] = ismember(Expected, Header);
if ~all(Found)
    error('BIGPN:InvalidDataset', '%s: missing columns: %s.', File, NameList(Expected(~Found)));
end
Extra = setdiff(Header, Expected);
if ~isempty(Extra)
    error('BIGPN:InvalidDataset', '%s: unexpected columns: %s.', File, NameList(Extra));
end
Column = reshape(Column, 1, []);

end

function [Header, Body] = CsvRead(File)

Handle = fopen(File, 'r', 'n', 'UTF-8');
if Handle < 0
    error('BIGPN:InvalidDataset', 'Cannot open %s.', File);
end
Cleanup = onCleanup(@() fclose(Handle));
Text = fread(Handle, [1, Inf], '*char');
if strncmp(Text, char(65279), 1)
    Text = Text(2:end);
elseif strncmp(Text, char([239, 187, 191]), 3)
    Text = Text(4:end);
end
Line = regexp(Text, '\r\n|\n|\r', 'split');
Last = find(~cellfun(@isempty, Line), 1, 'last');
if isempty(Last)
    error('BIGPN:InvalidDataset', '%s is empty.', File);
end
Line = Line(1:Last);
Blank = find(cellfun(@isempty, Line), 1);
if ~isempty(Blank)
    error('BIGPN:InvalidDataset', '%s: line %d is empty.', File, Blank);
end
Field = cell(numel(Line), 1);
for Number = 1:numel(Line)
    Field{Number} = LineSplit(File, Number, Line{Number});
end
Width = cellfun(@numel, Field);
Uneven = find(Width ~= Width(1), 1);
if ~isempty(Uneven)
    error('BIGPN:InvalidDataset', '%s: line %d has %d fields but the header has %d.', File, Uneven, Width(Uneven), Width(1));
end
Cell = regexprep(vertcat(Field{:}), '^[ \t]+|[ \t]+$', '');
Header = Cell(1, :);
Body = Cell(2:end, :);
if any(cellfun(@isempty, Header))
    error('BIGPN:InvalidDataset', '%s has an empty column name.', File);
end

end

function Field = LineSplit(File, Number, Line)

if ~any(Line == '"')
    Field = regexp(Line, ',', 'split');
    return
end
Field = cell(1, nnz(Line == ',') + 1);
Quote = find(Line == '"');
Length = numel(Line);
Position = 1;
Count = 0;
while true
    Count = Count + 1;
    if Position <= Length && Line(Position) == '"'
        Close = Position;
        Closed = false;
        while ~Closed
            Next = Quote(find(Quote > Close, 1));
            if isempty(Next)
                error('BIGPN:InvalidDataset', '%s: line %d has an unterminated quoted field.', File, Number);
            end
            Closed = Next == Length || Line(Next + 1) ~= '"';
            Close = Next + ~Closed;
        end
        Field{Count} = strrep(Line(Position + 1:Close - 1), '""', '"');
        Position = Close + 1;
        if Position <= Length && Line(Position) ~= ','
            error('BIGPN:InvalidDataset', '%s: line %d has text after a closing quote.', File, Number);
        end
    else
        Stop = find(Line(Position:end) == ',', 1);
        if isempty(Stop)
            Stop = Length - Position + 2;
        end
        Field{Count} = Line(Position:Position + Stop - 2);
        Position = Position + Stop - 1;
    end
    if Position > Length
        break
    end
    Position = Position + 1;
end
Field = Field(1:Count);

end

function Value = NumberParse(File, Cell, What)

if isempty(Cell)
    Value = zeros(size(Cell));
    return
end
Valid = ~cellfun(@isempty, regexp(Cell, '^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$', 'once'));
Value = str2double(Cell);
Bad = ~Valid | ~isfinite(Value);
if any(Bad(:))
    [Row, Col] = find(Bad, 1);
    error('BIGPN:InvalidDataset', '%s: %s value ''%s'' on line %d is not a finite number.', File, What, Cell{Row, Col}, Row + 1);
end

end

function Text = NameList(Name)

Name = reshape(cellstr(Name), 1, []);
Text = strjoin(Name(1:min(end, 5)), ', ');
if numel(Name) > 5
    Text = sprintf('%s and %d more', Text, numel(Name) - 5);
end

end
