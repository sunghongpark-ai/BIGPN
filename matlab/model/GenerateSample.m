function File = GenerateSample(Folder, Seed, NumSubject)

Root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 1 || isempty(Folder)
    Folder = fullfile(Root, 'dataset');
end
if nargin < 2
    Seed = 20261006;
end
if nargin < 3
    NumSubject = 906;
end
validateattributes(Seed, {'numeric'}, {'real','scalar','finite','integer','nonnegative','<',2^32});
validateattributes(NumSubject, {'numeric'}, {'real','scalar','finite','integer','positive'});
Counts = [113,98,78,57,18];
Protein = arrayfun(@(j) sprintf('SYN_P%03d',j),1:Counts(1),'UniformOutput',false);
Levels = cell(1,numel(Counts));
Levels{1} = Protein;
for Level = 1:numel(Counts)-1
    Levels{Level+1} = arrayfun(@(j) sprintf('SYN_L%d_%03d',Level,j),1:Counts(Level+1),'UniformOutput',false);
end
Targets = {'abt','gfa','nfl','tau'};
Stream = RandStream('mt19937ar', 'Seed', Seed);
Labels = double(rand(Stream, numel(Targets), NumSubject) < 0.5);
Features = 0.05 + 0.9 * rand(Stream, numel(Protein), NumSubject);
Pairs = nchoosek(1:numel(Protein),2);
[~, Order] = sort(rand(Stream,size(Pairs,1),1));
Selected = sortrows(Pairs(Order(1:1255),:));
Scores = 0.15 + 0.849 * rand(Stream, size(Selected,1), 1);
Links = cell(0,3);
ParentIndex = 0;
for Level = 1:numel(Levels)-1
    Candidates = Protein;
    Depth = zeros(1,numel(Protein));
    if Level > 1
        Candidates = [Candidates, Levels{Level}];
        Depth = [Depth, repmat(Level-1,1,numel(Levels{Level}))];
    end
    for Parent = Levels{Level+1}
        Count = 3 + double(ParentIndex < 121);
        [~, Order] = sort(rand(Stream,numel(Candidates),1));
        for Index = Order(1:Count)'
            Links(end+1,:) = {Candidates{Index},Depth(Index),Parent{1}};
        end
        ParentIndex = ParentIndex + 1;
    end
end
Present = Links(:,1);
for Level = 0:numel(Levels)-1
    for Node = Levels{Level+1}
        if ~ismember(Node{1},Present)
            Links(end+1,:) = {Node{1},Level,''};
        end
    end
end
Links = sortrows(Links,[2,1,3]);
if exist(Folder, 'dir') ~= 7
    mkdir(Folder);
end
File = fullfile(Folder, 'sample.csv');
Handle = fopen(File, 'w', 'n', 'UTF-8');
if Handle < 0
    error('BIGPN:OutputFile', 'Cannot write %s.', File);
end
Cleanup = onCleanup(@() fclose(Handle));
fprintf(Handle, '%s\n', strjoin([{'ID'}, strcat('Y', Targets), Protein], ','));
for Index = 1:NumSubject
    fprintf(Handle, 'SYN_BIGPN_%04d', Index);
    fprintf(Handle, ',%d', Labels(:, Index));
    fprintf(Handle, ',%.6f', Features(:, Index));
    fprintf(Handle, '\n');
end
clear Cleanup
Handle = fopen(fullfile(Folder, 'network.csv'), 'w', 'n', 'UTF-8');
if Handle < 0
    error('BIGPN:OutputFile', 'Cannot write network.csv.');
end
Cleanup = onCleanup(@() fclose(Handle));
fprintf(Handle, 'Protein1,Protein2,Score\n');
for Index = 1:size(Selected,1)
    fprintf(Handle, '%s,%s,%.6f\n', Protein{Selected(Index,1)}, Protein{Selected(Index,2)}, Scores(Index));
end
clear Cleanup
Handle = fopen(fullfile(Folder,'pathway.csv'), 'w', 'n', 'UTF-8');
if Handle < 0
    error('BIGPN:OutputFile', 'Cannot write pathway.csv.');
end
Cleanup = onCleanup(@() fclose(Handle));
fprintf(Handle,'Node,Level,Parent\n');
for Index = 1:size(Links,1)
    fprintf(Handle,'%s,%d,%s\n',Links{Index,:});
end

end
