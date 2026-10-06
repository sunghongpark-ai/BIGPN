function [Risk, Results] = RunBIGPN(Cohort, Options)

arguments
    Cohort {mustBeTextScalar, mustBeFile} = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'dataset', 'sample.csv')
    Options.MaxEpoch (1,1) double = 20
    Options.LearnRate (1,1) double = 0.001
    Options.RegGamma (1,1) double = 0.0001
    Options.Network {mustBeTextScalar, mustBeFile} = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'dataset', 'network.csv')
    Options.Pathway {mustBeTextScalar, mustBeFile} = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'dataset', 'pathway.csv')
    Options.NumIter (1,1) double = 1
    Options.NumFold (1,1) double = 5
    Options.Seed (1,1) double = 1
    Options.Targets {mustBeText} = {'abt', 'gfa', 'nfl', 'tau'}
    Options.EdgeThreshold (1,1) double = 0.15
    Options.Gradient {mustBeMember(Options.Gradient, {'exact', 'legacy'})} = 'exact'
    Options.OutputFile {mustBeTextScalar} = ""
    Options.UseParallel (1,1) logical = false
    Options.Verbose (1,1) logical = false
end


Data = DatasetRead(Cohort, Options.Network, Options.Pathway, Options.Targets);
Model = ModelInitialize(Data, rmfield(Options, {'Network', 'Pathway', 'Targets', 'OutputFile', 'UseParallel', 'Verbose'}));

NumModel = size(Model.CVlist, 1);
Fit = cell(NumModel, 1);
Verbose = Options.Verbose;
NumWorker = 0;
if Options.UseParallel
    NumWorker = Inf;
end

Clock = tic;
parfor (IdxModel = 1:NumModel, NumWorker)
    Split = DataIndexing(Model, IdxModel);
    Result = ParamTraining(Model, Split, ParamInitialize(Model, Split.IdxIter));
    Fit{IdxModel} = Result;
    if Verbose
        fprintf('Model %d/%d (iteration %d, test fold %d, validation fold %d): best epoch %d, validation loss %.6f\n', IdxModel, NumModel, Split.IdxIter, Split.FoldTest, Split.FoldValid, Result.BestEpoch, Result.LossValid(Result.BestEpoch));
    end
end
ElapsedTime = toc(Clock);

Fit = vertcat(Fit{:});
Probability = cat(3, Fit.Probability);
Risk = struct();
for IdxTarget = 1:Model.NumTarget
    Risk.(Model.Target{IdxTarget}) = reshape(Probability(IdxTarget, :, :), Model.NumSubj, NumModel)';
end

Results.Model = Model;
Results.Weight = [Fit.Weight]';
Results.BestEpoch = [Fit.BestEpoch]';
Results.LossTrain = [Fit.LossTrain]';
Results.LossValid = [Fit.LossValid]';
Results.TestRisk = TestRiskSummary(Model, Probability);
Results.ElapsedTime = ElapsedTime;
Results.Options = Options;

if strlength(Options.OutputFile) > 0
    RiskWrite(char(Options.OutputFile), Model, Results.TestRisk);
end

end

function TestRisk = TestRiskSummary(Model, Probability)

Total = zeros(Model.NumTarget, Model.NumSubj, Model.NumIter);
Count = zeros(1, Model.NumSubj, Model.NumIter);
for IdxModel = 1:size(Model.CVlist, 1)
    IdxIter = Model.CVlist(IdxModel, 1);
    Test = Model.CVdata(IdxIter, :) == Model.CVlist(IdxModel, 2);
    Total(:, Test, IdxIter) = Total(:, Test, IdxIter) + Probability(:, Test, IdxModel);
    Count(1, Test, IdxIter) = Count(1, Test, IdxIter) + 1;
end
Average = Total ./ Count;
TestRisk = struct();
for IdxTarget = 1:Model.NumTarget
    TestRisk.(Model.Target{IdxTarget}) = reshape(Average(IdxTarget, :, :), Model.NumSubj, Model.NumIter)';
end

end

function RiskWrite(File, Model, TestRisk)

Value = zeros(Model.NumSubj, Model.NumTarget);
for IdxTarget = 1:Model.NumTarget
    Value(:, IdxTarget) = mean(TestRisk.(Model.Target{IdxTarget}), 1)';
end
Subject = Model.Subject;
Quote = contains(Subject, {',', '"'});
Subject(Quote) = strcat('"', strrep(Subject(Quote), '"', '""'), '"');
Text = [Subject, num2cell(Value)]';
Handle = fopen(File, 'w', 'n', 'UTF-8');
if Handle < 0
    error('BIGPN:OutputFile', 'Cannot write %s.', File);
end
Cleanup = onCleanup(@() fclose(Handle));
fprintf(Handle, '%s\n', strjoin([{'ID'}, strcat('P', Model.Target')], ','));
fprintf(Handle, ['%s', repmat(',%.15g', 1, Model.NumTarget), '\n'], Text{:});

end
