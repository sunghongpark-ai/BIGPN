function Adam = AdamInitialize(NumVar, Alpha, Beta1, Beta2, Epsilon)

if nargin < 2
    Alpha = 1e-4;
end
if nargin < 3
    Beta1 = 0.9;
end
if nargin < 4
    Beta2 = 0.999;
end
if nargin < 5
    Epsilon = 1e-8;
end

Adam.Alpha = Alpha;
Adam.Beta1 = Beta1;
Adam.Beta2 = Beta2;
Adam.Epsilon = Epsilon;
Adam.Step = 0;
Adam.Moment1 = zeros(NumVar, 1);
Adam.Moment2 = zeros(NumVar, 1);

end
