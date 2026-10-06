function [Weight, Adam] = ParameterUpdate(Weight, Gradient, Adam)

Adam.Step = Adam.Step + 1;
Adam.Moment1 = Adam.Beta1 * Adam.Moment1 + (1 - Adam.Beta1) * Gradient;
Adam.Moment2 = Adam.Beta2 * Adam.Moment2 + (1 - Adam.Beta2) * (Gradient .^ 2);
Moment1Hat = Adam.Moment1 / (1 - Adam.Beta1 ^ Adam.Step);
Moment2Hat = Adam.Moment2 / (1 - Adam.Beta2 ^ Adam.Step);
Weight = Weight - Adam.Alpha * Moment1Hat ./ (sqrt(Moment2Hat) + Adam.Epsilon);

end
