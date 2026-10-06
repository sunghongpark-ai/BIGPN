function Loss = LossCalculation(Logit, Label)

Loss = sum(mean(max(Logit, 0) - Label .* Logit + log1p(exp(-abs(Logit))), 2));

end
