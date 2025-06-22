function param = ParameterUpdate(param)

param.AdamParam.t = param.AdamParam.t + 1;
param.AdamParam.m = param.AdamParam.beta1*param.AdamParam.m + (1-param.AdamParam.beta1)*param.Gradient;
param.AdamParam.v = param.AdamParam.beta2*param.AdamParam.v + (1-param.AdamParam.beta2)*(param.Gradient.^2);
param.AdamParam.m_hat = param.AdamParam.m / (1 - param.AdamParam.beta1^param.AdamParam.t);
param.AdamParam.v_hat = param.AdamParam.v / (1 - param.AdamParam.beta2^param.AdamParam.t);
param.WeightParam = param.WeightParam - param.AdamParam.alpha * param.AdamParam.m_hat ./ (sqrt(param.AdamParam.v_hat) + param.AdamParam.epsilon);