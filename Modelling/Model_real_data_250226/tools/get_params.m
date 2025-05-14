function [params] = get_params( modelID)
% Lookup table to get number of free parameters per model
% MKW 2018
% SL 2025 allow 3 k/beta

% generates and returns a list of parameter names (params) based on a specified model identifier (modelID). 
% This list of parameters is used later in model fitting procedures.

if contains(modelID, 'one_k')
    kparams = {'k'};
elseif contains(modelID, 'two_k')
    kparams = {'k_self', 'k_other'};
elseif contains(modelID, 'three_k')
    kparams = {'k_self', 'k_otherSame', 'k_otherDifferent'};
elseif ~contains(modelID, 'k')
    kparams = {};
else
    error(['Cant`t determine number of k parameters from model name: ', modelID])
end
    
if contains(modelID, 'one_beta')
    betaparams = {'beta'};
elseif contains(modelID, 'two_beta')
    betaparams = {'beta_self', 'beta_other'};
elseif contains(modelID, 'three_beta')
    betaparams = {'beta_self', 'beta_otherSame', 'beta_otherDifferent'};
else
    error(['Cant`t determine number of beta parameters from model name: ', modelID])
end
    
params = [kparams, betaparams];

end


