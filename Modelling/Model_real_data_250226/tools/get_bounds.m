function [bounds] = get_bounds(s, modelID, bounds)
% SL 2025: add 3 k/beta

if ~isfield(bounds, 'k') && contains(modelID, 'k')
bounds.k = [0, maxValue(s, 'k', modelID)];
end

if contains(modelID, 'one_k')
    lbk = bounds.k(1);
    ubk = bounds.k(2);
elseif contains(modelID, 'two_k')
    lbk = [bounds.k(1), bounds.k(1)];
    ubk = [bounds.k(2), bounds.k(2)];
elseif contains(modelID, 'three_k')
    lbk = [bounds.k(1), bounds.k(1), bounds.k(1)];
    ubk = [bounds.k(2), bounds.k(2), bounds.k(2)];
else
    error(['Cant`t determine number of k parameters from model name: ', modelID])
end

if contains(modelID, 'one_beta')
    lbbeta = bounds.beta(1);
    ubbeta = bounds.beta(2);
elseif contains(modelID, 'two_beta')
    lbbeta = [bounds.beta(1), bounds.beta(1)];
    ubbeta = [bounds.beta(2), bounds.beta(2)];
elseif contains(modelID, 'three_beta')
    lbbeta = [bounds.beta(1), bounds.beta(1), bounds.beta(1)];
    ubbeta = [bounds.beta(2), bounds.beta(2), bounds.beta(2)];
else
    error(['Cant`t determine number of beta parameters from model name: ', modelID])
end

bounds.lower = [lbk, lbbeta];   %lower bounds on parameters
bounds.upper = [ubk, ubbeta];   %upper bounds on parameters

end