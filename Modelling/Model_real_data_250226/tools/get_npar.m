function [ npar ] = get_npar( modelID)
% Lookup table to get number of free parameters per model
% JC 2022 (from MKW 2018)
% SL 2025 (add option 3 k/beta)

%%%%%
if       contains(modelID,'three_k'),             nk = 3;
elseif   contains(modelID,'two_k'),               nk = 2;
elseif   contains(modelID,'one_k'),               nk = 1;
end

if       contains(modelID,'three_beta'),          nb = 3;
elseif   contains(modelID,'two_beta'),            nb = 2;
elseif   contains(modelID,'one_beta'),            nb = 1;
end

npar = nk + nb;

end

