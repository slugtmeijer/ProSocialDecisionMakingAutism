function [fval,fit] = mod_ms_all(behavData, q, fitop, modelID, bounds, varargin)

% runs standard Prosocial motivation disocunt model
% P Lockwood modified 1 July 2019 from MK Wittmann, Oct 2018
% SL 2025: add 3 k/beta

% INPUT:    - behavData: behavioural input file
% OUTPUT:   - fval and fitted variables
%

%%
% -------------------------------------------------------------------------------------
% 1 ) Define free parameters
% -------------------------------------------------------------------------------------

if nargin > 5
    prior      = varargin{1};
end

params = get_params(['ms_', modelID]);

for p = 1:length(params)
    
    qt(p) = norm2positive(q(p), [bounds.lower(p), bounds.upper(p)]); % transform parameters from gaussian space to model space
    
end

all_prob = [];
all_V  = [];

%%% 0.) Load information for that subject:    % load in each subjects variables for the experiment
chosen   = behavData.choice';
effort   = behavData.effort';
reward   = behavData.reward';
agent    = behavData.agent';

% SL create an extra agent variable in which you take the two other groups
% together - so there is only self - other
agent2 = agent;  % Copy the original variable
agent2(agent2 == 3) = 2;  % Recode value 3 into 2


indmiss2=chosen==2;
indmiss=isnan(chosen);
indmiss = (indmiss + indmiss2) > 0;
chosen(indmiss) = [];
effort(indmiss) = [];
reward(indmiss) = [];
agent(indmiss)  = [];
agent2(indmiss)  = [];

% Define free parameters

if contains(modelID, 'one_k')
    %discount = qt(1);
    discount = (repelem(qt(1), length(agent), 1))';
    beta1 = 2;
elseif contains(modelID, 'two_k')
    discount = (agent2==1).*qt(1) + (agent2==2).*qt(2); % SL changed agent -> agent2
    beta1 = 3;
elseif contains(modelID, 'three_k')
    discount = (agent==1).*qt(1) + (agent==2).*qt(2) + (agent==3).*qt(3); % SL changed added 3rd option
    beta1 = 4; % SL change from 3 to 4
elseif ~contains(modelID, 'k')
    discount = [];
    beta1 = 1;
else
    error(['Cant`t determine number of k parameters from model name: ', modelID])
end

if contains(modelID, 'one_beta')
    %beta = qt(beta1);
    beta = (repelem(qt(beta1), length(agent), 1))';
elseif contains(modelID, 'two_beta')
    beta = (agent2==1).*qt(beta1) + (agent2==2).*qt(beta1+1); % SL changed agent -> agent2
elseif contains(modelID, 'three_beta')
    beta = (agent==1).*qt(beta1) + (agent==2).*qt(beta1+1) + (agent==3).*qt(beta1+2); % SL changed added 3rd option
else
    error(['Cant`t determine number of beta parameters from model name: ', modelID])
end

base = 1;
   
%%%% Model - devalue reward by effort
    
val = reward - (discount.*(effort.^2));
    
prob =  exp(val.*beta)./(exp(base*beta) + exp(beta.*val));
prob(~chosen) =  1 - prob(~chosen);
prob = prob(1,:);

% SL try out
epsilon = 1e-16;  % Prevent prob from becoming too small
prob = max(prob, epsilon);
%prob = min(prob, 1 - epsilon);  % Prevent prob from reaching exactly 1

%%% 4. now save stuff:
all_V      =  val;
all_prob   =  prob;


% all choice probablities
ChoiceProb=all_prob';
ChoiceProb(ChoiceProb < 0.0001) = 0.0001; % very small probabilities can cause problems

% -------------------------------------------------------------------------------------
% 4 ) Calculate model fit:
% -------------------------------------------------------------------------------------

nll =-nansum(log(ChoiceProb));                                              % the thing to minimize

if fitop.doprior == 0                                                               % NLL fit
    fval = nll;
elseif fitop.doprior == 1                                                           % EM-fit:   P(Choices | h) * P(h | O) should be maximised, therefore same as minimizing it with negative sign
    fval = -(-nll + prior.logpdf(q));
end

% % make sure f is not just low because of Nans in prob-variable:

sumofnans=sum(sum(isnan(chosen)));
if sum(isnan(ChoiceProb))~=sumofnans
    disp('ERROR NaNs in choice and choice prob dont agree');
    keyboard;
    return;
end

% -------------------------------------------------------------------------------------
% 5) Calculate additional Parameters and save:
% -------------------------------------------------------------------------------------

if fitop.dofit ==1
    
    fit         = struct;
    fit.xnames  = params;
    
    fit.choiceprob = [ChoiceProb];
    fit.mat    = [all_V];
    fit.names  = {'V'};
    
end

end
