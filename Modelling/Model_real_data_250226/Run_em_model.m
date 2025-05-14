%%%%%%%%%
%% Modelling for prosocial motivation task using expectation maximisation
%%%%%%%%%

% Fits models using expectation maximisation (em) approach and does model comparison
% Written by Patricia Lockwood, January 2020
% Based on code by MK Wittmann, October 2018
% Edited by Jo Cutler, August 2020
% Edited by Selma Lugtmeijer, Feb 2025  - for neurotypical / autism study -
% 3 k / 3 beta

%%%%%%%%%
% Step 1 - get data in the format of a varible 's' that contains a struct for each persons data
% Step 2 - run this script to fit models  
% Dependencies: tools subfolder containing required functions e.g. fit_PM_model
%               models subfolder containing various comp models you have made 
% Step 3 - compare the AIC's and BIC's using the script visualize_model_PM
% (see below)

%% Input for script
%       - Participants data file format saved in 's':

%% Output from script
%   - workspaces/EM_fit_results_[date] has all variables from script
%           - 's.PM.em' contains model results including the model parameters per ppt
%   - datafiles in specified output directory:
%       - PM_model_fit_statistics.csv - model comparison fit statistics
%       - EM_fit_parameters.csv - estimated parameters for each participant
%       - Compare_fit_between_groups.csv - median R^2 for each participant with group index

%% Prosocial motivation models based Lockwood et al. (2017)
% test different variations of discount rate (k) and beta parameters
% SL: for two k/beta two different other conditions taken together
%   - one_k_one_beta
%   - two_k_one_beta
%   - three_k_one_beta
%   - one_k_two_beta
%   - two_k_two_beta
%   - three_k_two_beta
%   - one_k_three_beta
%   - two_k_three_beta
%   - three_k_three_beta
% and shape of discounting:
%   - parabolic 

%%

cd /Users/selmalugtmeijer/Library/CloudStorage/OneDrive-UniversityofBirmingham/Birmingham/Projects/autism/code/Model_real_data_250226

%== -I) Prepare workspace: ============================================================================================
clear
clearvars
addpath('models');
addpath('tools');
% addpath('spm12');
setFigDefaults; % custom function - make sure it is in the folder

%== 0) Load and organise data: ==========================================================================================
% load data:
file_name = 'data_for_MLE_modelling.mat'; % specify data **
load(file_name); % .mat file saved from the behavioural script that contains all participants data in 's'

% change the organisation of the 2 different participant groups
% Assuming s.PM.group is a 1x60 cell array with 'ASD' or 'NEUROTYPICAL'
s.PM.groups = cellfun(@(x) strcmp(x, 'ASD') + 2 * strcmp(x, 'NEUROTYPICAL'), s.PM.group)';

% recode into same/different as own neurotype other? If not keep neurotypical and
% ASD other
recodeNeurotype = 1; 

if recodeNeurotype == 1
    s = recode(s);
end

s.PM.expname = 'neurotype';
s.PM.em = {};

bounds.beta = [0, 10];

output_dir = '../Results/'; % enter path to save output in **

% how to fit models:
M.dofit     = 1;                                                                            % whether to fit (1) or not (0)
M.doMC      = 1;                                                                            % whether to do model comparison or not
M.modid     = {'ms_one_k_one_beta', 'ms_two_k_one_beta', 'ms_three_k_one_beta'...
    'ms_one_k_two_beta', 'ms_two_k_two_beta', 'ms_three_k_two_beta'...
    'ms_one_k_three_beta', 'ms_two_k_three_beta', 'ms_three_k_three_beta'};

fitMeasures = {'lme','bicint','xp','pseudoR2','choiceProbMedianR2'}; % which fit measures to calculate **
criteria = 'bicint'; % of above, which to use to choose the best model *

% define experiment of interest:
e = 'PM';% **

%== I) RUN MODELS: ==========================================================================================

if M.dofit
    %%% EM fit %%%
    for im = 1:numel(M.modid) % for the number of models, can also fit in parallel with parfor
        if ~isfield(s.(e).em,(M.modid{im}))
            rng default % resets the randomisation seed to ensure results are reproducible (MATLAB 2019b)
            attempt = 0;
            max_attempts = 50;
            dotry=1;
            while dotry && attempt < max_attempts + 1
                try
                    close all;
                    allfits{im} = EMfit_ms_par(s.(e),M.modid{im},bounds,attempt,max_attempts);
                    dotry=0;
                catch
                    attempt = attempt + 1;
                    dotry=1; disp('caught');
                end
            end
            if attempt == max_attempts
                disp(['Model ' num2str(im) ' failed after ' num2str(max_attempts) ' attempts. Moving on...']);
            end
        end
        attempts(im) = attempt;
    end
    
    for im = 1:numel(M.modid) % for the number of models
%         try
        s.(e).em.(M.modid{im}) = allfits{1,im}.(M.modid{im});
%         catch
%         end
    end
    
    save(['workspaces/EM_fit_results_',date,'.mat'])  
    
%     M.modid = fieldnames(s.(e).em);
    
    %%% calc BICint for EM fit
    parfor im = 1:numel(M.modid)
        rng default % resets the randomisation seed to ensure results are reproducible (MATLAB 2019b)
        allbics{im} =  cal_BICint_ms(s.(e),M.modid{im},bounds);
    end
    
    for im = 1:numel(M.modid) % for the number of models
        s.(e).em.(M.modid{im}).fit.bicint = allbics{1,im};
    end
    
end

%== II) COMPARE MODELS: ==========================================================================================

if M.doMC
    rng default % resets the randomisation seed to ensure results are reproducible (MATLAB 2019b)
    s.(e) = EMmc_ms(s.(e),M.modid);
    
    % Calculate R^2 & extract model fit measures
    
    for im = 1:numel(M.modid) % for the number of models
        s.(e).em.(M.modid{im}).fit.pseudoR2 = pseudoR2(s.(e),M.modid{im},2,1);
        s.(e) = choiceProbR2(s.(e),M.modid{im},1);
    end
    [fits.(e),fitstab.(e)] = getfits(s.(e),fitMeasures,M.modid);
end

%== III) LOOK AT PARAMETERS: ==========================================================================================

switch criteria
    case {'xp', 'pseudoR2','choiceProbMedianR2'}
        bestmod = find(fitstab.(e).(criteria) == max(fitstab.(e).(criteria)));
    case {'lme', 'bicint'}
        bestmod = find(fitstab.(e).(criteria) == min(fitstab.(e).(criteria)));
end
bestname = M.modid{bestmod};
disp(['Extracting parameters from ', bestname,' based on best ',criteria])

for i=1:length(s.(e).ID)
    IDs(i, :)=s.(e).ID{1,i}.ID;
    try
        group(i, :)=s.(e).groups(i,1);
    catch
        group = [];
    end
end

params = getparams(s.(e), bestname, bounds, IDs, group);

% SL also save parameters 3k 3beta
disp(['Extracting parameters from ms_three_k_three_beta based on best ',criteria])
params2 = getparams(s.(e), 'ms_three_k_three_beta', bounds, IDs, group); % model preregistration

% SL also save parameters 3k 1beta
disp(['Extracting parameters from ms_three_k_one_beta based on best ',criteria])
params3 = getparams(s.(e), 'ms_three_k_one_beta', bounds, IDs, group); % model preregistration

%== IV) SAVE: ==========================================================================================

writetable(params.all_table,[output_dir,'EM_fit_parameters_',bestname,'.csv'],'WriteRowNames',true) 
writetable(params2.all_table,[output_dir,'EM_fit_parameters_ms_three_k_three_beta.csv'],'WriteRowNames',true) % model preregistration
writetable(params3.all_table,[output_dir,'EM_fit_parameters_ms_three_k_one_beta.csv'],'WriteRowNames',true) % model preregistration

save(['workspaces/EM_fit_results_',date,'.mat'])  

fit = fits.(e);
fit = [[1:numel(M.modid)]',fit];
fit(:,end+1) = fit(:,find(contains(fitMeasures, 'bicint'))+1) - min(fit(:,find(contains(fitMeasures, 'bicint'))+1));
fittabnum = cell2table(num2cell(fit), 'VariableNames', ['model', fitMeasures, 'relbic']);
writetable(fittabnum,[output_dir,e,'_model_fit_statistics','.csv'],'WriteRowNames',true)

%== V) COMPARE PARAMETERS BETWEEN GROUPS: ==========================================================================================

compareFitGroups = [(s.PM.groups),(s.PM.em.(bestname).fit.bic),(s.PM.em.(bestname).fit.eachSubProbMedianR2)];
compareFitTab = cell2table(num2cell(compareFitGroups), 'VariableNames', {'group', 'bic', 'R2'});
writetable(compareFitTab,[output_dir,'Compare_em_fit_between_groups_.csv'],'WriteVariableNames',true) % export file to analyse in R

doVBA = 1;

if doVBA == 1
    
    % optional - use the VBA toolbox - https://mbb-team.github.io/VBA-toolbox/
    % to calcuate the xp and expected frequencies in each group and test
    % whether the groups are different in model fit
    
    L1ind = find(s.PM.groups == 1); % group 1
    L2ind = find(s.PM.groups == 2); % group 2
    Lallind = [L1ind;L2ind]; % all groups
    
    for im = 1:numel(M.modid) % for the number of models
        Lall(im,1:length(Lallind)) = s.PM.em.(M.modid{im}).fit.lme(Lallind); % extract log model evidence for each participant
        L1(im,1:length(L1ind)) = s.PM.em.(M.modid{im}).fit.lme(L1ind);
        L2(im,1:length(L2ind)) = s.PM.em.(M.modid{im}).fit.lme(L2ind);
    end

    [posterior,out] = VBA_groupBMC(Lall); % all groups
    [posterior1, out1] = VBA_groupBMC(L1) ; % group 1
    [posterior2, out2] = VBA_groupBMC(L2) ; % group 2
    [h12, p12] = VBA_groupBMC_btwGroups({L1, L2}); % between groups comparison
    
else
end

% Extract groups and R2 values
group1_R2 = compareFitGroups(compareFitGroups(:,1) == 1, 3); % Select R2 values for Group 1
group2_R2 = compareFitGroups(compareFitGroups(:,1) == 2, 3); % Select R2 values for Group 2

% Compute mean and median for each group
mean_R2 = mean(compareFitGroups(:,3))
median_R2 = median(compareFitGroups(:,3))

mean_R2_group1 = mean(group1_R2)
median_R2_group1 = median(group1_R2)

mean_R2_group2 = mean(group2_R2)
median_R2_group2 = median(group2_R2)