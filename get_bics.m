load('/Users/selmalugtmeijer/Library/CloudStorage/OneDrive-UniversityofBirmingham/Birmingham/Projects/autism/code/Model_real_data_250226/workspaces/EM_fit_results_26-Feb-2025.mat')

M.modid     = {'ms_one_k_one_beta', 'ms_two_k_one_beta', 'ms_three_k_one_beta'...
    'ms_one_k_two_beta', 'ms_two_k_two_beta', 'ms_three_k_two_beta'...
    'ms_one_k_three_beta', 'ms_two_k_three_beta', 'ms_three_k_three_beta'};

num_models = length(M.modid); % Number of models (9)
num_participants = length(allfits{1,1}.(M.modid{1}).fit.bic); % Number of participants (60)

% Initialize a matrix to store BIC values
bic_matrix = NaN(num_participants, num_models);

% Loop through models and extract BIC values
for i = 1:num_models
    model_name = M.modid{i}; % Get the model name
    bic_matrix(:, i) = allfits{1, i}.(model_name).fit.bic; % Extract BIC values
end

data_to_save = [group, bic_matrix];  

column_names = ['Group', M.modid];  

% Convert column names to a string for writing
column_names_str = strjoin(column_names, ',');

% Define the filename
filename = 'bic_data.csv';

% Open file for writing
fid = fopen(filename, 'w');
fprintf(fid, '%s\n', column_names_str);  % Write column names
fclose(fid);

% Write the data
writematrix(data_to_save, filename, 'WriteMode', 'append');

disp('CSV file saved successfully.');

