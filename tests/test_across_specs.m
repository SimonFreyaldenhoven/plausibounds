%% Preliminaries
clear;
addpath('../lib', '../lib/external/1SimInferenceClass', '../lib/Wald_opt');
rng(03082024)
tic
fig_size = [0 0 200 100];
colors = {'k',...
    [0 0.7 0],...
    'r',...
    'b',...
    [1 0.8 0],...
    [0.8 0 0.8]};

%% Initialize parameter arrays
design_names = {'quadratic', 'wiggly', 'constant', 'no_flat'};
se_values = [0.014, 0.014*0.5^7, 0.014];
corr_values = [0.0, 0.0, 0.95];
corr_names = {'iid', 'large_n', 'strong_corr'};

% Create output directory if it doesn't exist
output_dir = "C:\Users\c1rak01\Documents\projects\WaldBounds\test_files";
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Loop through all combinations
for gg = 1:length(design_names)
    fprintf('Processing design: %s\n', design_names{gg});
    
    % Load true path for this design
    instub = sprintf('../../output/simulation/%s', design_names{gg});
    load(sprintf('%s/true_path.mat', instub), 'delta');
    p = length(delta);
    
    % Create output stub for this design
    outstub = sprintf('../../output/simulation/%s/single', design_names{gg});
    if ~exist(outstub, 'dir')
        mkdir(outstub);
    end
    
    % Set ylims (empty for automatic scaling)
    ylims = [];
    
    for hh = 1:length(corr_values)
        fprintf('  Processing correlation case: %s (%.3f)\n', corr_names{hh}, corr_values(hh));
        
        % Create variance matrix for this correlation and SE combination
        Vhat = create_var(corr_values(hh), se_values(hh), p);
        
        % Generate noisy observation
        dhat = delta + mvnrnd(zeros(p,1), Vhat)';
        
        % Create filenames with design and correlation identifiers
        delta_filename = sprintf('%s/delta_%s_%s.csv', output_dir, design_names{gg}, corr_names{hh});
        vhat_filename = sprintf('%s/vhat_%s_%s.csv', output_dir, design_names{gg}, corr_names{hh});
        dhat_filename = sprintf('%s/dhat_%s_%s.csv', output_dir, design_names{gg}, corr_names{hh});
        
        % Write to CSV files
        writematrix(delta, delta_filename);
        writematrix(Vhat, vhat_filename);
        writematrix(dhat, dhat_filename);
        
        % Save parameters in MAT format as well
        save(sprintf('%s/param_%s.mat', outstub, corr_names{hh}), 'dhat', 'Vhat');
        
        % Generate diagnostics without saving figure
        figname = sprintf('IR_%s_%s', design_names{gg}, corr_names{hh});
        [diagnostics, restricted_LB, restricted_UB, lb, ub] = full_eventplot_l2tf(dhat, Vhat, outstub, figname, ylims, delta, 0); % 0 = don't save figure
        
        % Save diagnostics to CSV
        diagnostics_filename = sprintf('%s/diagnostics_%s_%s.csv', output_dir, design_names{gg}, corr_names{hh});
        save_diagnostics_to_csv(diagnostics, diagnostics_filename);
        
        fprintf('    Saved files: %s, %s, %s, %s\n', delta_filename, vhat_filename, dhat_filename, diagnostics_filename);
    end
end

fprintf('All combinations processed successfully!\n');
toc

%% Helper function to save diagnostics to CSV
function save_diagnostics_to_csv(diagnostics, filename)
    % Extract only diagnostics that do NOT depend on truth or condtruth
    scalar_fields = {};
    scalar_values = {};
    
    % Width metrics (do not depend on truth)
    if isfield(diagnostics, 'pw_width')
        scalar_fields{end+1} = 'pw_width';
        scalar_values{end+1} = diagnostics.pw_width;
    end
    if isfield(diagnostics, 'supt_width')
        scalar_fields{end+1} = 'supt_width';
        scalar_values{end+1} = diagnostics.supt_width;
    end
    if isfield(diagnostics, 'restricted_width')
        scalar_fields{end+1} = 'restricted_width';
        scalar_values{end+1} = diagnostics.restricted_width;
    end
    if isfield(diagnostics, 'wald_width')
        scalar_fields{end+1} = 'wald_width';
        scalar_values{end+1} = diagnostics.wald_width;
    end
    
    % Model parameters (do not depend on truth)
    if isfield(diagnostics, 'df')
        scalar_fields{end+1} = 'df';
        scalar_values{end+1} = diagnostics.df;
    end
    if isfield(diagnostics, 'K')
        scalar_fields{end+1} = 'K';
        scalar_values{end+1} = diagnostics.K;
    end
    if isfield(diagnostics, 'lam1')
        scalar_fields{end+1} = 'lam1';
        scalar_values{end+1} = diagnostics.lam1;
    end
    if isfield(diagnostics, 'lam2')
        scalar_fields{end+1} = 'lam2';
        scalar_values{end+1} = diagnostics.lam2;
    end
    
    % Create table for scalar values
    if ~isempty(scalar_fields)
        scalar_table = table(scalar_fields', cell2mat(scalar_values)', ...
            'VariableNames', {'Field', 'Value'});
        writetable(scalar_table, filename);
    end
    
    % Save vector/matrix fields to separate files (only those not depending on truth)
    [filepath, name, ext] = fileparts(filename);
    
    % Save suptbands if it exists (does not depend on truth)
    if isfield(diagnostics, 'suptbands')
        suptbands_file = fullfile(filepath, [name '_suptbands' ext]);
        writematrix(diagnostics.suptbands, suptbands_file);
    end
    
    % Note: surrogate and surrogate_class may depend on condtruth, so excluded
end