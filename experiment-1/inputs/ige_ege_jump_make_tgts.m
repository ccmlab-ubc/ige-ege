% script to make .tgt files
clear all
close all
clc

namefile = 'ige_ege_jump';

tgt_diam = 6; %in mm
n_baseline = 70;
pert_start = n_baseline + 1;
n_trials_per_level = 100;

% Create perturbation schedule
rot = [-4; -2; 0; 2; 4]; % possible rotation sizes
rot = repmat(rot, n_trials_per_level, 1); 
rot = rot(randperm(length(rot)));

tgt_jump = [-4; -2; 2; 4];
tgt_jump = repmat(tgt_jump, n_trials_per_level, 1);
tgt_jump = tgt_jump(randperm(length(tgt_jump)));

n_training = (length(rot) + length(tgt_jump)) * 2;
total_trials = n_baseline + n_training;

pert_idx = pert_start:2:total_trials; 
pert_idx = pert_idx(randperm(length(pert_idx)));
rot_idx = pert_idx(randperm(length(pert_idx), length(rot)));
jump_idx = pert_idx(~ismember(pert_idx, rot_idx));

rotation_size = zeros(total_trials, 1);
rotation_size(rot_idx) = rot;
jump_size = zeros(total_trials, 1);
jump_size(jump_idx) = tgt_jump;

vis_fb = zeros(n_training / 2, 1);
vis_fb(randperm(length(vis_fb), length(vis_fb) / 2)) = 1;
online_fb = ones(total_trials, 1);
online_fb(pert_start + 1:2:end) = vis_fb;
endpoint_fb = online_fb;

for subnum = 1
    % Constant values
    T.trialnum = (1:total_trials)';
    T.tgt_dist = ones(total_trials, 1)*90;
    T.tgt_angle_1 = ones(total_trials,1)*90;
    T.tgt_angle_2 = T.tgt_angle_1 + jump_size;
    T.rotation = rotation_size;
    T.gain = ones(total_trials, 1);
    T.online_fb = online_fb;
    T.endpoint_fb = endpoint_fb;
    T.clamped_fb = zeros(total_trials, 1);
    T.between_blocks = zeros(total_trials, 1);
    T.tgtsize = ones(total_trials,1)*tgt_diam;
    
    % decide here where you want to place between block messages
    T.between_blocks([n_baseline, n_baseline + (100:100:n_training - 100)]) = 2;
    T.between_blocks(total_trials) = 5;
    
    dummy = struct2dataset(T); %%% These two lines are to convert structure into double
    set = double(dummy);
    
    set(:,1) = 1:size(set,1);
    %Add in last Betweenblocks
    %set(end,15) = 1;
    %Variables header file
    header = {'trialnum','tgt_distance','tgt_angle_1','tgt_angle_2','rotation',...
        'gain','online_fb', 'endpoint_feedback',...
        'clamped_fb','between_blocks','target_size'};
    
    filename = strcat(namefile,num2str(subnum),'.tgt');
    %If you ever need strings in here use this way
    fid = fopen(filename,'wt');
    [rows,cols] = size(set);
    fprintf(fid,'%s\t',header{:});
    for i = 1:rows
        fprintf(fid,'\n');
        for j = 1:cols
            fprintf(fid,'%3.2f\t',set(i,j));
        end
    end
    fclose(fid)   
end

