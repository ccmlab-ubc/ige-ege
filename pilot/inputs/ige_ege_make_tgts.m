% script to make .tgt files
clear all
close all
clc

namefile = 'ige_ege';

nBsl = 70;
pert_start = nBsl + 1;
nPert = 1200;
total_trials = nBsl + nPert;

tgt_diam = 6; %in mm

% Create perturbation schedule
pert = [-4; -2; 0; 2; 4]; % possible rotation sizes
pert = repmat(pert, 120, 1); 
pert = pert(randperm(600));

rotation_size = zeros(total_trials, 1);
rotation_size(pert_start:2:end) = pert;

vis_fb = zeros(600, 1);
vis_fb(randperm(600, 300)) = 1;
online_fb = ones(total_trials, 1);
online_fb(pert_start + 1:2:end) = vis_fb;
endpoint_fb = online_fb;

for subnum = 1
        
    % Constant values
    T.trialnum = (1:total_trials)';
    T.tgt_dist = ones(total_trials, 1)*90;
    T.tgt_angle = ones(total_trials,1)*90;
    T.rotation = rotation_size;
    T.gain = ones(total_trials, 1);
    T.online_fb = online_fb;
    T.endpoint_fb = endpoint_fb;
    T.clamped_fb = zeros(total_trials, 1);
    T.between_blocks = zeros(total_trials, 1);
    T.tgtsize = ones(total_trials,1)*tgt_diam;
    
    % decide here where you want to place between block messages
    T.between_blocks([nBsl, nBsl+100, nBsl+200, nBsl+300, nBsl+400, nBsl+500,...
        nBsl+600, nBsl+700, nBsl+800, nBsl+900, nBsl+1000, nBsl+1100]) = 2;
    T.between_blocks([total_trials], 1) = 5;
    
    dummy = struct2dataset(T); %%% These two lines are to convert structure into double
    set = double(dummy);
    
    set(:,1) = 1:size(set,1);
    %Add in last Betweenblocks
    %set(end,15) = 1;
    %Variables header file
    header = {'trialnum','tgt_distance','tgt_angle','rotation',...
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

