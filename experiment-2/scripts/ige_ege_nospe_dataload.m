% This function loads in all subject data. If the data were saved in rows
% (eg, x is a vector made up of n rows, 1 column), then you must transpose
% them.

clear all; close all

tic

baseDir= '/Users/hyosubkim/Documents/GitHub/Projects/ige-ege/experiment-2/scripts';
dataDir='/Users/hyosubkim/Library/CloudStorage/OneDrive-UBC/projects/ige-ege/experiment-2/';
analyzeDir='/Users/hyosubkim/Documents/GitHub/Projects/ige-ege/experiment-2/results';

cd(dataDir);


subj = {'s01_ige_ege_nospe', 's02_ige_ege_nospe', 's03_ige_ege_nospe', 's04_ige_ege_nospe',...
    's05_ige_ege_nospe', 's06_ige_ege_nospe', 's07_ige_ege_nospe', 's08_ige_ege_nospe',...
    's09_ige_ege_nospe', 's10_ige_ege_nospe', 's11_ige_ege_nospe', 's12_ige_ege_nospe'};

%%% look for table. if it's there, we will append stuff. look to see how
%%% many subjects are in existing table. 
% if (exist('PD_motorPlan_table.mat') & exist('PD_motorPlan_table.mat')) == 1
%     load('TgtSwitch_table.mat'); load('TgtSwitch_movefile_table.mat')
% 
%     T=T;
%     M=M;
%     newSubj = setdiff(subj,unique(T.id)); 
%     num_tested_subj=max(T.si);
%     if length(newSubj)==0
%         return
%     end
% else
%     T=[];
%     M=[];
%     newSubj = subj;
%     num_tested_subj = 0;
% end

T=[];
M=[];
newSubj = subj;
num_tested_subj = 0;

tpe = 1;    % trials per epoch
mm2pixel = 3.6137;  pixel2mm = 1/mm2pixel; pixel2cm = pixel2mm./10;
dt_tablet = 0.005;  % check to make sure this is correct
startofexpmt = 1;
maxReachTime = 600; %number of samples assuming 500 hz sampling rate

% Define some anonymous helper functions
NaNmask_ = @(x) (double(x)./double(x).*double(x));
Outlier_ = @(x,Niqr) (abs(x-nanmean(nanmedian(x))) > Niqr*nanmean(iqr(x)));

% loop through all subjects
for s = 1:length(newSubj)
    s
    newSubj{s}
    load(fullfile(dataDir,[newSubj{s}, '.mat']));
    
    %check for new var names bc of switch to dual disp var names
    if exist('tgtxSubj')
        tgtx = tgtxSubj;
        tgty = tgtySubj;
        xCenter = xCenterSubj;
        yCenter = yCenterSubj;
    end
    
    nt = max(trial_num);
    S = [];
    
    S.hand_x = nan(nt,maxReachTime);     % game loop ran at 500 hz
    S.hand_y = nan(nt,maxReachTime);
    S.hand_dist = nan(nt,maxReachTime);
    S.hand_v = nan(nt,maxReachTime);
    S.trialtime = nan(nt,maxReachTime);
    S.SN(1:nt,1) = s;
%     S.tgtpos(1:nt,1) = [tgtlocSubj(1,1)-xCenterSubj];
%     S.tgtpos(1:nt,2) = [yCenterSubj-tgtlocSubj(1,2)];
    D = [];
    V = [];
    Z = [];

    Z.move_trial = trial_move;
    Z.gamephase = gamephase_move;
    
    for i=1:nt
        
        %hand_dist = sqrt(hand_x.^2 + hand_y.^2);
        
        % trimming hand data
        idx1 = find(Z.move_trial==i & Z.gamephase==2);
        %idx2 = find(Z.move_trial==i+1 & Z.gamephase==0,100);
        idx = [idx1];
        if length(idx)<=maxReachTime
%             S.hand_x(i,1:length(idx)) = hand_x(idx);
%             S.hand_y(i,1:length(idx)) = hand_y(idx);
            % For flipped monitor:
            S.hand_x(i,1:length(idx)) = hand_x(idx)*(-1);
            S.hand_y(i,1:length(idx)) = hand_y(idx)*(-1);
            S.hand_dist(i,1:length(idx)) = hand_dist_all(idx);
            S.trialtime(i,1:length(idx)) = trial_time(idx);
        else
            S.hand_x(i,:) = hand_x(idx(1:maxReachTime));
            S.hand_y(i,:) = hand_y(idx(1:maxReachTime));
            S.hand_dist(i,:) = hand_dist_all(idx(1:maxReachTime));
            S.trialtime(i,:) = trial_time(idx(1:maxReachTime));
        end
        
    end
    
    [S.hx,S.hy,S.hdist,S.trialt] = deal(NaN(size(S.hand_x)));
    
    for k=1:nt
        
        % resample hand position based on samples where position
        % changes(this is ok because we're looking at data during movement
        % only)
        moveidx1 = abs(diff(S.hand_y(k,:))) > 0;  % creating logical to keep track of when there is movement
        moveidx2 = abs(diff(S.hand_x(k,:))) > 0;
        
        ii = [1, find(moveidx1+moveidx2)+1];    % indices of when there is new s
        Nii(k) = length(ii);
        S.hx(k,1:Nii(k)) = S.hand_x(k,ii)*pixel2mm;
        S.hy(k,1:Nii(k)) = S.hand_y(k,ii)*pixel2mm;
        S.hdist(k,1:Nii(k)) = S.hand_dist(k,ii)*pixel2mm;
        S.trialt(k,1:Nii(k)) = S.trialtime(k,ii);
        
    end
    
    hvx = []; hvy = [];
    hvx = diff(S.hx')'./dt_tablet;    % compute velocity based on a simple difference (this can be smoothed later if necessary)
    hvy = diff(S.hy')'./dt_tablet;
    S.absvel = sqrt(hvx.^2 + hvy.^2);
    S.radvel = diff(S.hdist')'/dt_tablet;
    S.absacc = diff(S.absvel')'/dt_tablet';     % compute acceleration
    
    % Remove small number of huge velocity spikes -- figure this out!
    zz = zeros(1,nt);
    q = Outlier_([zz; diff(S.hdist')],10)';
    S.hdist_SpikeNum = sum(sum(q));
    S.hdist = S.hdist.*NaNmask_(~q);
    
    q = Outlier_([zz; diff(S.absvel')],10)';
    S.absvel_SpikeNum = sum(sum(q));
    S.absvel = S.absvel.*NaNmask_(~q);
    
    q = Outlier_([zz; diff(S.radvel')],10)';
    S.radvel_SpikeNum = sum(sum(q));
    S.radvel = S.radvel.*NaNmask_(~q);
    
    q = Outlier_([zz; diff(S.absacc')],10)';
    S.absacc_SpikeNum = sum(sum(q));
    S.absacc = S.absacc.*NaNmask_(~q);
    
    
    S.xi = S.hx(:,1);   % save some basic info about the movement: initial position, MT, max y-vel
    S.yi = S.hy(:,1);
    S.MT = (Nii*dt_tablet)';    
    [S.absvelmax, S.absvelmax_idx] = max(S.absvel');
    S.absvelmax = S.absvelmax';
    S.absvelmax_idx = S.absvelmax_idx';
    
    [S.radvelmax, S.radvelmax_idx] = max(S.radvel');
    S.radvelmax = S.radvelmax';
    S.radvelmax_idx = S.radvelmax_idx';
    
    
    for k = 1:nt
        
%         S.xf(k,1) = S.hx(k,Nii(k));
%         S.yf(k,1) = S.hy(k,Nii(k));
        S.xf(k,1) = hX_ep(k);
        S.yf(k,1) = hY_ep(k);
        S.hand_max_dist(k,1) = sqrt((hX_ep(k) - xCenterSubj)^2 + (hY_ep(k) - yCenterSubj)^2);
        S.maxv_hand_ang(k,1) = atan2d(S.hy(k,S.absvelmax_idx(k)), S.hx(k,S.absvelmax_idx(k)));
        S.maxradv_hand_ang(k,1) = atan2d(S.hy(k,S.radvelmax_idx(k)), S.hx(k,S.radvelmax_idx(k)));
        
        rot_hand_theta(k) = atan2d(sind(hand_angle(k) - tgt_ang_1(k)), cosd(hand_angle(k) - tgt_ang_1(k)));
        xy_vec = [hX_ep(k) - xCenterSubj, hY_ep(k) - yCenterSubj];
        tgt_angle = tgt_ang_1(k) * pi/180;
        rot_angle = -tgt_angle;
        R = [cos(rot_angle), sin(rot_angle); -sin(rot_angle), cos(rot_angle)];
        rot_xy_vec = R * xy_vec';
        rot_hX(i, 1) = rot_xy_vec(1);
        rot_hY(i, 1) = rot_xy_vec(2);

        % calculate hand angle at peak radial velocity
        theta_maxradv(k) = atan2d(sind(S.maxradv_hand_ang(k) - tgt_ang_1(k)), cosd(S.maxradv_hand_ang(k) - tgt_ang_1(k)));
        
    end
    
%     % rotate target to zero to get hand angle 
%     for i = 1:nt
%         
%         rot_hand_theta(i) = atan2d(sind(hand_angle(i) - tgt_ang(i)), cosd(hand_angle(i) - tgt_ang(i)));
%         
%         xy_vec = [hX_ep(i) - xCenter, hY_ep(i) - yCenter];
%         tgt_angle = tgt_ang(i) * pi/180;
%         rot_angle = -tgt_angle;
%         R = [cos(rot_angle), sin(rot_angle); -sin(rot_angle), cos(rot_angle)];
%         rot_xy_vec = R * xy_vec';
%         rot_hX(i, 1) = rot_xy_vec(1);
%         rot_hY(i, 1) = rot_xy_vec(2);
%         
%     end
    
    movement_cycle = ceil([1:nt]/tpe)';
    
    if isempty(T)==1
        D.si(1:nt,1) = 1;
    else
        D.si(1:nt,1) = max(T.si)+1;
    end
    D.SN(1:nt,1) = s;
    D.id(1:nt,:) = {newSubj{s}(2:5)};
    D.tester(1:nt,1) = newSubj{s}(1);
    D.TN(1:nt,1) = (1:nt);
    D.move_cycle = movement_cycle;
    D.hX = (S.xf-xCenterSubj) * pixel2mm;
    D.hY = (S.yf-yCenterSubj) * pixel2mm;
    D.rot_hX = rot_hX * pixel2mm;
    D.rot_hY = rot_hY * pixel2mm;
    D.rotation = rotation;
    D.tgt_jump = tgt_ang_2 - tgt_ang_1;
    D.hand_max_dist = S.hand_max_dist * pixel2mm;
    D.radvelmax = S.radvelmax * pixel2cm;
    D.tgtX = (tgtx_1_subj-xCenterSubj) * pixel2mm;
    D.tgtY = (tgty_1_subj-yCenterSubj) * pixel2mm;
    D.rot_hand_theta = rot_hand_theta';
    D.theta_maxradv = theta_maxradv';
    D.raw_ep_hand_ang = hand_angle;
    D.tgt_ang = tgt_ang_1;
    D.tgt_dist = tgt_dist * pixel2mm;
    D.fbi = online_fb;
    D.MT = MTs';
    D.RT = RTs';
    D.ST = SearchTimes';
    D.radvelmax = S.radvelmax;
      
    V.hx = S.hx;
    V.hy = S.hy;
    V.absvel = S.absvel;
    V.absacc = S.absacc;
    V.hdist = S.hdist;
    V.radvel = S.radvel;
    
    movement_cycle = ceil([1:length(D.TN)]/tpe)';
    D.TN(1:length(D.TN)) = 1:length(D.TN);
    D.move_cycle = movement_cycle;
            
    temp = struct2table(D);
    tempmove = struct2table(V);
    
    % index starting with trial 33, since first 32 trials were practice
    T = [T;temp(startofexpmt:end,:)];
    M = [M;tempmove(startofexpmt:end,:)];
    
    elapsed_times(s) = elapsedTime;
    
    %need to clear some vars bc of switch to dual display and new var names
    clear tgtx tgty
end

mean_elapsed_experiment_time = mean(elapsed_times)/60;
minutes = floor(mean_elapsed_experiment_time)
seconds = 60*(mean_elapsed_experiment_time  - floor(mean_elapsed_experiment_time) )

cd(analyzeDir)

save('ige_ege_nospe_table.mat','T');
save('ige_ege_nospe_movefile_table.mat','M');
writetable(T, 'ige_ege_nospe.csv', 'Delimiter', ',')


toc