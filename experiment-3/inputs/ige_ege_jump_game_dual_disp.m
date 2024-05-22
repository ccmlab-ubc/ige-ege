function ige_ege_jump_game_dual_disp(name_prefix, tgt_file_name_prefix, tgt_set)

% Include this line so program can run
Screen('Preference', 'SkipSyncTests', 1);

Priority(1)

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Initializing for PsychPortAudio
InitializePsychSound(1)

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenSubj = 2;     %subject screen
screenInv = 1;     %investigator screen

% Subject monitor either in normal orientation or "flipped" (rotated 180 d)
flipped_monitor = 1; 

% make sure monitor refresh is 144 hz
desired_refresh = 140;
actual_refresh = Screen('FrameRate',screenSubj);

if actual_refresh < desired_refresh
    disp('Fix monitor refresh rate!')
    return
end

% Define black and white
white = WhiteIndex(screenSubj);
black = BlackIndex(screenSubj);

% Open an on screen window
[windowSubj, windowRectSubj] = PsychImaging('OpenWindow', screenSubj, black);
[windowInv, windowRectInv] = PsychImaging('OpenWindow', screenInv, black);
WinTabMex(0, windowSubj); %Initialize tablet driver, connect it to 'win'
ListenChar(2)% 

% Query the frame duration
ifi = Screen('GetFlipInterval', windowSubj);

% Setup the text type for the window
Screen('TextFont', windowSubj, 'Ariel');
Screen('TextSize', windowSubj, 28);

% Get the centre coordinate of the window
[xCenterSubj, yCenterSubj] = RectCenter(windowRectSubj);
[xCenterInv, yCenterInv] = RectCenter(windowRectInv);

% Enable alpha blending for anti-aliasing
% For help see: Screen BlendFunction?
% Also see: Chapter 6 of the OpenGL programming guide
% Shows up as rect instead of smooth dot without this line. -HK
Screen('BlendFunction', windowSubj, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('BlendFunction', windowInv, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Load cursor image
[cursor_img, ~, cursor_alpha] = imread('cursor.png');
cursor_img(:,:,4) = cursor_alpha(:,:);
mm2pixel = 3.6137;
cursortext = Screen('MakeTexture', windowSubj, cursor_img);
cursor_r = 1.5*mm2pixel;

% There should be 2540 lines per inch (lpi) on tablet, and tablet active
% area is 19.2 in x 12.0 in. Need to convert tablet lines into pixel space. 
% If there are 2540 lines per inch on tablet, that means a line is equal to
% 1/2540th of an inch OR 1/1000th of a cm OR 1/100th of a mm. Therefore, 
% conversion factor for 1 line requires mm2pixel * 0.01. Offsets are 
% defined as origin of tablet (bottom left) relative to origin of monitor
% (top left, when oriented same way as tablet).
if flipped_monitor == 1
    % For flipped monitor:
    tablet_x_scale = -0.036137; % mm2pixel*0.01
    tablet_x_offset = 19.28*2540; % first value is offset in inches
    tablet_y_scale = 0.036137;
    tablet_y_offset = 0.2*2540;
else
    % For monitor in standard orientation
    tablet_x_scale = 0.036137;
    tablet_x_offset = -0.48*2540;
    tablet_y_scale = -0.036137;
    tablet_y_offset = 12.218*2540;
end

% Load target file
cd('target-files')

% Read in target file
tgt_file = dlmread([tgt_file_name_prefix,tgt_set,'.tgt'], '\t', 1, 0); % start reading in from 2nd row, 1st column
trial_num = tgt_file(:,1);
tgt_dist = tgt_file(:,2).*mm2pixel;
tgt_ang_1 = tgt_file(:,3);
tgt_ang_2 = tgt_file(:, 4);
rotation = tgt_file(:, 5);
gain = tgt_file(:, 6);
fixed = [];
ccw = 0; 
online_fb = tgt_file(:, 7);
endpoint_fb = tgt_file(:, 8);
clamped_feedback = tgt_file(:, 9);
between_blocks = tgt_file(:, 10);
targetsize = tgt_file(1, 11);
numtrials = size(tgt_file, 1);
maxtrialnum = max(numtrials);
% Target locations for target jumps
if flipped_monitor
    tgtx_1_subj = xCenterSubj - tgt_dist.*cosd(tgt_ang_1);
    tgty_1_subj = yCenterSubj + tgt_dist.*sind(tgt_ang_1);
    tgtx_2_subj = xCenterSubj - tgt_dist.*cosd(tgt_ang_2);
    tgty_2_subj = yCenterSubj + tgt_dist.*sind(tgt_ang_2);
else
    tgtx_1_subj = xCenterSubj + tgt_dist.*cosd(tgt_ang_1);
    tgty_1_subj = yCenterSubj - tgt_dist.*sind(tgt_ang_1);
    tgtx_2_subj = xCenterSubj + tgt_dist.*cosd(tgt_ang_2);
    tgty_2_subj = yCenterSubj - tgt_dist.*sind(tgt_ang_2);
end
tgtloc_1_subj = [tgtx_1_subj tgty_1_subj];
tgtloc_2_subj = [tgtx_2_subj tgty_2_subj];
tgtx_1_inv = xCenterInv + tgt_dist.*cosd(tgt_ang_1);
tgty_1_inv = yCenterInv - tgt_dist.*sind(tgt_ang_1);
tgtx_2_inv = xCenterInv + tgt_dist.*cosd(tgt_ang_2);
tgty_2_inv = yCenterInv - tgt_dist.*sind(tgt_ang_2);
tgtloc_1_inv = [tgtx_1_inv tgty_1_inv];
tgtloc_2_inv = [tgtx_2_inv tgty_2_inv];

% Set game constants
start_tolerance = 10*mm2pixel;  
movement_dist_thresh = 10*mm2pixel;
movement_time_thresh = 0.5;
startcirclewidth = 6*mm2pixel;
rt_dist_thresh = 5*mm2pixel;
rt_vel_thresh_start = 10*mm2pixel;
rt_vel_thresh_stop = 10*mm2pixel;
targetsize = targetsize*mm2pixel;
halfway_point = tgt_dist(1) / 2;
red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
magenta = [1 0 1];
yellow = [1 1 0];
insidetime = 0;
curtime = 0;
endptfbtime = 0.5;
pausetime = 0.0;
wait_time = 0.3 + (0.5-0.3).*rand(maxtrialnum,1); %uniformly distributed between 300-500ms
rt = 0;
RTs = [];
mt = 0;
MTs = [];
searchtime = 0;
SearchTimes = [];
fb_time = 0;
hits = 0;
jump = 0;
data = [];

gamephase = 0;
trial = 1;
cursor = [];

% hit will be any part of cursor touching target
hit_tolerance = targetsize./2 + cursor_r;

WinTabMex(2); %Empties the packet queue in preparation for collecting actual data

% Set the mouse to the center of the screen to start with
HideCursor;

% Define the ESC key
KbName('UnifyKeynames');
esc = KbName('ESCAPE');
space = KbName('SPACE');

% Variables that store data by trial
hand_angle = nan(maxtrialnum,1);
hand_dist_trial = nan(maxtrialnum,1);
hX_ep = nan(maxtrialnum,1);
hY_ep = nan(maxtrialnum,1);

% Variables that store data by sample
MAX_SAMPLES=6e6; %about 1 hour @ 1.6kHz = 60*60*1600
% timevec=nan(MAX_SAMPLES,1);
% delay_calc_time=nan(MAX_SAMPLES,1);
gamephase_move=nan(MAX_SAMPLES,1);
tablet_queue_length=nan(MAX_SAMPLES,1);
thePoints=nan(MAX_SAMPLES,2);
cursorPoints=nan(MAX_SAMPLES,2);
tabletPoints=uint16(nan(MAX_SAMPLES/8,2)); %reduce # of samples since the tablet is sampled @ 200Hz
tabletTime=nan(MAX_SAMPLES/8,1);
% total_vel=nan(MAX_SAMPLES,1);
% total_displacement=nan(MAX_SAMPLES,1);
% index_of_point_shown = nan(MAX_SAMPLES,1);
% [deltax,deltay]=deal(nan(MAX_SAMPLES,1));
dt_all = nan(MAX_SAMPLES,1);
t = nan(MAX_SAMPLES,1);
trial_time = nan(MAX_SAMPLES,1);
trial_move = nan(MAX_SAMPLES,1);
start_x_move = nan(MAX_SAMPLES,1);
start_y_move = nan(MAX_SAMPLES,1);
rotation_move = nan(MAX_SAMPLES,1);
hand_dist_all = nan(MAX_SAMPLES,1);

% Set sampling rate of game loop
desiredSampleRate = 500;
k = 0;
tab_k = 15;

tic;
begintime = GetSecs;
t0 = begintime;
nextsampletime = begintime;
too_long = 0;

% Loop game until over or ESC key press
while trial <= maxtrialnum
    
    % Exits experiment when ESC key is pressed.
    [keyIsDown, secs, keyCode] = KbCheck;
    if keyIsDown
        if keyCode(esc)
            break
        end
    end
    
    k = k+1;    %keep track of each loop iteration   
    t(k) = GetSecs - begintime;
    dt = toc-curtime;
    dt_all(k) = dt;
    
    if k == 1;
        trial_time(k) = dt;
    else
        trial_time(k) = trial_time(k-1) + dt;
    end
    curtime = toc;
    
    % Flip to the screen
    % last argument - 1: synchronous screen flipping, 2:asynchronous screen flipping
    Screen('Flip', windowSubj, 0, 0, 2);
    Screen('Flip', windowInv, 0, 0, 2);

    % Record trial number
    trial_move(k) = trial;
    rotation_move(k) = rotation(trial,1);
    
    % Read information from the tablet
    pkt = WinTabMex(5); % reads the latest data point out of a tablet's event queue
    tablet_queue_length(k) = 0;
    
    while ~isempty(pkt) % makes sure data are in packet; once pkt is 'grabbed,' then rest of code executes
        tabletPoints(tab_k,1:2) = pkt(1:2)';    % placing x and y (pkt rows 1,2) into tabletPoints variable
        tabletTime(tab_k) = (pkt(6)-tabletTime(16))/1000;   % tab_k initialized to 15; giving a little buffer at start of game?
        tab_k = tab_k+1;    % now tab_k is just another iterating variable
        tablet_queue_length(k) = tablet_queue_length(k)+1;  % adding each loop through
        pkt = WinTabMex(5); % reads the latest data point out of a tablet's event queue
    end
    
    % HAND COORDINATES
    % x,y coordinates from WinTabMex pkt
    hX = (double(tabletPoints(tab_k-1,1))-tablet_x_offset)*tablet_x_scale;
    hY = (double(tabletPoints(tab_k-1,2))-tablet_y_offset)*tablet_y_scale;
    hX_cur = (double(tabletPoints(tab_k-1,1))-tablet_x_offset)*tablet_x_scale*(-1);
    hY_cur = (double(tabletPoints(tab_k-1,2))-tablet_y_offset)*tablet_y_scale*(-1);
    
    thePoints(k,:) = [hX hY]; % record full precision points
    
    %store radial hand distance
    hand_dist = sqrt((hX-xCenterSubj)^2 + (hY-yCenterSubj)^2);
    hand_dist_all(k) = hand_dist;
    
    %use initialized value for first sample
    if k > 1
        delta_hand_dist = hand_dist_all(k) - hand_dist_all(k-1);
    else
        delta_hand_dist = 0;
    end
    
    %check for fresh sample; if there, calculate hand velocity and restart
    %'stopwatch'
    if  abs(delta_hand_dist) > 0
        t_between_samples = t(k) - t0;
        hand_vel = delta_hand_dist / t_between_samples;
        t0 = t(k);
    end
    
    % ROTATED CURSOR (including clamp)
    if clamped_feedback(trial,1) == 1
        % Clamped fb location:
        if flipped_monitor
            rcX = xCenterSubj - gain(trial,1)*hand_dist.*cosd(tgt_ang_1(trial,1) + rotation(trial,1)); % may need to subtract rotation in order to make (+) clamp CCW
            rcY = yCenterSubj + gain(trial,1)*hand_dist.*sind(tgt_ang_1(trial,1) + rotation(trial,1));
        else
            rcX = xCenterSubj + gain(trial,1)*hand_dist.*cosd(tgt_ang_1(trial,1) + rotation(trial,1)); % may need to subtract rotation in order to make (+) clamp CCW
            rcY = yCenterSubj - gain(trial,1)*hand_dist.*sind(tgt_ang_1(trial,1) + rotation(trial,1));
        end
    else
        [rcX_rotated, rcY_rotated] = rotatexy(round(hX)-xCenterSubj,(round(hY)-yCenterSubj),rotation(trial,1),gain(trial,1));
        rcX = rcX_rotated + xCenterSubj;
        rcY = rcY_rotated + yCenterSubj;

        [rcX_cur_rotated, rcY_cur_rotated] = rotatexy(round(hX_cur)-xCenterSubj,(round(hY_cur)-yCenterSubj),rotation(trial, 1),gain(trial, 1));
        rcX_cur = rcX_cur_rotated + xCenterSubj;
        rcY_cur = rcY_cur_rotated + yCenterSubj;
    end

    % Calculate distance of cursor from home
    cur_dist = sqrt((rcX-xCenterSubj)^2 + (rcY-yCenterSubj)^2);
    
    % Draw home position.
    if gamephase == 0   % Searching for start location
        if trial == 1
            [rcX_rotated, rcY_rotated] = rotatexy(round(hX)-xCenterSubj,(round(hY)-yCenterSubj),rotation(trial,1),gain(trial,1));
            rcX = rcX_rotated + xCenterSubj;
            rcY = rcY_rotated + yCenterSubj;
        else
            [rcX_rotated, rcY_rotated] = rotatexy(round(hX)-xCenterSubj,(round(hY)-yCenterSubj),rotation(trial - 1, 1),gain(trial - 1, 1));
            rcX = rcX_rotated + xCenterSubj;
            rcY = rcY_rotated + yCenterSubj;
        end

        searchtime = searchtime + dt;
        SearchTimes(trial) = searchtime;
        
        %draw start circle
        Screen('DrawDots', windowSubj, [xCenterSubj yCenterSubj], startcirclewidth, yellow, [], 2);
        Screen('DrawDots', windowInv, [xCenterInv yCenterInv], startcirclewidth, yellow, [], 2);

%         % Check whether trial had online fb or not
%         if (trial > 1) && (online_fb(trial - 1, 1) == 1)
%             visible = 1;
%         else
%             %if hand is within start tolerance, show cursor
%             if hand_dist < start_tolerance
%                 visible = 1;
%             else
%                 visible = 1; % Setting to always visible on trial 1
%             end
%     end

        % Check whether trial had online fb or not
        if (trial == 1) || (online_fb(trial - 1, 1) == 1)
            visible = 1;
        else
            if hand_dist < start_tolerance
                visible = 1;
            else
                visible = 0;
            end
        end
        
        %if hand is inside start position, start timer
        if cur_dist < startcirclewidth/2
            inside = 1;
            insidetime = insidetime + dt;
        else
            inside = 0;
            insidetime = 0;
        end
        
        %check to see if hand has held within start position appropriate
        %amount of time
        if inside == 1 && insidetime > wait_time(trial, 1)
            gamephase = 1;
            insidetime = 0;
            t0 = t(k);
        end
        
    elseif gamephase == 1  % Show target
        visible = online_fb(trial,1);
        rt = rt + dt;
        
        % Draw target
        Screen('DrawDots', windowSubj, tgtloc_1_subj(trial,:), targetsize, green, [], 2);      
        Screen('DrawDots', windowInv, tgtloc_1_inv(trial,:), targetsize, green, [], 2);
        % Draw start circle
        Screen('DrawDots', windowSubj, [xCenterSubj yCenterSubj], startcirclewidth, yellow, [], 2);
        Screen('DrawDots', windowInv, [xCenterInv yCenterInv], startcirclewidth, yellow, [], 2);

        % Check for movement initiation
        if (hand_vel >= rt_vel_thresh_start) && (hand_dist >= rt_dist_thresh)
            RTs(trial) = rt;
            gamephase = 2;
        end
        
    elseif gamephase == 2  % Moving towards target
        visible = online_fb(trial,1);
        mt = mt + dt;
        MTs(trial) = mt;

        % Jump the target when hand crosses half reach distance
        if (hand_dist <= halfway_point) && (jump == 0)
            Screen('DrawDots', windowSubj, tgtloc_1_subj(trial,:), targetsize, green, [], 2);
            Screen('DrawDots', windowInv, tgtloc_1_inv(trial,:), targetsize, green, [], 2);
        else
            Screen('DrawDots', windowSubj, tgtloc_2_subj(trial,:), targetsize, green, [], 2);
            Screen('DrawDots', windowInv, tgtloc_2_inv(trial,:), targetsize, green, [], 2);
            jump = 1;
        end

        %check for movement end
        if (hand_vel <= rt_vel_thresh_stop) && (hand_dist >= movement_dist_thresh)
            fb_angle = atan2d(rcY-yCenterSubj, rcX-xCenterSubj);
            fb_x = gain(trial,1)*hand_dist_all(k)*cosd(fb_angle) + xCenterSubj;
            fb_y = gain(trial,1)*hand_dist_all(k)*sind(fb_angle) + yCenterSubj;
            if flipped_monitor
                hand_angle(trial,1) = atan2d(hY-yCenterSubj, (hX-xCenterSubj)*-1);
            else
                hand_angle(trial,1) = atan2d((hY-yCenterSubj)*(-1), hX-xCenterSubj);
            end
            hand_dist_trial(trial,1) = hand_dist_all(k);
            hX_ep(trial,1) = hX;
            hY_ep(trial,1) = hY;
            gamephase = 3;
            if MTs(trial) >= movement_time_thresh
                too_long = 1;
                gamephase = 3; % If we want to automatically progress to
                               % next gamephase when MT is too long
            end
        end
        
    elseif gamephase == 3  % Endpoint feedback
        visible = endpoint_fb(trial,1);
        
        if too_long
            Screen('DrawDots', windowSubj, tgtloc_2_subj(trial,:), targetsize, red, [], 2);
            Screen('DrawDots', windowInv, tgtloc_2_inv(trial,:), targetsize, red, [], 2);
        else
            Screen('DrawDots', windowSubj, tgtloc_2_subj(trial,:), targetsize, green, [], 2);
            Screen('DrawDots', windowInv, tgtloc_2_inv(trial,:), targetsize, green, [], 2);
        end
        
        if fb_time <= endptfbtime
            fb_time = fb_time + dt;
        else
            gamephase = 4;
            visible = 0;
        end
        
    elseif gamephase == 4  % Between Blocks Message
        trial_time(k) = 0;
        if between_blocks(trial) == 1
            Screen('DrawText', windowSubj, 'Great job!' , xCenterSubj-100, yCenterSubj, white);
            Screen('DrawText', windowInv, 'Great job!' , xCenterInv-100, yCenterInv, white);
        elseif between_blocks(trial) == 2
            DrawFormattedText(windowSubj, 'Reach to the target as quickly and accurately as possible.', 'center', 'center', white, [], 1, 1);
            DrawFormattedText(windowInv, 'Reach to the target as quickly and accurately as possible.', 'center', 'center', white);
        elseif between_blocks(trial) == 3
            Screen('DrawText', windowSubj, 'You will no longer get to see the white cursor during your reach.' , xCenterSubj-300, yCenterSubj, white);
            Screen('DrawText', windowSubj, 'Continue to reach quickly and accurately to the target.' , xCenterSubj-300, yCenterSubj, white);
            Screen('DrawText', windowInv, 'You will no longer get to see the white cursor during your reach.' , xCenterInv-300, yCenterInv, white);
            Screen('DrawText', windowInv, 'Continue to reach quickly and accurately to the target.' , xCenterInv-300, yCenterInv, white);
        elseif between_blocks(trial) == 4
            Screen('DrawText', windowSubj, 'Every trial will start by bringing the white cursor to the center start position.' , xCenterSubj-430, yCenterSubj, white);
            Screen('DrawText', windowSubj, 'The white cursor will only appear when you are close to the start position.' , xCenterSubj-430, yCenterSubj, white);
            Screen('DrawText', windowInv, 'Every trial will start by bringing the white cursor to the center start position.' , xCenterInv-430, yCenterInv, white);
            Screen('DrawText', windowInv, 'The white cursor will only appear when you are close to the start position.' , xCenterInv-430, yCenterInv, white);
        elseif between_blocks(trial) == 5
            DrawFormattedText(windowSubj, 'Thank you for participating in our study!', 'center', 'center', white, [], 1, 1);
            DrawFormattedText(windowInv, 'Thank you for participating in our study!' ,  'center', 'center', white);
        else
            gamephase = 0;
            fb_time = 0;
            searchtime = 0;
            rt = 0;
            mt = 0;
            jump = 0;
            trial_time(k) = 0;
            trial = trial + 1;
            too_long = 0;
        end
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(space)
                trial = trial + 1;
                gamephase = 0;
                fb_time = 0;
                searchtime = 0;
                rt = 0;
                mt = 0;
                trial_time(k) = 0;
                too_long = 0;
            end
        end
    end
    
    % Draw actual hand position all phases on investigator screen
    if flipped_monitor
        cursorDistX = (hX - xCenterSubj)*-1;
        cursorDistY = (hY - yCenterSubj)*-1;
    else
        cursorDistX = hX - xCenterSubj;
        cursorDistY = hY - yCenterSubj;
    end
    cursorInv = [xCenterInv+cursorDistX yCenterInv+cursorDistY];
    Screen('DrawDots', windowInv, cursorInv, cursor_r*2, red, [], 2);
    
    % Draw Cursor
    if visible
        if (gamephase == 0)
            cursor = [(rcX - cursor_r) (rcY - cursor_r) (rcX + cursor_r) (rcY + cursor_r)];
            cursorPoints(k,:) = [rcX rcY]; % record full precision points
            cursor2 = [(xCenterInv+(rcX-xCenterSubj)*(-1) - cursor_r) (yCenterInv+(rcY-yCenterSubj)*(-1) - cursor_r) (xCenterInv+(rcX-xCenterSubj)*(-1) + cursor_r) (yCenterInv+(rcY-yCenterSubj)*(-1) + cursor_r)];            
            Screen('DrawTexture', windowSubj, cursortext, [], cursor, [], [], [],[255 255 255]);
            Screen('DrawTexture', windowInv, cursortext, [], cursor2, [], [], [],[255 255 255]);
        elseif (gamephase == 1 || gamephase == 2 || gamephase == 3)
            cursor = [(rcX - cursor_r) (rcY - cursor_r) (rcX + cursor_r) (rcY + cursor_r)];
            cursorPoints(k,:) = [rcX rcY]; % record full precision points
            cursor2 = [(xCenterInv+(rcX-xCenterSubj)*(-1) - cursor_r) (yCenterInv+(rcY-yCenterSubj)*(-1) - cursor_r) (xCenterInv+(rcX-xCenterSubj)*(-1) + cursor_r) (yCenterInv+(rcY-yCenterSubj)*(-1) + cursor_r)];            
            Screen('DrawTexture', windowSubj, cursortext, [], cursor, [], [], [],[255 255 255]);
            Screen('DrawTexture', windowInv, cursortext, [], cursor2, [], [], [],[255 255 255]);
%         elseif  (gamephase == 3)
%             cursor = [fb_x - cursor_r, fb_y - cursor_r, fb_x + cursor_r, fb_y + cursor_r];
%             cursorPoints(k,:) = [fb_x fb_y]; % record full precision points
%             cursor2 = [xCenterInv+(fb_x-xCenterSubj)*(-1) - cursor_r, yCenterInv+(fb_y-yCenterSubj)*(-1) - cursor_r, xCenterInv+(fb_x-xCenterSubj)*(-1) + cursor_r, yCenterInv+(fb_y-yCenterSubj)*(-1) + cursor_r];
%             Screen('DrawTexture', windowSubj, cursortext, [], cursor, [], [], [],[255 255 255]);
%             Screen('DrawTexture', windowInv, cursortext, [], cursor2, [], [], [],[255 255 255]);
        end
    else
        cursorPoints(k,:) = [NaN NaN];
    end

    gamephase_move(k) = gamephase;    
    sampletime(k) = GetSecs;
    nextsampletime = nextsampletime + 1/desiredSampleRate;
    
    while GetSecs < nextsampletime
    end
    
end

endtime = GetSecs
elapsedTime = endtime - begintime
numberOfSamples = k
actualSampleRate = 1/(elapsedTime / numberOfSamples)
thePoints(1:k,:);

ShowCursor;
% Clear the screen
sca;
WinTabMex(3); % Stop/Pause data acquisition.
WinTabMex(1); % Shutdown driver.
ListenChar(0);

% Game file
hand_angle = hand_angle;

% Movement file
trial_move = trial_move;
gamephase_move = gamephase_move;
t = t;
dt_all = dt_all;
rotation_move = rotation_move;
start_x_move = ones(k,1).*xCenterSubj;
start_y_move = ones(k,1).*yCenterSubj;
hand_x(:,1) = thePoints(:,1)-xCenterSubj;
hand_y(:,1) = (thePoints(:,2)-yCenterSubj)*(-1); % adjust for other monitor points
cursor_x(:,1) = cursorPoints(:,1) - xCenterSubj;
cursor_y(:,1) = (cursorPoints(:,2) - yCenterSubj)*(-1); % adjust for other monitor points

%
cd('..\..\data')

% Save data
name_prefix_all = [name_prefix, '_',tgt_file_name_prefix];
disp('Saving...')
if ~exist([name_prefix_all,'.mat'],'file')
    datafile_name = [name_prefix_all,'.mat'];
    
elseif ~exist([name_prefix_all,'_a.mat'],'file'), datafile_name = [name_prefix_all,'_a.mat'];
elseif ~exist([name_prefix_all,'_b.mat'],'file'), datafile_name = [name_prefix_all,'_b.mat'];
else
    char1='c';
    while exist([name_prefix_all,'_',char1,'.mat'],'file'), char1=char(char1+1); end
    datafile_name = [name_prefix_all,'_',char1,'.mat'];
end
save(datafile_name);
disp(['Saved ', datafile_name]);
%end


function [rx, ry] = rotatexy(x,y,phi,gain)
% phi is in degrees
phi=phi*pi/180;
[theta r]=cart2pol(x,y);
[rx ry]=pol2cart(theta-phi,gain*r);
return

