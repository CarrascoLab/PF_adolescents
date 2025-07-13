function myscreen = DPF_ExoAttnPractice(observer,IndContrast,Eye)

global stimulus;
global MGL;

thisdir = pwd;
% make a data directory if necessary
if ~isdir(fullfile(thisdir,'data'))
    disp('Making data directory');
    mkdir('data');
end

% make an observer directory if necessary
datadirname = fullfile(thisdir,'data',observer);
if ~isdir(datadirname);
    disp(sprintf('Making observer directory %s',datadirname));
    mkdir(datadirname);
end

disp(sprintf('[ DATA ] saving data in: %s',datadirname));


stimulus = [];
stimulus.contrasts=IndContrast;


% clearing old variables:
clear task myscreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initalize the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.EyeTrack=Eye;
myscreen.datadir = datadirname;
myscreen.allowpause = 0;
myscreen.saveData = -2;
myscreen.background=.5;
myscreen.keyboard.nums = [84 85];

% initalize the screen
myscreen = initScreen(myscreen);
if stimulus.EyeTrack
    myscreen = eyeCalibDisp(myscreen);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task{1}.waitForBacktick = 1;
task{1}.segmin =     [1 Inf .06 .06 .12 .04 .66 5 5];
task{1}.segmax =     [1 Inf .06 .06 .12 .04 .66 5 5];
task{1}.getResponse = [0 0 0 0 0 0 0 1 0];

n_repeats = 2; % 16 trials per block

[contrast, location, ExoCueCondition, repeat] = ndgrid(1, 1:4,1:2,1:n_repeats);
task{1}.numTrials = length(location(:));
random_order = randperm(task{1}.numTrials);
task{1}.randVars.contrast = contrast(random_order);
task{1}.randVars.ExoCueCondition = ExoCueCondition(random_order);
task{1}.randVars.targetLocation = location(random_order);
task{1}.randVars.uniform.targetOrientation = 1:2;
task{1}.randVars.uniform.distractorOrientation1 = 1:2;
task{1}.randVars.uniform.distractorOrientation2 = 1:2;
task{1}.randVars.uniform.distractorOrientation3 = 1:2;
task{1}.randVars.len_ = task{1}.numTrials;
stimulus.trialend = 0;
stimulus.trialnum=1;
stimulus.FixationBreak=zeros(1,length(location(:)));
stimulus.LocationIndices=unique(location);

task{1}.random = 1;
[task{1}, myscreen] = initTask(task{1},myscreen,@StartSegmentCallback,@DrawStimulusCallback,@responseCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myscreen = initStimulus('stimulus',myscreen);
stimulus = myInitStimulus(stimulus,myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    % runs automatically the task, you only need to change: StartSegmentCallback,DrawStimulusCallback,responseCallback
    [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end
clear stimulus.tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = endTask(myscreen,task);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = StartSegmentCallback(task, myscreen);
% segments: 1:ITI,   2:fixation,    3:stimulus, 4:response
global stimulus;

if (task.thistrial.thisseg == 1) % ITI
    stimulus.trialend = stimulus.trialend + 1;
elseif (task.thistrial.thisseg == 2) % fixation
    stimulus.tmp.targetLocation  = stimulus.eccentricity*[stimulus.locations{task.thistrial.targetLocation}];
    stimulus.tmp.distractorIndices=stimulus.LocationIndices(stimulus.LocationIndices~=task.thistrial.targetLocation);
    for Locs=1:length(stimulus.tmp.distractorIndices)
        stimulus.tmp.distractorLocations{Locs}= stimulus.eccentricity*[stimulus.locations{stimulus.tmp.distractorIndices(Locs)}];
    end
    stimulus.FixationStarted=0;
    %response cue
    stimulus.tmp.respcueLocation=stimulus.EndocueLocation{task.thistrial.targetLocation};
    stimulus.CueSounded=0;
elseif (task.thistrial.thisseg == 9) % response
    stimulus.trialnum = stimulus.trialnum + 1;
end
mglClearScreen(stimulus.grayColor);
setGammaTable(1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = DrawStimulusCallback(task, myscreen);
global stimulus;

mglClearScreen(stimulus.grayColor);%###

if (task.thistrial.thisseg == 1) % ITI
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
        end
    end
elseif (task.thistrial.thisseg == 2) % FIXATION
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
        end
    end
    if stimulus.EyeTrack
        ep=myscreen.eyetracker.eyepos;
        if (sqrt(ep(end,1)^2+ep(end,2)^2))<=stimulus.TrialStartFixDist && ~stimulus.FixationStarted
            stimulus.FixationStart=mglGetSecs;
            stimulus.FixationStarted=1;
        elseif (sqrt(ep(end,1)^2+ep(end,2)^2))>=stimulus.TrialStartFixDist && stimulus.FixationStarted
            stimulus.FixationStarted=0;
        elseif (sqrt(ep(end,1)^2+ep(end,2)^2))<=stimulus.TrialStartFixDist
            stimulus.FixationDur=mglGetSecs(stimulus.FixationStart);
            if stimulus.FixationDur >=stimulus.TrialStartFixDur
                task = jumpSegment(task);
            end
        end
    else
        mglWaitSecs(.25)
        task = jumpSegment(task);
    end
    
elseif (task.thistrial.thisseg == 3) % Exo cue
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
        end
    end
    
    if task.thistrial.ExoCueCondition==2
        for Gabor=stimulus.LocationIndices
            for corner=2:5
                mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.ExoCueSize,stimulus.white)
            end
        end
    elseif task.thistrial.ExoCueCondition==1
            for Gabor=task.thistrial.targetLocation
                for corner=2:5
                    mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.ExoCueSize,stimulus.white)
                end
            end
    end
    
elseif (task.thistrial.thisseg == 4) % ISI
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
        end
    end
    
elseif (task.thistrial.thisseg == 5) % STIMULUS
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
        end
    end
    if stimulus.EyeTrack
        ep=myscreen.eyetracker.eyepos;
        if (sqrt(ep(end,1)^2+ep(end,2)^2))>stimulus.TrialStartFixDist
            stimulus.FixationBreak(stimulus.trialnum)=1;
        end
    end
    drawGabor(stimulus.contrasts(task.thistrial.contrast),stimulus.tmp.targetLocation, stimulus.rotation(task.thistrial.targetOrientation), 1);
    for Loc=1:length(stimulus.tmp.distractorIndices)
        eval(sprintf('drawGabor(stimulus.contrasts(task.thistrial.contrast),stimulus.tmp.distractorLocations{Loc}, stimulus.rotation(task.thistrial.distractorOrientation%g), 1);',Loc));
    end
    
elseif (task.thistrial.thisseg == 6) % ISI-2
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
        end
    end
    
    
elseif (task.thistrial.thisseg == 7) % RESPONSE CUE
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
        end
    end
    mglLines2(stimulus.tmp.respcueLocation(1), stimulus.tmp.respcueLocation(2), stimulus.tmp.respcueLocation(3), stimulus.tmp.respcueLocation(4),2.5,stimulus.white);
    
elseif (task.thistrial.thisseg == 8) % RESPONSE WINDOW
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
        end
    end
    mglLines2(stimulus.tmp.respcueLocation(1), stimulus.tmp.respcueLocation(2), stimulus.tmp.respcueLocation(3), stimulus.tmp.respcueLocation(4),2.5,stimulus.white);
    
    if ~stimulus.CueSounded
        mglPlaySound(stimulus.CueSound);
        stimulus.CueSounded=1;
    end
    
elseif (task.thistrial.thisseg == 9)
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
        end
    end
    if stimulus.trialnum<=task.numTrials% End of block Feedback
        task = jumpSegment(task);
    else
        CorrectCount=0;
        IncorrectCount=0;
        for Stars=1:length(stimulus.starColorFeedback)
            if stimulus.starColorFeedback(Stars)==2
                CorrectCount=CorrectCount+1;
            elseif stimulus.starColorFeedback(Stars)==3
                IncorrectCount=IncorrectCount+1;
            end
        end
        PercentCorrect=([num2str(round(100*(CorrectCount/(task.numTrials)))) num2str('%')]);
        NumFixBreak=num2str(sum(stimulus.FixationBreak));
        mglTextSet('Courier',50,[0 0 0]);
        eval(sprintf('mglTextDraw(''%s correct.'',[0 10]);',PercentCorrect));
        eval(sprintf('mglTextDraw(''%s fixation breaks.'',[0 -10]);',NumFixBreak));
    end
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get the observer's response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = responseCallback(task, myscreen)
global stimulus;
mglClearScreen(stimulus.grayColor); %###
if ~task.thistrial.gotResponse
    
    % check response correct or not
    stimulus.tmp.response = [task.thistrial.whichButton == (task.thistrial.targetOrientation)]; %1 for left and 2 for right
    
    % give feeback:
    if stimulus.tmp.response
        mglPlaySound(stimulus.CorrectSound);
        stimulus.starColorFeedback(stimulus.trialnum)=2;
        task = jumpSegment(task);
    else
        mglPlaySound(stimulus.IncorrectSound);
        stimulus.starColorFeedback(stimulus.trialnum)=3;
        task = jumpSegment(task);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to draw the gabor stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawGabor(desiredContrast,position,orientation,sf);
% drawGaborPedCdeltaC
%
%        $Id: drawGabor.m, v 1 2007/01/18 19:40:56 ?? ?? $
%      usage: drawGabor(desiredContrast,position,orientation,sf)
%    purpose: draw a gabor stimulus on the screen with a specified contrast
%             within a given clut range (it finds the closest contrast to
%             the requested one within the available clut range)

global stimulus;
% now find closest matching contrast we can display with this gamma table
displayContrastNum = min(round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast),stimulus.nDisplayContrasts);
% disp(sprintf('Desired contrast: %0.4f Actual contrast: %0.4f Difference: %0.4f',desiredContrast,stimulus.currentMaxContrast*(displayContrastNum/stimulus.nDisplayContrasts),desiredContrast-stimulus.currentMaxContrast*(displayContrastNum/stimulus.nDisplayContrasts)));
if round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>stimulus.nDisplayContrasts
    disp(sprintf('[drawGabor] Desired contrast out of range %0.2f > %0.2f',desiredContrast,stimulus.currentMaxContrast));
    keyboard
end

mglBltTexture(stimulus.tex{sf}(displayContrastNum+1),position,0,0,orientation); %mglBltTexture(texture,position,hAlignment,vAlignment,rotation)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to create a gamma table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTable(maxContrast)
global stimulus;

% set the reserved colors
gammaTable(1:size(stimulus.reservedColors,1),1:size(stimulus.reservedColors,2))=stimulus.reservedColors;
% create the gamma table
cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
luminanceVals = cmin:((cmax-cmin)/(stimulus.nGratingColors-1)):cmax;
% now get the linearized range
redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');

% check the table
% plot(stimulus.linearizedGammaTable.redTable,'k');
% hold on
% plot(256*(0.25:0.5/250:0.75),redLinearized,'r');

% set the gamma table
gammaTable((stimulus.minGratingColors+1):256,:)=[redLinearized;greenLinearized;blueLinearized]';
% set the gamma table
mglSetGammaTable(gammaTable);
% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task)
global MGL;

% let's get the linearized gamma table
stimulus.linearizedGammaTable = mglGetGammaTable;
stimulus.linearizedGammaTable.redTable(1:3) = 0; % this is just to provisionally deal with what appears to be some bug: the first value in each of these gamma tables is a NaN
stimulus.linearizedGammaTable.greenTable(1:3) = 0;
stimulus.linearizedGammaTable.blueTable(1:3) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimulus parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gabors
stimulus.width = 4*.8;%stimulus.gaussSdx*7;             % in deg
stimulus.height = 4*.8;%stimulus.gaussSdy*7;            % in deg
% stimulus.gaussSdx = 1; %0.8;  %0.3; %0.5; %1; % stimulus.width/7;                % in deg
% stimulus.gaussSdy = 1; %0.8;  %0.3; %0.5; %1; % stimulus.height/7;               % in deg
stimulus.gaussSdx = stimulus.width/7;                % in deg
stimulus.gaussSdy = stimulus.height/7;               % in deg

stimulus.Tilt=20; % degrees from vertical
stimulus.rotation = [stimulus.Tilt -stimulus.Tilt]; % this is the tilt orientation of the gabor stimulus from vertical in Degrees
stimulus.DistractorRotation = [0];
stimulus.init = 1;

stimulus.sf = 4;    %2; %5;            % in cpd
stimulus.orientation = 90;      % in deg
stimulus.phase = 0;             % in deg
stimulus.eccentricity = 8*.8;  % in deg

% fixation
stimulus.FCwidth = 0.5;
stimulus.FClinewidth = 1.5;
stimulus.TrialStartFixDist=2; %2 degree radius in which to fixate befire trial starts
stimulus.TrialStartFixDur=.25;
stimulus.cornerDist=(stimulus.width/2);
stimulus.placeHolderSize=.05;%presents stimuli after this duration when fixation detected
stimulus.ExoCueDist=(stimulus.width/2)+.5;
stimulus.ExoCueSize=.16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frames and locations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.locations = {[cosd(90), sind(90)],[cosd(0), sind(0)],[cosd(-90), sind(-90)],[cosd(-180), sind(-180)]}; %4 locations: x,y coordinates are specified here. starts in N and moves clockwise
for i=1:length(stimulus.locations)
    stimulus.placeholders{i}= [stimulus.eccentricity*stimulus.locations{i}];
    stimulus.placeholders{i}(2,:)=[stimulus.placeholders{i}(1,1:2)]+[0 stimulus.cornerDist];
    stimulus.placeholders{i}(3,:)=[stimulus.placeholders{i}(1,1:2)]+[stimulus.cornerDist 0];
    stimulus.placeholders{i}(4,:)=[stimulus.placeholders{i}(1,1:2)]+[0 -stimulus.cornerDist];
    stimulus.placeholders{i}(5,:)=[stimulus.placeholders{i}(1,1:2)]+[-stimulus.cornerDist 0];
end
% Note: locations should specify a vector of length 1, i.e.,
%                       sum(locations.^2)==1
stimulus.frameThick = .08;
stimulus.reservedColors = [0 0 0; 1 1 1; 0 .6 0];

stimulus.precue.size  = 1.0;% for the peripheral cues
stimulus.precue.color = [1 1 1]; % white



% Response cue
stimulus.rcue.XLocation{1} = [-.6; 0; -0.4; 0];
stimulus.rcue.YLocation{1} = [ 0; .6; 0; -.6];
stimulus.rcue.XLocation{2} = [.7; 0; 0.4; 0];
stimulus.rcue.YLocation{2} = [ 0; .6; 0; -.6];
stimulus.rcue.color        = [0; .6; 0]; % green
row=[4 5 2 3];
signs{1}=[-1 -1];signs{2}=[-1 1];signs{3}=[1 1];signs{4}=[1 -1];
for i=1:length(stimulus.LocationIndices)
    stimulus.EndocueLocation{i}(1:2)=.8*[stimulus.locations{i}];
    stimulus.EndocueLocation{i}(3:4)=1.6*[stimulus.locations{i}];
    stimulus.NeutralcueLocation{i}(1:2)=.8*[stimulus.locations{i}];
    stimulus.NeutralcueLocation{i}(3:4)=1.2*[stimulus.locations{i}];
    stimulus.ExocueLocation{i}=(stimulus.eccentricity+stimulus.ExoCueDist)*stimulus.locations{i};
end
stimulus.ExocueLocation{length(stimulus.LocationIndices)+1}=[0 0];
stimulus.respcueLocation{1}=[0.8785 0;1.1714 0];
stimulus.respcueLocation{2}=[-0.8785 0;-1.1714 0];
stimulus.respcueLocation{3}=[0.8785 0;1.1714 0];
stimulus.respcueLocation{4}=[-0.8785 0;-1.1714 0];

% tentative_threshold = 0.1;
% log_steps = 1; % i.e., step=1/2 means every 2 steps is equal to doubling the contrast. step=1 means every step doubles the contrast, step=1/3 means every 3 steps double the contrast!
% % step and number of contrast levels determine the range. We want a range
% % that covers very low (i.e., performance is at chance) and very high
% % (i.e., performance at ceiling).
%
% % % calculate a logarithmic scale, centered around the average threshold
% stimulus.contrasts = 2.^((-3:3)*log_steps)*tentative_threshold;
% MG changed: stimulus.contrasts =[0.10 0.20 0.30 0.40 0.60 0.80 0.90];
%stimulus.contrasts =.3; %changed to use 1 contrast level

% some computations:
% for contrast
% idea: for low contrast stimuli, take advantage of the color-resolution of the
% colormap (10 bits) instead of the limited resolution of the grayspace (8 bit)
% that is available. This minimized the quatization error (causes
% artifacts, in particular in the periphery of the Gabor) by a factor of 4
% (2 bits). The same idea is achievable by dithering (adding random noise)
stimulus.nReservedColors=size(stimulus.reservedColors,1);
stimulus.nGratingColors = 256-(2*floor(stimulus.nReservedColors/2)+1);
stimulus.minGratingColors = 2*floor(stimulus.nReservedColors/2)+1;
stimulus.midGratingColors = stimulus.minGratingColors+floor(stimulus.nGratingColors/2);
stimulus.maxGratingColors = 255;
stimulus.deltaGratingColors = floor(stimulus.nGratingColors/2);

% to set up color values
stimulus.black = [0 0 0];
stimulus.white = [1/255 1/255 1/255];
stimulus.green = [0 6/255 0];
stimulus.grey = [.025 .025 .025];
stimulus.background = [stimulus.midGratingColors/255 stimulus.midGratingColors/255 stimulus.midGratingColors/255];
stimulus.fixation.color = [0; .6; 0]'; % green

% calculate a grating, a gaussian envelope (gaussian is in the alpha
% channel), and a mask (for now, just high-contrast random noise)
for thisSF = 1:length(stimulus.sf)      %only one spatial frequency
    gratingMatrix{thisSF} = mglMakeGrating(stimulus.width,stimulus.height,stimulus.sf(thisSF),stimulus.orientation,stimulus.phase);
end

grating(:,:,4) = 255*mglMakeGaussian(stimulus.width,stimulus.height,stimulus.gaussSdx,stimulus.gaussSdy);

% making the texture for all the Gabor stimuli:
disppercent(-inf,'Calculating gabors');
for thisSF = 1:length(stimulus.sf)
    for thisContrast = 0:stimulus.deltaGratingColors
        %% stimulus.texture
        grating(:,:,1) = stimulus.midGratingColors+gratingMatrix{thisSF}*thisContrast;
        grating(:,:,2) = grating(:,:,1);
        grating(:,:,3) = grating(:,:,1);
        stimulus.tex{thisSF}(thisContrast+1) = mglCreateTexture(grating);
        disppercent(thisContrast/stimulus.deltaGratingColors);
    end
end
disppercent(inf);
stimulus.nDisplayContrasts = stimulus.deltaGratingColors;
disppercent(inf);


% calculate gray color
stimulus.grayColor = stimulus.background; %stimulus.midGratingColors/255;

% sounds
stimulus.CueSound = find(strcmp(MGL.soundNames,'Ping'));
stimulus.CorrectSound = find(strcmp(MGL.soundNames,'Submarine'));
stimulus.IncorrectSound = find(strcmp(MGL.soundNames,'Basso'));
