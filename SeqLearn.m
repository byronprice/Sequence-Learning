function [] = SeqLearn(AnimalName,Day,holdTime)
%SeqLearn.m
%  Display static full-field sinusoidal gratings for 150ms
%   12 unique elements (0 to 165 in increments of 15 degrees)
%   each element will be displayed for 150 ms
%   followed by a 0.5 to 1.5 second pause (uniform random distribution)

% could also try with the Berry Patch stimulus (choose say 10 unique
%  images and then you have 90 possible combinations)

% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%        Day - experimental day
%
%        Optional- 
%        holdTime - amount of time to wait between blocks of about 50 stimuli
% 
%
% OUTPUT: a file with stimulus parameters named SeqStimDate_AnimalName
%           e.g. SeqStim20160708_12345.mat to be saved in CloudStation's 
%           Seq folder under '~/CloudStation/ByronExp/SEQ'
% Created: 2018/03/03 at 24 Cummington, Boston
%  Byron Price
% Updated: 2018/03/03
%  By: Byron Price

cd('~/CloudStation/ByronExp/SEQ');

blocks = 4;
repsPerBlock = 50;
numElements = 2;
% orientations = [15,30,45,60,75,105,120,135,150,165].*pi./180;
% orientations = orientations(randperm(length(orientations),numElements));
orientations = [60,150]; % .*pi/180
numOrient = numElements;
stimTimes = 150/1000;
ISI = [0.5,1.5];
spatFreq = 0.05;
DistToScreen = 25;
gama = 2.1806;
degreeRadius = 179;
stimOnTime = 50/1000;

directory = '~/Documents/MATLAB/Byron/Sequence-Learning';

if nargin < 3
    holdTime = 30;
end

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date = str2double(Date);
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

usb = usb1208FSPlusClass;
display(usb);

WaitSecs(1);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% 
% % Open a fullscreen onscreen window on that display, choose a background
% % color of 127 = gray with 50% max intensity; 0 = black;255 = white
background = 127;
[win,~] = Screen('OpenWindow', screenid,background);

gammaTable = makeGrayscaleGammaTable(gama,0,255);
Screen('LoadNormalizedGammaTable',win,gammaTable);

% Switch color specification to use the 0.0 - 1.0 range
Screen('ColorRange', win, 1);

% Query window size in pixels
[w_pixels, h_pixels] = Screen('WindowSize', win);

% Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);

dgshader = [directory '/SequenceStim.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [directory '/SequenceStim.frag.txt'] }, 1);
gratingTex = Screen('SetOpenGLTexture', win, [], 0, GL.TEXTURE_3D,w_pixels,...
    h_pixels, 1, GratingShader);

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;
mmPerPixel = conv_factor;
conv_factor = 1/conv_factor;

% perform unit conversions
Radius = (tan((degreeRadius/2)*pi/180)*(DistToScreen*10*2))*conv_factor; % get number of pixels
     % that degreeRadius degrees of visual space will occupy
temp = (tan(((1/spatFreq)/2)*pi/180)*(DistToScreen*10*2))*conv_factor;
spatFreq = 1/temp;clear temp;

centerVals = [w_pixels/2,h_pixels/2];

if Day<5
    waitTimes = ISI(1)+(ISI(2)-ISI(1)).*rand([repsPerBlock*blocks,1]);
    stimParams = cell(blocks,5);
    
    for ii=1:blocks
        stimParams{ii,1} = numElements;
        stimParams{ii,2} = orientations;
        stimParams{ii,3} = stimTimes;
        stimParams{ii,4} = pi.*ones(numElements,1);%2*pi*rand([numEl,1]);
        stimParams{ii,5} = [1,2];
    end
    
    offsetGrey = numOrient+1;
    
    estimatedTime = ((mean(ISI)+mean(stimTimes))*repsPerBlock*blocks+blocks*holdTime+2)/60;
    fprintf('\nEstimated time: %3.2f minutes',estimatedTime);
    
    % Define first and second ring color as RGBA vector with normalized color
    % component range between 0.0 and 1.0, based on Contrast between 0 and 1
    % create all textures in the same window (win), each of the appropriate
    % size
    Grey = 0.5;
    Black = 0;
    White = 1;
    
    Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    
    % Perform initial flip to gray background and sync us to the retrace:
    Priority(9);
    
    usb.startRecording;WaitSecs(1);usb.strobeEventWord(0);
    WaitSecs(holdTime);
    
    % Animation loop
    count = 1;
    vbl = Screen('Flip',win);
    for yy=1:blocks
        numEl = stimParams{yy,1};
        currentOrient = stimParams{yy,2};
        currentPause = stimParams{yy,3};
        currentPhase = stimParams{yy,4};
        currentEvent = stimParams{yy,5};
        
        vbl = Screen('Flip',win);
        zz = 0;
        while zz < repsPerBlock
            ww = 0;
            while ww<numEl
                % ELEMENT on
                Screen('DrawTexture', win,gratingTex, [],[],...
                    [],[],[],[Grey Grey Grey Grey],...
                    [], [],[White,Black,...
                    Radius,centerVals(1),centerVals(2),spatFreq,currentOrient(ww+1),currentPhase(ww+1)]);
                % Request stimulus onset
                vbl = Screen('Flip',win,vbl-ifi/2+(currentPause-stimOnTime));
                usb.strobeEventWord(currentEvent(ww+1));
                vbl = Screen('Flip',win,vbl-ifi/2+stimOnTime);
                ww = ww+1;
            end
            usb.strobeEventWord(offsetGrey);
            vbl = Screen('Flip',win,vbl-ifi/2+stimOnTime);
            vbl = Screen('Flip',win,vbl-ifi/2+waitTimes(count));
            zz = zz+1;count = count+1;
        end
        if yy~=blocks
            timeIncrement = 1;
            totalTime = timeIncrement;
            while totalTime<holdTime
                usb.strobeEventWord(0);
                vbl = Screen('Flip',win,vbl-ifi/2+timeIncrement);
                totalTime = totalTime+timeIncrement;
            end
        end
    end
    WaitSecs(1);
    usb.stopRecording;
    Priority(0);
    
elseif Day==5
    conditions = 6;
    waitTimes = ISI(1)+(ISI(2)-ISI(1)).*rand([conditions*repsPerBlock*blocks,1]);
    stimParams = cell(conditions,blocks,5);
    
    A = orientations(1);B = orientations(2);C = orientations(1)-35;
    D = orientations(2)+20;
    % multiply by pi/180 later
    
    order = randperm(conditions,conditions);
    for jj=1:blocks
        stimParams{order(1),jj,1} = 2;
        stimParams{order(1),jj,2} = [A,B];
        stimParams{order(1),jj,3} = stimTimes;
        stimParams{order(1),jj,4} = pi.*ones(numElements,1);%2*pi*rand([numEl,1]);
        stimParams{order(1),jj,5} = [1,2];
    end
    
    for jj=1:blocks
        stimParams{order(2),jj,1} = 2;
        stimParams{order(2),jj,2} = [B,A];
        stimParams{order(2),jj,3} = stimTimes;
        stimParams{order(2),jj,4} = pi.*ones(numElements,1);%2*pi*rand([numEl,1]);
        stimParams{order(2),jj,5} = [2,1];
    end
    
    for jj=1:blocks
        stimParams{order(3),jj,1} = 2;
        stimParams{order(3),jj,2} = [C,D];
        stimParams{order(3),jj,3} = stimTimes;
        stimParams{order(3),jj,4} = pi.*ones(numElements,1);%2*pi*rand([numEl,1]);
        stimParams{order(3),jj,5} = [3,4];
    end
    
    for jj=1:blocks
        stimParams{order(4),jj,1} = 1;
        stimParams{order(4),jj,2} = A;
        stimParams{order(4),jj,3} = stimTimes;
        stimParams{order(4),jj,4} = pi.*ones(numElements,1);%2*pi*rand([numEl,1]);
        stimParams{order(4),jj,5} = 1;
    end
    
    for jj=1:blocks
        stimParams{order(5),jj,1} = 1;
        stimParams{order(5),jj,2} = B;
        stimParams{order(5),jj,3} = stimTimes;
        stimParams{order(5),jj,4} = pi.*ones(numElements,1);%2*pi*rand([numEl,1]);
        stimParams{order(5),jj,5} = 2;
    end
    
    for jj=1:blocks
        stimParams{order(6),jj,1} = 1;
        stimParams{order(6),jj,2} = C;
        stimParams{order(6),jj,3} = stimTimes;
        stimParams{order(6),jj,4} = pi.*ones(numElements,1);%2*pi*rand([numEl,1]);
        stimParams{order(6),jj,5} = 3;
    end
    
    offsetGrey = 5;
    
    estimatedTime = ((mean(ISI)+mean(stimTimes))*repsPerBlock*blocks*conditions+...
        conditions*blocks*holdTime+2)/60;
    fprintf('\nEstimated time: %3.2f minutes',estimatedTime);
    
    % Define first and second ring color as RGBA vector with normalized color
    % component range between 0.0 and 1.0, based on Contrast between 0 and 1
    % create all textures in the same window (win), each of the appropriate
    % size
    Grey = 0.5;
    Black = 0;
    White = 1;
    
    Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    
    % Perform initial flip to gray background and sync us to the retrace:
    Priority(9);
    
    usb.startRecording;WaitSecs(1);usb.strobeEventWord(0);
    WaitSecs(holdTime);
    
    % Animation loop
    count = 1;
    for xx=1:conditions
        vbl = Screen('Flip',win);
        for yy=1:blocks
            numEl = stimParams{xx,yy,1};
            currentOrient = stimParams{xx,yy,2};
            currentPause = stimParams{xx,yy,3};
            currentPhase = stimParams{xx,yy,4};
            currentEvent = stimParams{xx,yy,5};
            
            vbl = Screen('Flip',win);
            zz = 0;
            while zz < repsPerBlock
                ww = 0;
                while ww<numEl
                    % ELEMENT on
                    Screen('DrawTexture', win,gratingTex, [],[],...
                        [],[],[],[Grey Grey Grey Grey],...
                        [], [],[White,Black,...
                        Radius,centerVals(1),centerVals(2),spatFreq,currentOrient(ww+1),currentPhase(ww+1)]);
                    % Request stimulus onset
                    vbl = Screen('Flip',win,vbl-ifi/2+(currentPause-stimOnTime));
                    usb.strobeEventWord(currentEvent(ww+1));
                    vbl = Screen('Flip',win,vbl-ifi/2+stimOnTime);
                    ww = ww+1;
                end
                usb.strobeEventWord(offsetGrey);
                vbl = Screen('Flip',win,vbl-ifi/2+stimOnTime);
                vbl = Screen('Flip',win,vbl-ifi/2+waitTimes(count));
                zz = zz+1;count = count+1;
            end

            timeIncrement = 1;
            totalTime = timeIncrement;
            while totalTime<holdTime
                usb.strobeEventWord(0);
                vbl = Screen('Flip',win,vbl-ifi/2+timeIncrement);
                totalTime = totalTime+timeIncrement;
            end
        end
    end
    WaitSecs(1);
    usb.stopRecording;
    Priority(0);
    
end

cd('~/CloudStation/ByronExp/SEQ');
fileName = sprintf('SeqStim%d-%d_%d.mat',Day,Date,AnimalName);
save(fileName,'repsPerBlock','blocks','stimParams','stimTimes',...
    'w_pixels','h_pixels','spatFreq','mmPerPixel','waitTimes','holdTime',...
    'DistToScreen','orientations','offsetGrey','Day')
% Close window
Screen('CloseAll');

end

function gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
% Generates a 256x3 gamma lookup table suitable for use with the
% psychtoolbox Screen('LoadNormalizedGammaTable',win,gammaTable) command
% 
% gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
%
%   gamma defines the level of gamma correction (1.8 or 2.2 common)
%   blackSetPoint should be the highest value that results in a non-unique
%   luminance value on the monitor being used (sometimes values 0,1,2, all
%   produce the same black pixel value; set to zero if this is not a
%   concern)
%   whiteSetPoint should be the lowest value that returns a non-unique
%   luminance value (deal with any saturation at the high end)
% 
%   Both black and white set points should be defined on a 0:255 scale

gamma = max([gamma 1e-4]); % handle zero gamma case
gammaVals = linspace(blackSetPoint/255,whiteSetPoint/255,256).^(1./gamma);
gammaTable = repmat(gammaVals(:),1,3);
end