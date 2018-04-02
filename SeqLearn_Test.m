function [] = SeqLearn_Test(AnimalName,holdTime)
%SeqLearn_Test.m
%  Display static full-field sinusoidal gratings for 100-300ms
%   12 unique elements (0 to 165 in increments of 15 degrees)
%   each element will be displayed for 100 to 300 msec (randomly chosen),
%   followed by a 0.5 to 1.5 second pause (uniform random distribution)

%   2-element sequences will also be presented (132 possible combinations)

% could also try with the Berry Patch stimulus (choose say 10 unique
%  images and then you have 90 possible combinations)

% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%
%        Optional- 
%        holdTime - amount of time to wait between blocks of about 50 stimuli
%
% OUTPUT: a file with stimulus parameters named SeqStimDate_AnimalName
%           e.g. SeqTestStim20160708_12345.mat to be saved in CloudStation's 
%           Seq folder under '~/CloudStation/ByronExp/Seq'
% Created: 2018/02/20 at 24 Cummington, Boston
%  Byron Price
% Updated: 2018/04/02
%  By: Byron Price

blocks = 100;
repsPerBlock = 10;
numElements = [1,2];
orientations = (0:15:165).*pi./180;numOrient = length(orientations);
stimTimes = (100:16.66666667:300)./1000;
ISI = [0.5,1.5];
spatFreq = 0.05;
DistToScreen = 25;
gama = 2.1806;
degreeRadius = 179;
radianRadius = degreeRadius*pi/180;
stimOnTime = 50/1000;

directory = '~/Documents/MATLAB/Byron/Sequence-Learning';

if nargin < 2
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

dgshader = [directory '/SequenceStim2.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [directory '/SequenceStim2.frag.txt'] }, 1);
gratingTex = Screen('SetOpenGLTexture', win, [], 0, GL.TEXTURE_3D,w_pixels,...
    h_pixels, 1, GratingShader);

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;
mmPerPixel = conv_factor;

spatFreq = spatFreq*180/pi;
DistToScreenPix = DistToScreen*10/mmPerPixel;

centerVals = [w_pixels/2,90/mmPerPixel];
centerPos = [0,0].*pi/180;

waitTimes = ISI(1)+(ISI(2)-ISI(1)).*rand([repsPerBlock*blocks,1]);
stimParams = cell(repsPerBlock*blocks,5);

for ii=1:blocks
   ind1 = binornd(1,0.9)+1;
   numEl = numElements(ind1);
   stimParams{ii,1} = numEl; % 1 or 2 element sequence
   ind2 = random('Discrete Uniform',numOrient,[numEl,1]);
   stimParams{ii,2} = orientations(ind2); % stim orientation
   ind3 = randperm(length(stimTimes),1);
   stimParams{ii,3} = stimTimes(ind3); % stimulus timing
   stimParams{ii,4} = pi/3.*ones(numEl,1);%2*pi*rand([numEl,1]); % phase
   stimParams{ii,5} = ind2;
end

offsetGrey = numOrient+1;

estimatedTime = ((mean(ISI)+mean(stimTimes)+2*stimOnTime)*repsPerBlock*blocks+5*holdTime+5)/60;
fprintf('\nEstimated time: %3.2f minutes\n',estimatedTime);

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
                    radianRadius,centerVals(1),centerVals(2),spatFreq,currentOrient(ww+1),...
                    currentPhase(ww+1),DistToScreenPix,centerPos(1),centerPos(2),0]);
            % Request stimulus onset
            vbl = Screen('Flip',win,vbl-ifi/2+(currentPause-stimOnTime));
            usb.strobeEventWord(currentEvent(ww+1));
            vbl = Screen('Flip',win,vbl-ifi/2+stimOnTime);
            ww = ww+1;
        end
        usb.strobeEventWord(offsetGrey);
        vbl = Screen('Flip',win,vbl-ifi/2+waitTimes(count));
        zz = zz+1;count = count+1;
    end
    if mod(yy,blocks/4)==0
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

cd('~/CloudStation/ByronExp/SEQ');
fileName = sprintf('SeqTestStim%d_%d.mat',Date,AnimalName);
save(fileName,'repsPerBlock','blocks','stimParams','stimTimes',...
    'w_pixels','h_pixels','spatFreq','mmPerPixel','waitTimes','holdTime',...
    'DistToScreen','orientations','offsetGrey')
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
