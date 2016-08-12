function [] = SequenceStim(AnimalName,holdTime)
%SequenceStim.m
%  Display a sequence of sinusoidal white/black circles on a gray background,
%   set around the center of mass of the already-recorded retinotopy of an LFP
%   recording electrode. This code will coordinate with the Retinotopy.m
%   code and the saved file RetinoMapAnimalName.mat, e.g.
%   RetinoMap26881.mat .
%  Each circle will occupy a 5-degree radius of visual space
% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%
%        Optional- 
%        holdTime - amount of time to wait between blocks of about 50 stimuli
% 
%        see file SequenceVars.mat for other changeable presets
%
% OUTPUT: a file with stimulus parameters named SeqStimDate_AnimalName
%           e.g. SeqStim20160708_12345.mat to be saved in CloudStation's 
%           SeqExp folder under '~/CloudStation/ByronExp/SeqExp'
% Created: 2016/07/25 at 24 Cummington, Boston
%  Byron Price
% Updated: 2016/08/04
%  By: Byron Price

cd('~/CloudStation/ByronExp/RetinoExp');
load(sprintf('RetinoMap%d.mat',AnimalName));

cd('~/CloudStation/ByronExp/SeqExp');
load('SequenceVars.mat');

directory = '/home/jglab/Documents/MATLAB/Byron/Sequence-Learning';

if nargin < 2
    holdTime = 30;
end
reps = reps-mod(reps,blocks);

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date = str2double(Date);
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

usb = usb1208FSPlusClass;
display(usb);

WaitSecs(10);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% Open a fullscreen onscreen window on that display, choose a background
% color of 127 = gray with 50% max intensity; 0 = black
colorRange = 0:1:255;
[~,index] = min(abs(colorRange.^gama-255^gama/2));
background = colorRange(index)-6;
[win,~] = Screen('OpenWindow', screenid,background);

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
Radius = (tan(degreeRadius*pi/180)*(DistToScreen*10))*conv_factor; % get number of pixels
     % that degreeRadius degrees of visual space will occupy
Radius = round(Radius);
temp = (tan((1/spatFreq)*pi/180)*(DistToScreen*10))*conv_factor;
spatFreq = 1/temp;clear temp;

% extract the angle of the 2-D Gaussian to target the stimulus to the
% strongest responding locations within the retinotopic map
[V,D] = eig(squeeze(Sigma(Channel,:,:)));
[~,index] = max(max(D));
offset = atan2(V(2,index),V(1,index));

if offset < 0 
    offset = offset+2*pi;
end

centerVals = zeros(numElements,2);
degreeDiv = (2*pi)/numElements;
centerVals(1,1) = centerMass(Channel,1);centerVals(1,2) = centerMass(Channel,2);
for ii=2:numElements
    centerVals(ii,1) = centerMass(Channel,1)+cos(degreeDiv*(ii-2)+offset)*2*Radius;
    centerVals(ii,2) = centerMass(Channel,2)+sin(degreeDiv*(ii-2)+offset)*2*Radius;
end
temp2 = centerVals(2,:);temp3 = centerVals(3,:);
centerVals(2,:) = centerVals(1,:);
centerVals(3,:) = temp2;
centerVals(1,:) = temp3;

estimatedTime = ((stimTime*numElements+waitTime)*reps+blocks*holdTime)/60;
display(sprintf('\nEstimated time: %3.2f minutes',estimatedTime));

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Grey = 0.5^(1/gama);
Black = 0;
White = 1;

Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

orient = [150,60,30,120].*pi./180;
% Perform initial flip to gray background and sync us to the retrace:
Priority(9);

usb.startRecording;WaitSecs(1);usb.strobeEventWord(0);
WaitSecs(holdTime);

% Animation loop
for yy=1:blocks
    vbl = Screen('Flip',win);
    for zz = 1:reps/blocks
        for ii=1:numElements
            % Draw the procedural texture as any other texture via 'DrawTexture'
            Screen('DrawTexture', win,gratingTex, [],[],...
                [],[],[],[Grey Grey Grey Grey],...
                [], [],[White,Black,...
                Radius,centerVals(ii,1),centerVals(ii,2),spatFreq,orient(ii),gama]);
            % Request stimulus onset
            vbl = Screen('Flip', win,vbl+ifi/2);usb.strobeEventWord(ii);
            vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
        end
        usb.strobeEventWord(5);
        vbl = Screen('Flip',win,vbl-ifi/2+waitTime);
    end
    if yy ~= blocks
        usb.strobeEventWord(0);
        vbl = Screen('Flip',win,vbl-ifi/2+holdTime);
    end
end
WaitSecs(2);
usb.stopRecording;
Priority(0);

cd('~/CloudStation/ByronExp/SeqExp');
fileName = sprintf('SeqStim%d_%d.mat',Date,AnimalName);
save(fileName,'centerVals','Radius','reps','stimTime','numElements',...
    'w_pixels','h_pixels','spatFreq','mmPerPixel','waitTime','holdTime',...
    'DistToScreen','orient')
% Close window
Screen('CloseAll');

end
