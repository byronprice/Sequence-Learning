function [] = SequenceStim(AnimalName,numElements,DistToScreen,degreeRadius)
%SequenceStim.m
%  Display a sequence of white circles on a black background, set around
%   the center of mass of the already-recorded retinotopy of an LFP
%   recording electrode. This code will coordinate with the Retinotopy.m
%   code and the saved file RetinoMapAnimalName.mat, e.g.
%   RetinoMap26881.mat .
%  Each circle will occupy a 5-degree radius of visual space
% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%
%        Optional- 
%        numElements - number of elements in the sequence
%        DistToScreen - physical distance of observer from the screen, in
%           units of cm
%        degreeRadius - degrees of visual field that radius of square will occupy
%        spatFreq - spatial frequency of oriented sinusoidal grating
%
% OUTPUT: a file with stimulus parameters named SeqStimDate_AnimalName
%           e.g. SeqStim20160708_12345.mat to be saved in the RetinoExp
%           folder under '/MATLAB/Byron/SeqExp'
% Created: 2016/07/25 at 24 Cummington, Boston
%  Byron Price
% Updated: 2016/07/26
%  By: Byron Price

cd('/home/jglab/Documents/MATLAB/Byron/RetinoExp/');
load(strcat('RetinoMap',num2str(AnimalName),'.mat'));

directory = '/home/jglab/Documents/MATLAB/Byron/Sequence-Learning/';
if nargin < 2
    numElements = 4;
    DistToScreen = 25;
    degreeRadius = 5;
    reps = 200;
    stimLen = 150/1000;
    waitTime = 1.5;
    startPause = 30; % 120 seconds of silence before commencing
    spatFreq = 0.1;
elseif nargin < 3
    DistToScreen = 25;
    degreeRadius = 5;
    reps = 200;
    stimLen = 150/1000;
    waitTime = 1.5;
    startPause = 30; % 120 seconds of silence before commencing
    spatFreq = 0.1;
end

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');
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
% color of 128 = gray with 50% max intensity; 0 = black
[win,~] = Screen('OpenWindow', screenid,0);

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
spatFreq = 1/temp;

distToMass = 135;
centerVals = zeros(numElements,2);
degreeDiv = (2*pi)/numElements;
offset = (2*pi)/3;
for ii=1:numElements
    centerVals(ii,1) = round(centerMass(Channel,1)+cos(degreeDiv*(ii-1)+offset)*distToMass);
    centerVals(ii,2) = round(centerMass(Channel,2)+sin(degreeDiv*(ii-1)+offset)*distToMass);
end
temp = centerVals(2,:);
centerVals(2,:) = centerVals(3,:);
centerVals(3,:) = temp;
 
% for ii=1:4
%     for jj=ii+1:4
%         dist = sqrt((centerVals(ii,1)-centerVals(jj,1)).^2+(centerVals(ii,2)-centerVals(jj,2)).^2)
%     end
% end

estimatedTime = ((stimLen*numElements+waitTime)*reps+3*30+startPause)/60;
display(sprintf('Estimated time: %3.2f minutes',estimatedTime));

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Grey = 0.5;
Black = 0;
White = 1;

%orientation = (rand([numStimuli,1])*360)*pi/180;
% Perform initial flip to gray background and sync us to the retrace:
Priority(9);

usb.startRecording;
WaitSecs(startPause);

% Animation loop
for yy=1:4
    for zz = 1:reps/4
        vbl = Screen('Flip', win);
        for ii=1:numElements
            %         orient = rand*2*pi;
            % Draw the procedural texture as any other texture via 'DrawTexture'
            Screen('DrawTexture', win,gratingTex, [],[],...
                [],[],[],[Grey Grey Grey Grey],...
                [], [],[White,Black,...
                Radius,centerVals(ii,1),centerVals(ii,2),spatFreq,0,0]);
            % Request stimulus onset
            vbl = Screen('Flip', win, vbl + ifi/2);usb.strobe;
            vbl = Screen('Flip',win,vbl+ifi/2+stimLen);
        end
        vbl = Screen('Flip',win,vbl+ifi/2+waitTime);
    end
    if yy ~= 4
        WaitSecs(30);
    end
end
WaitSecs(2);
usb.stopRecording;
Priority(0);

cd('~/Documents/MATLAB/Byron/SeqExp')
fileName = strcat('SeqStim',Date,'_',num2str(AnimalName),'.mat');
save(fileName,'centerVals','Radius','reps','stimLen','startPause',...
    'numElements','w_pixels','h_pixels','spatFreq','mmPerPixel','waitTime')
% Close window
Screen('CloseAll');

end
