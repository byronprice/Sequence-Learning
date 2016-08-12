function [] = SequenceTest(AnimalName,holdTime)
%SequenceTest.m
%  Test after running SequenceStim.m for three days.
% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%
%        Optional- 
%        holdTime - amount of time to wait between blocks of about 50 stimuli
% 
%        see file SequenceVars.mat for other changeable presets
%
% OUTPUT: a file with stimulus parameters named SeqTestDate_AnimalName
%           e.g. SeqTest20160708_12345.mat to be saved in CloudStation's SeqExp
%           folder under '~/CloudStation/ByronExp/SeqExp'
% Created: 2016/08/04 at 5920 Colchester Road, Fairfax, VA
%  Byron Price
% Updated: 2016/08/10
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

% usb = ttlInterfaceClass.getTTLInterface;
% if ~usb.validateInterface
%     handleWarning('TTL Interface validation failed',true,'TTL Warning',[],true);
% end
% usb.stopRecording; % for some reason, doesn't start unless previously stopped

% sio = screenInterfaceClass.returnInterface;
% sio.openScreen;
% win = sio.window;
% ifi = sio.slack*2;
% mp = sio.getMonitorProfile;
% w_pixels = mp.cols;
% h_pixels = mp.rows;
% w_mm = 1e3*mp.screen_width;
% h_mm = 1e3*mp.screen_height;

% % Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));
% 
% % Open a fullscreen onscreen window on that display, choose a background
% % color of 127 = gray with 50% max intensity; 0 = black
background = 127; % gray, mean luminance
[win,~] = Screen('OpenWindow', screenid,background);

% Switch color specification to use the 0.0 - 1.0 range
Screen('ColorRange', win, 1);

% % Query window size in pixels
[w_pixels, h_pixels] = Screen('WindowSize', win);
% 
% % Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);

dgshader = [directory '/SequenceTest.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [directory '/SequenceTest.frag.txt'] }, 1);
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

numTests = 4;
% extract the angle of the 2-D Gaussian to target the stimulus to the
% strongest responding locations within the retinotopic map
[V,D] = eig(squeeze(Sigma(Channel,:,:)));
[~,index] = max(max(D));
offset = atan2(V(2,index),V(1,index));

if offset < 0 
    offset = offset+2*pi;
end

centerVals = zeros(numTests,numElements,2);
degreeDiv = (2*pi)/numElements;
centerVals(1,1,1) = centerMass(Channel,1);centerVals(1,1,2) = centerMass(Channel,2);
for ii=2:numElements
    centerVals(1,ii,1) = centerMass(Channel,1)+cos(degreeDiv*(ii-2)+offset)*2*Radius;
    centerVals(1,ii,2) = centerMass(Channel,2)+sin(degreeDiv*(ii-2)+offset)*2*Radius;
end
temp2 = centerVals(1,2,:);temp3 = centerVals(1,3,:);
centerVals(1,2,:) = centerVals(1,1,:);
centerVals(1,3,:) = temp2;
centerVals(1,1,:) = temp3;

% first test, original sequence
% second test, original sequence, but missing second element
% third test, reversed positions, original orientation order
% fourth test, reveresed orientations, original position order
alpha = ones(numTests,numElements);
alpha(2,2) = 0; % alpha mixing for second test, second element

orient = zeros(numTests,numElements);
temp = [150,60,30,120].*pi./180;
orient(1,:) = temp;
orient(2,:) = temp;
orient(3,:) = temp;
orient(4,:) = fliplr(temp);
centerVals(2,:,1) = centerVals(1,:,1);
centerVals(2,:,2) = centerVals(1,:,2);
centerVals(3,:,1) = fliplr(centerVals(1,:,1));
centerVals(3,:,2) = fliplr(centerVals(1,:,2));
centerVals(4,:,1) = centerVals(1,:,1);
centerVals(4,:,2) = centerVals(1,:,2);

estimatedTime = (numTests*((stimTime*numElements+waitTime)*reps+blocks*holdTime))/60;
display(sprintf('\nEstimated time: %3.2f minutes',estimatedTime));

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Grey = 0.5;

Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

% Perform initial flip to gray background and sync us to the retrace:
Priority(9);

elemNums = zeros(numTests,numElements);
for ii=1:numTests
    for jj=1:numElements
    mystr = sprintf('%d%d',ii,jj);
    elemNums(ii,jj) = str2double(mystr);
    end
end

usb.startRecording;
WaitSecs(1);
usb.strobeEventWord(0);
WaitSecs(holdTime);

% Animation loop
for ii=1:numTests
    for jj=1:blocks
        vbl = Screen('Flip',win);
        for kk= 1:reps/blocks
            for ll=1:numElements
                % Draw the procedural texture as any other texture via 'DrawTexture'
                Screen('DrawTexture', win,gratingTex,[],[],...
                    [],[],[],[Grey Grey Grey Grey],...
                    [], [],[alpha(ii,ll),0,...
                    Radius,centerVals(ii,ll,1),centerVals(ii,ll,2),spatFreq,orient(ii,ll),gama]);
                % Request stimulus onset
                vbl = Screen('Flip', win,vbl+ifi/2);usb.strobeEventWord(elemNums(ii,ll));
                vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
            end
            greyNum = ii;
            usb.strobeEventWord(greyNum);
            vbl = Screen('Flip',win,vbl-ifi/2+waitTime);
        end
        usb.strobeEventWord(0);
        vbl = Screen('Flip',win,vbl-ifi/2+holdTime);
    end
end
WaitSecs(2);
usb.stopRecording;
Priority(0);

Test = struct('name',cell(numTests,1));
Test(1).name = 'Original';
Test(2).name = 'Blank second element';
Test(3).name = 'Same orientations, reversed positions';
Test(4).name = 'Same positions, reversed orientations';

cd('~/CloudStation/ByronExp/SeqExp');
fileName = sprintf('SeqTest%d_%d.mat',Date,AnimalName);
save(fileName,'centerVals','Radius','reps','stimTime','numElements',...
    'w_pixels','h_pixels','spatFreq','mmPerPixel','waitTime','holdTime',...
    'DistToScreen','numTests','Test','orient')
% Close window
Screen('CloseAll');

end