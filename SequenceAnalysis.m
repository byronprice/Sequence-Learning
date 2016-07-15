function [Statistic,Response] = SequenceAnalysis(AnimalName,Date,statFun,Chans,Day)
%SequenceAnalysis.m
%   Analyze data from one day of experiments for sequence learning.  Mice
%   are shown ~200 sequences of four elements. Each element is an oriented
%   grating that displays for 160 ms.  The stimulus lasts, therefore, a
%   total of 160*4 = 640 ms, followed by a blank grey screen for 1.5
%   seconds and then a repetition of the four-element sequence.
%
% Reference: Gavornik & Bear, 2014

%INPUT: AnimalName - unique ID for the animal, e.g. 12340.  The first four
%        digits are its cage code and the last is the value of the individual
%        animal (starting at 0 and increasing)
%       Date - experimental date, e.g. 20160711 for July 11, 2016, always
%        use two digits for days and months, i.e. June is 06, the ninth day
%        of June is 09
%
%        OPTIONAL
%       Chans - a vector of channel numbers used to record neural activity,
%        defaults to channels 6 and 8, e.g. [6,8]
%       Day - which day of the experiment as a number
%        (they usually run 5 days), e.g. 1
%OUTPUT: Statistic -
%        Parameters - 
%        stdErrors - 
%        estCurve - 
%        stimLength -
%
%Created: 2016/07/11
% Byron Price
%Updated: 2016/07/13
%  By: Byron Price

% read in the .plx file

cd('/Users/byronprice/Documents/Current-Projects/ExperimentData/');
EphysFileName = strcat('SeqData',num2str(Date),'_',num2str(AnimalName));

if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
    MyReadall(EphysFileName);
end

%StimulusFileName = strcat('SequenceStim',num2str(Date),'_',num2str(AnimalName),'.mat');
EphysFileName = strcat(EphysFileName,'.mat');
load(EphysFileName)

if nargin < 3
    statFun = @(x) abs(min(x));
    Chans = [6,8];
    Day = 1;
elseif nargin < 4
    Chans = [6,8];
    Day = 1;
elseif nargin < 5
    Day = 1;
end

stimTime = 0.167;
sampleFreq = adfreq;
strobeStart = 33;

% convert allad data to millivolts, then lowpass filter the data
dataLength = length(allad{1,strobeStart+Chans(1)-1});
numChans = length(Chans);
ChanData = zeros(dataLength,numChans);
preAmpGain = 1/1000;

for ii=1:numChans
    voltage = ((allad{1,strobeStart+Chans(ii)-1}).*SlowPeakV)./(0.5*(2^SlowADResBits)*adgains(strobeStart+Chans(ii)-1)*preAmpGain);
    n = 30;
    lowpass = 100/(sampleFreq/2); % fraction of Nyquist frequency
    blo = fir1(n,lowpass,'low',hamming(n+1));
    ChanData(:,ii) = filter(blo,1,voltage);
end

timeStamps = 0:1/sampleFreq:dataLength/sampleFreq-1/sampleFreq;

if length(timeStamps) ~= dataLength
    display('Error: Review allad cell array and timing')
    return;
end
strobeData = tsevs{1,strobeStart};
strobeDiff = strobeData(2:end)-strobeData(1:end-1);

elementStrobes = [];
for ii=1:length(strobeDiff)
    if strobeDiff(ii) > (stimTime-0.1) && strobeDiff(ii) < (stimTime+0.1)
        elementStrobes = [elementStrobes,strobeData(ii)];
    end
end

% COLLECT LFP RESPONSE TO STIMULI IN ONE MATRIX
numElements = 4;
numStimuli = length(elementStrobes)/numElements;

stimLength = round(stimTime*sampleFreq); % 167ms per sequence element

Response = zeros(numChans,numStimuli,numElements,stimLength);
Statistic = zeros(numChans,numElements,4);
alpha = 0.05;
for ii=1:numChans
    for jj=1:numStimuli
        check = (jj-1)*numElements+1:jj*numElements;
        for kk=1:numElements
            stimOnset = elementStrobes(check(kk));
            [~,index] = min(abs(timeStamps-stimOnset));
            temp = ChanData(index:index+stimLength-1,ii);
            Response(ii,jj,kk,:) = temp;
        end
        clear check temp;
    end
    
    % BOOTSTRAP ERROR ON TEST STATISTIC
    for jj=1:numElements
        n = 5000;
        Tboot = zeros(n,1);
        for kk=1:n
            indeces = random('Discrete Uniform',numStimuli,[numStimuli,1]);
            temp = squeeze(Response(ii,:,jj,:));temp = temp(indeces,:);
            meanResponse = mean(temp,1);
            Tboot(kk) = statFun(meanResponse); 
        end
        Statistic(ii,jj,1) = statFun(mean(squeeze(Response(ii,:,jj,:)),1));Statistic(ii,jj,2) = std(Tboot);
        Statistic(ii,jj,3) = quantile(Tboot,alpha/2);Statistic(ii,jj,4) = quantile(Tboot,1-alpha/2);
    end
%     figure();
%     for ll=1:numElements
%         subplot(ceil(numElements/2),2,ll);histogram(squeeze(Statistic(ii,:,ll)));
%         title(sprintf('Histogram of Response to Sequence Element %d',ll));
%         ylabel('Counts');xlabel( ...
%             sprintf('Test Statistic (\\muV) [max(LFP)-min(LFP) in the %0.2f seconds after sequence element onset]',stimTime));
%     end
    
end


end

