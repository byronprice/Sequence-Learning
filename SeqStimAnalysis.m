function [Stat,Data,Latency,sampleFreq ] = SeqStimAnalysis(AnimalName,Date)
%SeqStimAnalysis.m
%   Analyze data from one day of experiments for sequence learning.  Mice
%   are shown ~200 sequences of four elements. Each element is an oriented
%   grating that displays for 150 ms. The elements are sinusoidal grating with a
%   a circular 2D Gaussian kernel overlying to create a circular image that
%   decays at the edges. Each circle occupies about a 5-degree radius of visual
%   space (see SequenceStim.m).  The stimulus lasts a
%   total of 150*4 = 600 ms, followed by a blank grey screen for 1.5
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
%OUTPUT: Statistic - stats on VEP magnitude in response to sequence
%         elements
%        Response - the LFP signal of the response to sequence elements
%        Latency - latency of minimum and maximum mean LFP after each
%         element onset (in seconds)
%        sampleFreq - sampling frequency
%
%Created: 2016/07/11
% Byron Price
%Updated: 2016/08/17
%  By: Byron Price

% read in the .plx file

cd('~/CloudStation/ByronExp/Seq/');
EphysFileName = sprintf('SeqData%d_%d',Date,AnimalName); % don't want the
                      % file identifier at this point as MyReadall.m does
                      % that part

if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
    MyReadall(EphysFileName);
end

StimulusFileName = sprintf('SeqStim%d_%d.mat',Date,AnimalName);
EphysFileName = strcat(EphysFileName,'.mat');
load(EphysFileName)
load(StimulusFileName)

sampleFreq = adfreq;
strobeStart = 33;

Chans = find(~cellfun(@isempty,allad));numChans = length(Chans);

% convert allad data to millivolts, then lowpass filter the data
dataLength = length(allad{1,Chans(1)});
ChanData = zeros(dataLength,numChans);
preAmpGain = 1;
for ii=1:numChans
    voltage = 1000.*((allad{1,Chans(ii)}).*SlowPeakV)./(0.5*(2^SlowADResBits)*adgains(Chans(ii))*preAmpGain);
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

strobeTimes = tsevs{1,strobeStart};

% COLLECT LFP RESPONSE TO STIMULI IN ONE MATRIX
stimLen = round((stimTime+0.1)*sampleFreq); % 150ms per sequence element but we'll take 
          % 200 ms because the peak sometimes occurs beyond 150ms
minStart = round(0.05*sampleFreq);minEnd = round(0.15*sampleFreq);
maxStart = round(.1*sampleFreq);maxEnd = round(0.225*sampleFreq);
minWin = minStart:1:minEnd;
maxWin = maxStart:1:maxEnd;
smoothKernel = 4;

Data = struct;
Data.VEP = zeros(numChans,numElements,reps,stimLen);
Data.meanVEP = zeros(numChans,numElements,stimLen);
Latency = struct;
Latency.minVal = zeros(numChans,numElements);
Latency.minTime = zeros(numChans,numElements);
Latency.maxVal = zeros(numChans,numElements);
Latency.maxTime = zeros(numChans,numElements);
Stat = struct;
Stat.VEPsize = zeros(numChans,numElements);
Stat.VEPsem = zeros(numChans,numElements);
Stat.lbound = zeros(numChans,numElements);
Stat.ubound = zeros(numChans,numElements);

alpha = 0.05;
for ii=1:numChans
    for jj=1:numElements
        elemStrobes = strobeTimes(svStrobed==jj);
        %check = (jj-1)*numElements+1:jj*numElements;
        for kk=1:reps
            stimOnset = elemStrobes(kk);
            [~,index] = min(abs(timeStamps-stimOnset));
            temp = ChanData(index:index+stimLen-1,ii);
            Data.VEP(ii,jj,kk,:) = temp;
        end
        clear check temp;
    end
    
    % BOOTSTRAP ERROR ON TEST STATISTIC
    for jj=1:numElements
        n = 5000;
        Tboot = zeros(n,1);
        for kk=1:n
            indeces = random('Discrete Uniform',reps,[reps,1]);
            temp = squeeze(Data.VEP(ii,jj,:,:));temp = temp(indeces,:);
            meanResponse = mean(temp,1);
            Tboot(kk) = max(meanResponse(maxWin))-min(meanResponse(minWin)); 
        end
        meanResponse = mean(squeeze(Data.VEP(ii,jj,:,:)),1);
        meanResponse = smooth(meanResponse,smoothKernel);
        Data.meanVEP(ii,jj,:) = meanResponse;
        Stat.VEPsize(ii,jj) = max(meanResponse(maxWin))-min(meanResponse(minWin));Stat.VEPsem(ii,jj) = std(Tboot);
        Stat.lbound(ii,jj) = quantile(Tboot,alpha/2);Stat.ubound(ii,jj) = quantile(Tboot,1-alpha/2);
        [minVal,minInd] = min(meanResponse(minWin));Latency.minTime(ii,jj) = minStart/sampleFreq+minInd/sampleFreq;
        Latency.minVal(ii,jj) = minVal;
        [maxVal,maxInd] = max(meanResponse(maxWin));Latency.maxTime(ii,jj) = maxStart/sampleFreq+maxInd/sampleFreq;
        Latency.maxVal(ii,jj) = maxVal;
        clear meanResponse;
    end
end

trueStimLen = mean(strobeTimes(svStrobed==2)-strobeTimes(svStrobed==1));
trueStimLen = round(trueStimLen*sampleFreq);
Data.VEP = Data.VEP(:,:,:,1:trueStimLen);
Data.meanVEP = Data.meanVEP(:,:,1:trueStimLen);
end

