function [Stat,Data,Latency] = SeqTestAnalysis(AnimalName,Date)
%SeqTestAnalysis.m
%   Analyze data from the test day, after 3 days of experiments.
%   SequenceStim.m is shown to the mice on three consecutive days, with
%   analysis performed by SeqStimAnalysis.m and SeqAnWrapper.m . Here, the
%   stimulus file is SequenceTest.m . 
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
%        Latency - time (in seconds) to 
%
%Created: 2016/08/10, 24 Cummington, Boston
% Byron Price
%Updated: 2016/08/17
%  By: Byron Price

% read in the .plx file

cd('~/CloudStation/ByronExp/Seq/');
EphysFileName = sprintf('SeqTestData%d_%d',Date,AnimalName); % don't want the
                      % file identifier at this point as MyReadall.m does
                      % that part

if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
    MyReadall(EphysFileName);
end

StimulusFileName = sprintf('SeqTest%d_%d.mat',Date,AnimalName);
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
          % 250 ms because the peak sometimes occurs beyond 150ms
minStart = round(0.05*sampleFreq);minEnd = round(0.15*sampleFreq);
maxStart = round(.1*sampleFreq);maxEnd = round(0.225*sampleFreq);
minWin = minStart:1:minEnd;
maxWin = maxStart:1:maxEnd;
smoothKernel = 4;

Data = struct;
Data.VEP = zeros(numChans,numTests,numElements,reps,stimLen);
Data.meanVEP = zeros(numChans,numTests,numElements,stimLen);
Latency = struct;
Latency.minVal = zeros(numChans,numTests,numElements);
Latency.minTime = zeros(numChans,numTests,numElements);
Latency.maxVal = zeros(numChans,numTests,numElements);
Latency.maxTime = zeros(numChans,numTests,numElements);
Stat = struct;
Stat.VEPsize = zeros(numChans,numTests,numElements);
Stat.VEPsem = zeros(numChans,numTests,numElements);
Stat.lbound = zeros(numChans,numTests,numElements);
Stat.ubound = zeros(numChans,numTests,numElements);
alpha = 0.05;
for ii=1:numChans
    for jj=1:numTests
        for kk=1:numElements
            mystr = sprintf('%d%d',jj,kk);
            elemNum = str2double(mystr);
            elemStrobes = strobeTimes(svStrobed==elemNum);
            %check = (jj-1)*numElements+1:jj*numElements;
            for ll=1:reps
                stimOnset = elemStrobes(ll);
                [~,index] = min(abs(timeStamps-stimOnset));
                temp = ChanData(index:index+stimLen-1,ii);
                Data.VEP(ii,jj,kk,ll,:) = temp;
            end
            clear check temp;
        end
    end
    % BOOTSTRAP ERROR ON TEST STATISTIC
    n = 5000;
    for jj=1:numTests
        for kk=1:numElements
            Tboot = zeros(n,1);
            for ll=1:n
                indeces = random('Discrete Uniform',reps,[reps,1]);
                temp = squeeze(Data.VEP(ii,jj,kk,:,:));temp = temp(indeces,:);
                meanResponse = mean(temp,1);
                Tboot(ll) = max(meanResponse(maxWin))-min(meanResponse(minWin));
            end
            meanResponse = mean(squeeze(Data.VEP(ii,jj,kk,:,:)),1);
            meanResponse = smooth(meanResponse,smoothKernel);
            Data.meanVEP(ii,jj,kk,:) = meanResponse;
            Stat.VEPsize(ii,jj,kk) = max(meanResponse(maxWin))-min(meanResponse(minWin));Stat.VEPsem(ii,jj,kk) = std(Tboot);
            Stat.lbound(ii,jj,kk) = quantile(Tboot,alpha/2);Stat.ubound(ii,jj,kk) = quantile(Tboot,1-alpha/2);
            [minVal,minInd] = min(meanResponse(minWin));Latency.minTime(ii,jj,kk) = minStart/sampleFreq+minInd/sampleFreq;
            Latency.minVal(ii,jj,kk) = minVal;
            [maxVal,maxInd] = max(meanResponse(maxWin));Latency.maxTime(ii,jj,kk) = maxStart/sampleFreq+maxInd/sampleFreq;
            Latency.maxVal(ii,jj,kk) = maxVal;
            clear meanResponse;
        end
    end
end
trueStimLen = mean(strobeTimes(svStrobed==2)-strobeTimes(svStrobed==1));
trueStimLen = round(trueStimLen*sampleFreq);
Data.VEP = Data.VEP(:,:,:,:,1:trueStimLen);

ConvFileName = sprintf('SeqTestConv%d_%d.mat',Date,AnimalName);
save(ConvFileName,'Stat','Data','Test','Latency','sampleFreq');

latStart = 0:trueStimLen:trueStimLen*(numElements-1);
for ii=1:numChans
    figure();plotRows = ceil(numTests/2);
    for jj=1:numTests
        meanRes = [];
        stdRes = [];
        %     lq = [];
        %     uq = [];
        for kk = 1:numElements
            meanRes = [meanRes,squeeze(Data.meanVEP(ii,jj,kk,:))];
            stdRes = [stdRes,std(squeeze(Data.VEP(ii,jj,kk,:,:)),0,1)];
        end
        stdRes = 2.*stdRes./sqrt(reps);
        subplot(plotRows,2,jj);
        boundedline(1:trueStimLen*numElements,meanRes,stdRes,'alpha');
        title(strcat(sprintf('Mean VEP: Channel %d, Test- ',ii),Test(jj).name));
        ylabel('LFP Voltage (\muV)');xlabel('Time (milliseconds)');
        axis([0 trueStimLen*numElements -500 500]);
        hold on; plot(latStart'+squeeze(Latency.minTime(ii,jj,:))*sampleFreq,squeeze(Latency.minVal(ii,jj,:)),'vr');
        plot(latStart'+squeeze(Latency.maxTime(ii,jj,:))*sampleFreq,squeeze(Latency.maxVal(ii,jj,:)),'^k');
        hold off;
    end
end

for ii=1:numChans
    figure();
    plotRows = ceil(numElements/2);
    means = zeros(numTests,numElements);
    stds = zeros(numTests,numElements);
    for jj=1:numTests
        means(jj,:) = squeeze(Stat.VEPsize(ii,jj,:));
        stds(jj,:) = squeeze(Stat.VEPsem(ii,jj,:));
    end
    for kk=1:numElements
        subplot(plotRows,2,kk);
        errorbar(1:numTests,means(:,kk),...
            stds(:,kk),'LineWidth',2);
        title(sprintf('VEP Magnitude with Bootstrap Standard Error for Channel %d, Element %d',ii,kk));
        ylabel('VEP Magnitude (\muV)');xlabel('Test #');
        axis([0 numTests+1 0 500]);
    end
end

end

