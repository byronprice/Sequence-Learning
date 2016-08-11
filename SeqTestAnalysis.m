function [Statistic,Response] = SeqTestAnalysis(AnimalName,Date)
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
%
%Created: 2016/08/10, 24 Cummington, Boston
% Byron Price
%Updated: 2016/08/10
%  By: Byron Price

% read in the .plx file

cd('~/CloudStation/ByronExp/SeqExp/');
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
    temp = smooth(voltage,0.013*sampleFreq);
    n = 30;
    lowpass = 100/(sampleFreq/2); % fraction of Nyquist frequency
    blo = fir1(n,lowpass,'low',hamming(n+1));
    ChanData(:,ii) = filter(blo,1,temp);
end

timeStamps = 0:1/sampleFreq:dataLength/sampleFreq-1/sampleFreq;

if length(timeStamps) ~= dataLength
    display('Error: Review allad cell array and timing')
    return;
end

strobeTimes = tsevs{1,strobeStart};

% strobeDiff = strobeTimes(2:end)-strobeTimes(1:end-1);
% 
% elementStrobes = [];
% for ii=1:length(strobeDiff)
%     if strobeDiff(ii) > (stimTime-0.1) && strobeDiff(ii) < (stimTime+0.1)
%         elementStrobes = [elementStrobes,strobeTimes(ii)];
%     end
% end

% COLLECT LFP RESPONSE TO STIMULI IN ONE MATRIX
stimLen = round(0.2*sampleFreq); % 150ms per sequence element but we'll take 
          % 200 ms because the peak sometimes occurs beyond 150ms
minWin = round(0.04*sampleFreq):1:round(0.1*sampleFreq);
maxWin = round(.1*sampleFreq):1:round(0.2*sampleFreq);

Response = zeros(numChans,numTests,numElements,reps,stimLen);
Statistic = zeros(numChans,numTests,numElements,4);
alpha = 0.05;
for ii=1:numChans
    for jj=1:numTests
        for kk=1:numElements
            elemNum = (jj-1)*(numElements+1)+kk;
            elemStrobes = strobeTimes(svStrobed==elemNum);
            %check = (jj-1)*numElements+1:jj*numElements;
            for ll=1:reps
                stimOnset = elemStrobes(ll);
                [~,index] = min(abs(timeStamps-stimOnset));
                temp = ChanData(index:index+stimLen-1,ii);
                Response(ii,jj,kk,ll,:) = temp;
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
                temp = squeeze(Response(ii,jj,kk,:,:));temp = temp(indeces,:);
                meanResponse = mean(temp,1);
                Tboot(ll) = max(meanResponse(maxWin))-min(meanResponse(minWin));
            end
            meanResponse = mean(squeeze(Response(ii,jj,kk,:,:)),1);
            Statistic(ii,jj,kk,1) = max(meanResponse(maxWin))-min(meanResponse(minWin));Statistic(ii,jj,kk,2) = std(Tboot);
            Statistic(ii,jj,kk,3) = quantile(Tboot,alpha/2);Statistic(ii,jj,kk,4) = quantile(Tboot,1-alpha/2);
        end
    end
%     figure();
%     for ll=1:numElements
%         subplot(ceil(numElements/2),2,ll);histogram(squeeze(Statistic(ii,:,ll)));
%         title(sprintf('Histogram of Response to Sequence Element %d',ll));
%         ylabel('Counts');xlabel( ...
%             sprintf('Test Statistic (\\muV) [max(LFP)-min(LFP) in the %0.2f seconds after sequence element onset]',stimTime));
%     end
    
end
trueStimLen = mean(strobeTimes(svStrobed==2)-strobeTimes(svStrobed==1));
trueStimLen = round(trueStimLen*sampleFreq);
Response = Response(:,:,:,:,1:trueStimLen);

ConvFileName = sprintf('SeqTestConv%d_%d.mat',Date,AnimalName);
save(ConvFileName,'Statistic','Response','Test');

for ii=1:numChans
    figure();plotRows = ceil(numTests/2);
    for jj=1:numTests
        meanRes = [];
        stdRes = [];
        %     lq = [];
        %     uq = [];
        for kk = 1:numElements
            meanRes = [meanRes,mean(squeeze(Response(ii,jj,kk,:,:)),1)];
            stdRes = [stdRes,std(squeeze(Response(ii,jj,kk,:,:)),0,1)];
        end
        stdRes = 2.*stdRes./sqrt(reps);
        subplot(plotRows,2,jj);
        boundedline(1:trueStimLen*numElements,meanRes,stdRes,'alpha');
        title(strcat(sprintf('Mean VEP: Channel %d, Test- ',ii),Test(jj).name));
        ylabel('LFP Voltage (\muV)');xlabel('Time (milliseconds)');
        axis([0 trueStimLen*numElements -500 500]);
    end
end

for ii=1:numChans
    figure();
    plotRows = ceil(numElements/2);
    means = zeros(numTests,numElements);
    stds = zeros(numTests,numElements);
    for jj=1:numTests
        means(jj,:) = squeeze(Statistic(ii,jj,:,1));
        stds(jj,:) = squeeze(Statistic(ii,jj,:,2));
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

