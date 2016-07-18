function [] = SeqAnWrapper(AnimalName,firstDay,numDays,Chans)
%SeqAnWrapper.m
%   SequenceAnalysis.m runs functions on data from individual days, but
%    this wrapper will combine data from multiple days
%
%
%INPUT: AnimalName - unique ID for the animal, e.g. 12340.  The first four
%        digits are its cage code and the last is the value of the individual
%        animal (starting at 0 and increasing)
%       firstDay - experimental start day, e.g. 20160711 for July 11, 2016, always
%        use two digits for days and months, i.e. June is 06, the ninth day
%        of June is 09, 
%       numDays - the number of days of the experiment, will assume
%           consecutive days
%        OPTIONAL
%       Chans - a vector of channel numbers used to record neural activity,
%        defaults to channels 6 and 8, e.g. [6,8]
%OUTPUT:
%
%Created: 2016/07/12
% Byron Price
%Updated: 2016/07/13
%  By: Byron Price

if nargin < 4
    Chans = [6,8];
end
Days = zeros(numDays,1);
Days(1) = firstDay;
firstDay = num2str(firstDay);
x = datenum(str2double(firstDay(1:4)),str2double(firstDay(5:6)),str2double(firstDay(7:8)));

for ii=2:numDays
    x = x+1;
    result = datetime(x,'ConvertFrom','datenum','Format','yyyyMMdd');
    Days(ii) = str2double(datestr(result,'yyyymmdd'));
end

numChans = length(Chans);

%statFun = @(x) trapz(abs(x));
%statFun = @(x) abs(min(x));
statFun = @(x) max(x)-min(x);

Stats = cell(numDays,1);
Responses = cell(numDays,1);
count = 1;
for ii = 1:numDays
    [Statistic,Response] = SequenceAnalysis(...
        AnimalName,Days(ii),statFun,Chans,ii);
    Stats{ii} = Statistic;
    Responses{ii} = Response;
    count = count+1;
end
numStimuli = size(Responses{1},2);
numElements = size(Responses{1},3);
stimLength = size(Responses{1},4);

for ii=1:numChans
    figure();plotRows = ceil(numDays/2);
    for jj=1:numDays
        meanRes = [];
        stdRes = [];
        %     lq = [];
        %     uq = [];
        for zz = 1:numElements
            meanRes = [meanRes,mean(squeeze(Responses{jj}(ii,:,zz,:)),1)];
            stdRes = [stdRes,std(squeeze(Responses{jj}(ii,:,zz,:)),0,1)];
            %         lq = [lq,quantile(squeeze(Response(ii,:,zz,:)),alpha/2,1)];
            %         uq  = [uq,quantile(squeeze(Response(ii,:,zz,:)),1-alpha/2,1)];
        %     lq = meanRes-lq;
        %     uq = uq-meanRes; [lq',uq']
        end
        stdRes = 2*stdRes./sqrt(numStimuli);
        subplot(plotRows,2,jj);
        boundedline(1:stimLength*numElements,meanRes,stdRes,'alpha');
        title(sprintf('Mean VEP with 95%% Confidence Interval: Channel %d, Day %d',Chans(ii),jj));
        ylabel('LFP Voltage (\muV)');xlabel('Time (milliseconds)');
        axis([0 stimLength*numElements -250 250]);
    end
end

for ii=1:numChans
    figure();
    plotRows = ceil(numElements/2);
    means = zeros(numDays,numElements);
    stds = zeros(numDays,numElements);
    for jj=1:numDays
        means(jj,:) = squeeze(Stats{jj}(ii,:,1));
        stds(jj,:) = squeeze(Stats{jj}(ii,:,2));
    end
    for kk=1:numElements
        subplot(plotRows,2,kk);
        errorbar(1:numDays,means(:,kk),...
            stds(:,kk),'LineWidth',2);
        title(sprintf('Test Statistic (measure of VEP magnitude) with Bootstrap Standard Error for Channel %d, Element %d',Chans(ii),kk));
        ylabel('VEP Magnitude (\muV)');xlabel('Experimental Day');
        axis([0 numDays+1 0 300]);
    end
end

% WALD TEST
alpha = 0.05/numElements;
for ii=1:numChans
    for jj=1:numElements
        W = (Stats{numDays}(ii,jj,1)-Stats{1}(ii,jj,1))/(sqrt(Stats{numDays}(ii,jj,2)^2+Stats{1}(ii,jj,2)^2));
        c = norminv(1-alpha,0,1);
        if abs(W) > c
            display(sprintf('Wald Test for Channel %d,Element %d rejects null, mean VEP magnitude differs between first and last day',Chans(ii),jj));
        else
            display(sprintf('Wald Test for Channel %d, Element %d retains null',Chans(ii),jj))
        end
    end
end

numBases = 50;
Basis = zeros(stimLength*numElements,numBases);
s = 10;
% 50 bases and s = 10 works well, so does 30 and 15
% Gaussian radial basis functions
for ii=1:numBases
    Basis(:,ii) = exp((-((1:stimLength*numElements)-stimLength*numElements*(ii-1)/numBases).^2)./(2*s^2));
%     plot(Basis(:,ii));hold on;
end

Parameters = zeros(numChans,numDays,numBases);
stdErrors = zeros(numChans,numDays,numBases);
estCurve = zeros(numChans,numDays,stimLength*numElements,3);

for ii=1:numChans
    figure();plotRows = ceil(numDays/2);
    for jj=1:numDays
        Y = zeros(stimLength*numStimuli*numElements,1);
        BigBasis = [];
        
        for kk=1:numStimuli
            for ll=1:numElements
                indeces = (ll-1)*stimLength+1:ll*stimLength;
                indeces = indeces+(kk-1)*stimLength*numElements;
                Y(indeces) = squeeze(Responses{jj}(ii,kk,ll,:));
            end
            BigBasis = [BigBasis;Basis];
        end
        
        [b,~,stats] = glmfit(BigBasis,Y,'normal','constant','off');
        [yhat,lBound,uBound] = glmval(b,Basis,'identity',stats,'confidence',1-alpha,'constant','off');
        
        Parameters(ii,jj,:) = b; stdErrors(ii,jj,:) = stats.se;
        estCurve(ii,jj,:,1) = yhat; estCurve(ii,jj,:,2) = lBound; estCurve(ii,jj,:,3) = uBound;
        x = 1:(stimLength*numElements);
        subplot(plotRows,2,jj);boundedline(x,yhat,[lBound,uBound],'alpha');
        title(sprintf('VEP Regression Fit Using %d Gaussian Radial Basis Functions: Channel %d, Day %d',numBases,Chans(ii),jj));
        ylabel('LFP Voltage (\muV)');xlabel('Time (milliseconds)');
        axis([0 stimLength*numElements -250 250]);
        %     figure();boundedline(1:numBases,b,2*stats.se,'alpha');
    end
end
end

