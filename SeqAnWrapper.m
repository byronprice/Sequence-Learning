function [] = SeqAnWrapper(AnimalName)
%SeqAnWrapper.m
%    This wrapper will compare data from the sequence stimulus protocol 
%     across days for the same animal (see SequenceStim.m and readme). 
%     It assumes the data is all in SeqConv files, e.g. SeqConv20160701_34561.mat
%     DailySequences.m is set to automatically run every day on my computer
%     and convert the SeqData ... .plx files to SeqConv ... .mat files. 
%     DailySequences.m calls on SequenceAnalysis.m, which extracts the LFP for 
%     200 ms after each sequence element onset and calculates the VEP
%     magnitude, and its standard error (using bootstrapping). 
%
%INPUT: AnimalName - unique ID for the animal, e.g. 12340.  The first four
%        digits are its cage code and the last is the value of the individual
%        animal (starting at 0 and increasing)

%OUTPUT: figures
%
%Created: 2016/07/12, 24 Cummington Mall
% Byron Price
%Updated: 2016/08/03
%  By: Byron Price

cd('~/CloudStation/ByronExp/SeqExp');

fileStart = sprintf('SeqConv*_%d*',AnimalName);
fileList = dir(fileStart);
numDays = length(fileList);

Stats = cell(numDays,1);
Responses = cell(numDays,1);
for ii = 1:numDays
    load(fileList(ii).name);
    Stats{ii} = Statistic;
    Responses{ii} = Response;
    clear Statistic Response;
end
numChans = size(Responses{1},1);
numElements = size(Responses{1},2);
reps = size(Responses{1},3);
stimLen = size(Responses{1},4);

for ii=1:numChans
    figure();plotRows = ceil(numDays/2);
    for jj=1:numDays
        meanRes = [];
        stdRes = [];
        %     lq = [];
        %     uq = [];
        for zz = 1:numElements
            meanRes = [meanRes,mean(squeeze(Responses{jj}(ii,zz,:,:)),1)];
            stdRes = [stdRes,std(squeeze(Responses{jj}(ii,zz,:,:)),0,1)];
            %         lq = [lq,quantile(squeeze(Response(ii,:,zz,:)),alpha/2,1)];
            %         uq  = [uq,quantile(squeeze(Response(ii,:,zz,:)),1-alpha/2,1)];
        %     lq = meanRes-lq;
        %     uq = uq-meanRes; [lq',uq']
        end
        stdRes = 2*stdRes./sqrt(reps);
        subplot(plotRows,2,jj);
        boundedline(1:stimLen*numElements,meanRes,stdRes,'alpha');
        title(sprintf('Mean VEP with 95%% Confidence Interval: Channel %d, Day %d',ii,jj));
        ylabel('LFP Voltage (\muV)');xlabel('Time (milliseconds)');
        axis([0 stimLen*numElements -500 500]);
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
        title(sprintf('VEP Magnitude with Bootstrap Standard Error for Channel %d, Element %d',ii,kk));
        ylabel('VEP Magnitude (\muV)');xlabel('Experimental Day');
        axis([0 numDays+1 0 500]);
    end
end

% WALD TEST
alpha = 0.05/numElements;
for ii=1:numChans
    for jj=1:numElements
        W = (Stats{numDays}(ii,jj,1)-Stats{1}(ii,jj,1))/(sqrt(Stats{numDays}(ii,jj,2)^2+Stats{1}(ii,jj,2)^2));
        c = norminv(1-alpha,0,1);
        if abs(W) > c
            display(sprintf('Wald Test for Channel %d,Element %d rejects null, mean VEP magnitude greater last day than first day',ii,jj));
        else
            display(sprintf('Wald Test for Channel %d, Element %d retains null',ii,jj))
        end
    end
end

numBases = 50;
Basis = zeros(stimLen*numElements,numBases);
s = 10;
% 50 bases and s = 10 works well, so does 30 and 15
% Gaussian radial basis functions
for ii=1:numBases
    Basis(:,ii) = exp((-((1:stimLen*numElements)-stimLen*numElements*(ii-1)/numBases).^2)./(2*s^2));
%     plot(Basis(:,ii));hold on;
end

Parameters = zeros(numChans,numDays,numBases);
stdErrors = zeros(numChans,numDays,numBases);
estCurve = zeros(numChans,numDays,stimLen*numElements,3);

for ii=1:numChans
    figure();plotRows = ceil(numDays/2);
    for jj=1:numDays
        Y = zeros(stimLen*reps*numElements,1);
        BigBasis = [];
        
        for kk=1:reps
            for ll=1:numElements
                indeces = (ll-1)*stimLen+1:ll*stimLen;
                indeces = indeces+(kk-1)*stimLen*numElements;
                Y(indeces) = squeeze(Responses{jj}(ii,ll,kk,:));
            end
            BigBasis = [BigBasis;Basis];
        end
        
        [b,~,stats] = glmfit(BigBasis,Y,'normal','constant','off');
        [yhat,lBound,uBound] = glmval(b,Basis,'identity',stats,'confidence',1-alpha,'constant','off');
        
        Parameters(ii,jj,:) = b; stdErrors(ii,jj,:) = stats.se;
        estCurve(ii,jj,:,1) = yhat; estCurve(ii,jj,:,2) = lBound; estCurve(ii,jj,:,3) = uBound;
        x = 1:(stimLen*numElements);
        subplot(plotRows,2,jj);boundedline(x,yhat,[lBound,uBound],'alpha');
        title(sprintf('VEP Regression Fit Using %d Gaussian Radial Basis Functions: Channel %d, Day %d',numBases,ii,jj));
        ylabel('LFP Voltage (\muV)');xlabel('Time (milliseconds)');
        axis([0 stimLen*numElements -500 500]);
        %     figure();boundedline(1:numBases,b,2*stats.se,'alpha');
    end
end
end

