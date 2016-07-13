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

stats = zeros(numDays,numChans,4);
count = 1;
for ii = 1:numDays
    [Statistic,Parameters,stdErrors,estCurve,stimLength] = SequenceAnalysis(...
        AnimalName,Days(ii));
    stats(count,:,:) = Statistic;
    count = count+1;
end

for ii=1:numChans
    figure();errorbar(1:numDays,abs(squeeze(stats(:,ii,1))),...
        squeeze(stats(:,ii,2))./2,'LineWidth',2);
    title('Test Statistic [abs(min of mean VEP after first element)] with Bootstrap Standard Error');
    ylabel('VEP Magnitude (\muV)');xlabel('Experimental Day');
    axis([0 numDays+1 40 100]);
end
end

