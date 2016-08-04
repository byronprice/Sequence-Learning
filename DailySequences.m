%DailySequences.m
%   Used along with a shell script to automatically run all the animals
%    from the current day who have seen the sequence stimulus
%
%INPUT: Will run through all of the files in the CloudStation folder
%        ByronExp/SeqExp and see if any new files exist, if they do, it
%        will run SeqAnWrapper.m on those
%OUTPUT:
%
% Created: 2016/08/02, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/08/03
%  By: Byron Price 

cd('~/CloudStation/ByronExp/SeqExp');

today = datetime('today','Format','yyyy-MM-dd');
today = char(today); today = strrep(today,'-','');
today = str2double(today);

fileStart = 'SeqData*.plx';

fileList = dir(fileStart);
numFiles = size(fileList,1);

datelen = 8;
idlen = 5;
for ii=1:numFiles
    index = regexp(fileList(ii).name,'_');
    Date = str2double(fileList(ii).name(index-datelen:index-1));
    if Date == today
        AnimalName = str2double(fileList(ii).name(index+1:index+idlen));
        [Statistic,Response] = SequenceAnalysis(AnimalName,Date);
        save(sprintf('SeqConv%d_%d.mat',Date,AnimalName),...
            'Statistic','Response');
    end
    clear Statistic Response;
end
