% BGB trial function for D111: variable length

function [trl, event] = BGBtrialD111(cfg)

%% the first part is common to all trial functions
% read the header (needed for the samping rate) and the events
hdr        = ft_read_header(cfg.headerfile);
event      = ft_read_event(cfg.headerfile);

%% from here on it becomes specific to the experiment and the data format
% for the events of interest, find the sample numbers (these are integers)
% for the events of interest, find the trigger values (these are strings in the case of BrainVision)
EVsample   = [event.sample]';
EVvalue    = {event.value}';

% select the target stimuli
Word = find(strcmp('D111', EVvalue)==1);
% define:
trialNO=0; % trial index value
PostTrig=1000;
% for each word find the condition
begsample=[];
endsample=[];
offset = [];
RP =[];

%load data;
newSub =cfg.headerfile(1:3); % to get the sub code/ cannot use "Subject" directly
loaddata = ['load  /mnt/v7k/home/luy/BGBdata/forData/' newSub '/data.mat']; % here can be "data" or "dataD113"
eval(loaddata);
oriTr=data.cfg.trl(:,6);
oriRT=data.cfg.trl(:,5);
task=data.cfg.trl(:,4);
for i = 1:length(oriTr)
    trialNO=trialNO+1;
    if oriTr(i)<35 % be careful, may change here (depend on subject) % to deal with extremum
        begsample(trialNO,1)=EVsample(Word(oriTr(i)))-1000;
        endsample(trialNO,1)=EVsample(Word(oriTr(i)))+ oriRT(i) + PostTrig;
    else
        begsample(trialNO,1)=EVsample(Word(oriTr(i)-1))-1000;
        endsample(trialNO,1)=EVsample(Word(oriTr(i)-1))+ oriRT(i) + PostTrig;
    end
    offset(trialNO,1)= - oriRT(i) -1000;
    RP(trialNO,1)=oriRT(i) ; % Response time
    
end

%% the last part is again common to all trial functions
% return the trl matrix (required) and the event structure (optional)
trl = [begsample endsample offset task RP data.cfg.trl(:,6) data.cfg.trl(:,7)]; % original code
end % function