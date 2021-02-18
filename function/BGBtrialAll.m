% BGB trial fun new: variable length

function [trl, event] = BGBtrialAll(cfg)

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
Word = find(strcmp('D113', EVvalue)==1);

% define:
trialNO=0; % trial index value
PostTrig=1000;

% for each word find the condition
for w = 1:length(Word)
    % code for the judgement task: 1 => Aha; 2 => Control;
    if strcmp('D113', EVvalue{Word(w)}) == 1
        if Word(w)+1<=length(EVvalue)
            if strcmp('D251', EVvalue{Word(w)+1}) == 1 || strcmp('D252', EVvalue{Word(w)+1}) == 1
                trialNO=trialNO+1;
                %trialEVsample(trialNO,1)=EVsample(Word(w)+1);
                begsample(trialNO,1)=EVsample(Word(w))-1000;
                endsample(trialNO,1)=EVsample(Word(w)+1)+PostTrig;
                offset(trialNO,1)=EVsample(Word(w))-EVsample(Word(w)+1)-1000;
                RP(trialNO,1)=EVsample(Word(w)+1)-EVsample(Word(w)); % Response time
                if strcmp('D251', EVvalue{Word(w)+1}) == 1
                    task(trialNO,1)=1;
                elseif strcmp('D252', EVvalue{Word(w)+1}) == 1
                    task(trialNO,1)=2;
                end
            end
        end
    end
end

%% the last part is again common to all trial functions
% return the trl matrix (required) and the event structure (optional)
trl = [begsample endsample offset task RP];

end % function