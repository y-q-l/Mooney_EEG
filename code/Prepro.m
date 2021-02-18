%% Analysis based on FieldTrip

%% Pre-processing
filename = {......}; % use your own subjects code
% read the file name:
list = dir('R*.*'); % read file names with  "R......"   % use your own path /EEG data
trialcheck = [];
for sub = 1:length(filename)
Subject = filename{sub};
cfg = [];
cfg.headerfile = list(sub,1).name;
headerfile{sub,1} = list(sub,1).name; % headerfile is R*.*
cfg.trialfun = 'BGBtrialAll'; % private function  % Using BGBtrialD111 for the 1st stage data
cfg = ft_definetrial(cfg);
loaddata = ['load ....../' Subject '_piclistreact_new.mat']; % use your own path / behavioral data
eval(loaddata);
cfg.trl = [cfg.trl piclistreact(:,7:8)];
% remove the trials that don't consistent: 251~1 or 252~2
wrong=[];
for tr=1:length(cfg.trl)
    if cfg.trl(tr,7)==0 || cfg.trl(tr,5)>7000
        wrong=[wrong tr];
    end
end
cfg.trl(wrong,:) = [];
cfg.channel = {'all'}; % and use '-Exx' to remove bad channels
cfg.reref = 'yes';
cfg.refchannel = {'all'};
cfg.demean = 'no';
dataD113 = ft_preprocessing(cfg);  % D113 represents the 3rd-stage epochs / D111 represents the 1st-stage epochs
trialcheck(sub,1) = length(dataD113.cfg.trl);
savedata = ['save ....../' Subject '/dataD113' '.mat' ' dataD113']; eval(savedata); % use your own path
clearvars -except filename sub headerfile list trialcheck;
end

%% Separate RP and OS (for D111 the principle conventions as for D113)
%% For RP
load dataD113;
dataD113_rp=dataD113;
dataD113_rp.time=[];
dataD113_rp.trial=[];
for i=1:length(dataD113.time)
    dataD113_rp.time{1,i}=dataD113.time{1,i}(1,end-2000:end-500);
end
for j=1:length(dataD113.trial)
    dataD113_rp.trial{1,j}=dataD113.trial{1,j}(:,end-2000:end-500);
end
dataD113_rp.sampleinfo(:,1)=dataD113_rp.sampleinfo(:,2)-2000;
dataD113_rp.sampleinfo(:,2)=dataD113_rp.sampleinfo(:,2)-500;
save dataD113_rp dataD113_rp;
%% For OS:
dataD113_os=dataD113;
dataD113_os.time=[];
dataD113_os.trial=[];
for i=1:length(dataD113.time)
    dataD113_os.time{1,i}=dataD113.time{1,i}(1,1+500:1+2000);
    dataD113_os.time{1,i}=dataD113_os.time{1,i}+ (length(dataD113.time{1,i})-2001)/1000  ; % i.e., add RT, so the time axis of RP moves to OS
end
for j=1:length(dataD113.trial)
    dataD113_os.trial{1,j}=dataD113.trial{1,j}(:,1+500:1+2000);
end
dataD113_os.sampleinfo(:,2)=dataD113_os.sampleinfo(:,1)+2000;
dataD113_os.sampleinfo(:,1)=dataD113_os.sampleinfo(:,1)+500;
save  dataD113_os dataD113_os;

%%
cfg = [];
cfg.layout = 'GSN-HydroCel-129.sfp';
lay = ft_prepare_layout(cfg);

%% ICA
cfg = [];
cfg.method ='fastica';
ic_data = ft_componentanalysis(cfg,data);
%%
% plot the components for visual inspection
figure,
cfg = [];
cfg.layout = lay; % specify the layout file for plotting
cfg.comment = 'no';
cfg.component = 1:20; % and more /  % specify the component(s) that should be plotted
ft_topoplotIC(cfg, ic_data);
%%
% remove the bad components and backproject the data
cfg = [];
cfg.component = [......]; % to remove component(s)
data = ft_rejectcomponent(cfg, ic_data, data);
%% Visual inspection
cfg=[];
%cfg.channel='all';
cfg.ylim = [-200 200];
ft_databrowser(cfg,data); % "data" is your imported data
%%
cfg = [];
cfg.method = 'summary';
data_clean = ft_rejectvisual(cfg, data);

%% channel repair and save
elec = ft_read_sens('GSN-HydroCel-129.sfp');
load  ....../label.mat; % use your own path
cfg = [];
cfg.method = 'spline';
cfg.missingchannel = setdiff(label,data_clean.label);
cfg.elec = elec;
data_cr = ft_channelrepair(cfg, data_clean);
save data_cr data_cr;



