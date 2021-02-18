%% Analysis based on FieldTrip

%% Time-frequency analysis
load data; % use your own path /EEG data of each subject
cfg.method    = 'mtmconvol';
cfg.output    = 'pow';
cfg.taper     = 'hanning';
cfg.foi       = 8:1:12; % e.g. alpha
cfg.toi       = -0.75:0.005:0.25;  % or -0.25:0.005:0.75
cfg.t_ftimwin = 5./cfg.foi;  
cfg.trials = tr; % selected trials
data_fq   = ft_freqanalysis(cfg, data); % same method for 1st-stage, 3rd-stage, stimulus-locked and response-locked epochs
save data_fq data_fq;

%% Cluster statistics:
load ....../neighboursEEG128.mat;  % template from your own path
cfg = [];
cfg.channel          = 'all';
cfg.latency          =[-0.5 0];
cfg.frequency        = [8 12];  % e.g. alpha
cfg.avgoverfreq = 'yes';
cfg.method           = 'montecarlo';
cfg.numrandomization = 1000;  
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clustertail      = 0;
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.alpha            = 0.025;
cfg.neighbours       = neighbours;
 
subj = 22;
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,(subj+1):2*subj) = 2;
cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat] = ft_freqstatistics(cfg, GA1,  GA2); % GA1(Aha), GA2(Ctrl) from ft_freqgrandaverage; cfg.keepindividual = 'yes'

%% Topoplot for Aha, Ctrl, Diff
cfg = [];
cfg.xlim         = [-0.3 -0.2]; % select a time window
cfg.ylim         = [8 12]; % frequency
cfg.zlim         = [-3 3];  
cfg.layout       = 'GSN-HydroCel-129.sfp';
cfg.style        = 'straight';
cfg.markersymbo  =   '.';
cfg.comment      = 'no';
figure,
subplot(1,3,1)
ft_topoplotTFR(cfg, GA1avg);title('Aha'); % GA1avg: averaged grandaverage, from ft_freqgrandaverage; cfg.keepindividual = 'no'
subplot(1,3,2)
ft_topoplotTFR(cfg, GA2avg);title('Ctrl');
subplot(1,3,3)
ft_topoplotTFR(cfg, GAdifavg);title('Diff');
%% Topoplot for Diff cluster
cfg = [];
cfg.parameter = 'stat';
cfg.zlim = [-3 3];
cfg.layout = 'GSN-HydroCel-129.sfp';
cfg.highlightsymbolseries = ['.','.','.','.','.'];
cfg.highlightcolorneg = 'w';
cfg.style = 'straight';
cfg.marker = 'on';
cfg.markersymbol = '.';
cfg.highlightsizeseries = [15 15 15 15 15];
ft_clusterplot(cfg, stat);
 