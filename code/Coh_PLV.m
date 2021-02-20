%% Analysis based on FieldTrip

%% Calculate coherence/PLV
% Coherence/PLV analysis shares the same procedure.
% Only different in ft_connectivityanalysis: cfg.method = 'coh' or 'plv'.

load label; % 128 channels, from template

% Calculate freq first:
load data; % use your own path /EEG data of each subject
cfg.method    = 'mtmconvol';
cfg.output    = 'fourier';
cfg.taper     = 'hanning';
cfg.foi       = 4:1:12;  % e.g. theta and alpha
cfg.toi       = -0.75:0.01:0.25;  % or -0.25:0.01:0.75
cfg.t_ftimwin = 5./cfg.foi;
cfg.keeptrials = 'yes'; %
cfg.trials = tr; % selected trials
data_fourier  = ft_freqanalysis(cfg, data);
save data_fourier data_fourier;

% Calculate coherence/PLV:
cfg = [];
cfg.method = 'coh'; % or cfg.method = 'plv', others same as coherence analysis;
data_coh = ft_connectivityanalysis(cfg, data_fourier);

%% An example of changing the channel order (to ensure that data arrayed from E1 to E128):
chanorder=zeros(128,1);
for i=1:128
    chanorder(i) = find(ismember(D111_OS_L_fourier.label,channel(i)));
end
z = D111_OS_L_coh_1.cohspctrm;
z2=zeros(128,128,9,101); % new coh matrix
OScoh128=zeros(128,128,9,101);
for i=1:128
    z2(i,:,:,:)=z(chanorder(i),:,:,:); % correct dimension 1
end
for i=1:128
    OScoh128(:,i,:,:)=z2(:,chanorder(i),:,:); % correct dimmension 2;
end
D111_1_OScoh128=OScoh128;
save D111_1_OScoh128 D111_1_OScoh128; % e.g. coherence data of 1st stage, Aha, and stimulus-locked

%% Grand analysis
% filename={......}; % all subjects
allD111_1_OScoh128 = zeros(22,128,128,9,101); % e.g. coherence data of 1st stage, Aha, and stimulus-locked; Others same
for sub=1:length(filename)
    Subject=filename{sub};
    loaddata = ['load  ....../D111_1_OScoh128' '.mat']; eval(loaddata);
    allD111_1_OScoh128(sub,:,:,:,:) = D111_1_OScoh128;
end

x=allD113_1_OScoh128;%
y=allD111_1_OScoh128(:,:,:,:,26:76);%
y2=nanmean(y,5); %
y3=repmat(y2,[1 1 1 1 101]); %
chanOScoh128_1_cor = x-y3; %
chanOScoh128_1_cor_mean = squeeze(nanmean(chanOScoh128_1_cor,1));

%% Coh for topoplot and matrix plot
% get the average of alpha coh (-0.5~0 sec or 0~0.5 sec):
a_aha=nanmean(nanmean(chanOScoh128_1_cor_mean(:,:,5:9,26:76),3),4);   %alpha: 5:9
a_ctrl=nanmean(nanmean(chanOScoh128_2_cor_mean(:,:,5:9,26:76),3),4);
a_diff=nanmean(nanmean(diff(:,:,5:9,26:76),3),4);

%% Topoplot
cfg = [];
cfg.layout = 'GSN-HydroCel-129.sfp';
lay = ft_prepare_layout(cfg);
nochan = [8,14,17,21,25,125:128]; % not plotted channels
chan = 128;
tr = 0.05; % transparent value
lw = 1; % linewidth
colr = jet(100); % colormap
high = 0.09; % data value area % adjust according to the data
low = -high;

figure,
for i=1:length(lay.outline)
    if ~isempty(lay.outline{i})
        X = lay.outline{i}(:,1);
        Y = lay.outline{i}(:,2);
        h = line(X, Y);
        set(h, 'color', 'k');
        set(h, 'linewidth', 2);
    end
end
hold on
scatter(lay.pos(4:131,1),lay.pos(4:131,2),1,'MarkerEdgeColor',[1 1 1]);
axis auto
axis equal
axis off
hold on
data= a_aha; % e.g. here is alpha coh/or PLV for Aha, same method for theta, and for Ctrl and Diff
for m = 1:chan
    for n = 1:chan
        if m>n
            if ~ismember(m, nochan) && ~ismember(n, nochan)
                x=[lay.pos(m+3,1) lay.pos(n+3,1)]; % +3 means start with E1 (at raw 4th, from 4 to 131 (E1 - E128))
                y=[lay.pos(m+3,2) lay.pos(n+3,2)];
                if data(m,n)>high
                    data(m,n)=high-0.001; % to avoid extremum that exceed the bound
                elseif data(m,n)<low
                    data(m,n)=low+0.001;
                end
                patchline(x, y,'edgecolor', colr(ceil(100*(data(m,n)-low)/(high-low)),:),'linewidth',lw,'edgealpha',tr);  % tr means transparent value ******************
                % patchline is a new function, I put it in "function" file
                % patchline(https://www.mathworks.com/matlabcentral/fileexchange/36953-patchline), by Brett Shoelson
                hold on
            else
            end
        else
        end
    end
end
colormap(colr);
colorbar('FontWeight','Bold','FontSize',12);
title('Aha','FontName','Helvetica','FontWeight','Bold','fontsize',12);
caxis([low, high]);
hold off

%% LI analysis
% To obtain the following variables:
% D111_1_L, D111_1_R
% D111_2_L, D111_2_R
% D113_1_L, D113_1_R
% D113_2_L, D113_2_R
%
% L: left; R: right
% Here the variables are from coh/PLV that from Fp1, Fp2, F7, F3, Fz, F4, F8,
% T3, C3, CPz, C4, T4, T5, P3, Pz, P4, T6, O1, Oz, and O2.

freqbin=5:9;
for i=1:4
    if i==1
        mtxin = allD111_1_OScohnew;
    elseif i==2
        mtxin = allD111_2_OScohnew;
    elseif i==3
        mtxin = allD113_1_OScohnew;
    elseif i==4
        mtxin = allD113_2_OScohnew;
    else
    end
    
    % left side pairs:
    mtxout_L = zeros(60,22,9,201); %  60 pairs; 22 subjects
    cp =0;
    % Part 1:
    for m=[1 3 4 8 9 13 14 18]
        for n=[1 3 4 8 9 13 14 18]
            if m>n
                cp=cp+1;
                mtxout_L(cp,:,:,:)=squeeze(mtxin(:,m,n,:,:));
            end
        end
    end
    % Part 2:
    for m=[1 3 4 8 9 13 14 18]
        for n=[5 10 15 19]
            cp=cp+1;
            mtxout_L(cp,:,:,:)=squeeze(mtxin(:,m,n,:,:));
        end
    end
    
    % right side pairs:
    mtxout_R = zeros(60,22,27,201);
    cp =0;
    % Part 1:
    for m=[2 6 7 11 12 16 17 20]
        for n=[2 6 7 11 12 16 17 20]
            if m>n
                cp=cp+1;
                mtxout_R(cp,:,:,:)=squeeze(mtxin(:,m,n,:,:));
            end
        end
    end
    % Part 2:
    for m=[2 6 7 11 12 16 17 20]
        for n=[5 10 15 19]
            cp=cp+1;
            mtxout_R(cp,:,:,:)=squeeze(mtxin(:,m,n,:,:));
        end
    end
    
    if i==1
        D111_1_L =  squeeze(nanmean(nanmean(mtxout_L(:,:,freqbin,:),3),1));
        D111_1_R =  squeeze(nanmean(nanmean(mtxout_R(:,:,freqbin,:),3),1));
    elseif i==2
        D111_2_L =  squeeze(nanmean(nanmean(mtxout_L(:,:,freqbin,:),3),1));
        D111_2_R =  squeeze(nanmean(nanmean(mtxout_R(:,:,freqbin,:),3),1));
    elseif i==3
        D113_1_L =  squeeze(nanmean(nanmean(mtxout_L(:,:,freqbin,:),3),1));
        D113_1_R =  squeeze(nanmean(nanmean(mtxout_R(:,:,freqbin,:),3),1));
    elseif i==4
        D113_2_L =  squeeze(nanmean(nanmean(mtxout_L(:,:,freqbin,:),3),1));
        D113_2_R =  squeeze(nanmean(nanmean(mtxout_R(:,:,freqbin,:),3),1));
    else
    end
end
%% LI
LI1131 = (D113_1_R - D113_1_L)./(D113_1_R + D113_1_L); % Aha
LI1132 = (D113_2_R - D113_2_L)./(D113_2_R + D113_2_L); % Ctrl
LI1111 = (D111_1_R - D111_1_L)./(D111_1_R + D111_1_L); % Aha
LI1112 = (D111_2_R - D111_2_L)./(D111_2_R + D111_2_L); % Ctrl
% remove baseline for Aha
x=LI1131;
y=LI1111(:,51:151);
y2=nanmean(y,2);
y3=repmat(y2,[1 201]);
LI_1_cor = x-y3;
% same method for Ctrl

%% Wilcoxon signrank test
lg=201;
PPP=NaN(1,lg);
HHH=NaN(1,lg);
STATz=NaN(1,lg);
for j=1:lg
    t1=squeeze(LI_1_cor(:,j));
    t2=squeeze(LI_2_cor(:,j));
    t3=t1-t2;
    t3=t3(~isnan(t3));
    if isempty(t3) || length(t3)<2
        P=NaN;
        H=NaN;
        STATz(1,j)=NaN;
    else
        [P,H,stats] = signrank(t1,t2);
        STATz(1,j)=stats.zval;
    end
    PPP(1,j)=P;
    HHH(1,j)=H;
end

%% Cluster statistics
repeat = 1000;
total = zeros(1,repeat);
samp=22;
bins=65:95;  % the bins for test, e.g.,65:95, depend on the bins of interest
lb=length(bins);
sample1=data1(:,bins);
sample2=data2(:,bins);
for k=1:repeat
    new_1 = NaN(samp,lb);
    new_2 = NaN(samp,lb);
    
    for sub=1:samp
        pick = randi([1,2],1);
        if pick==1
            new_1(sub,:) = sample1(sub,:);
            new_2(sub,:) = sample2(sub,:);
        elseif pick==2
            new_1(sub,:) = sample2(sub,:);
            new_2(sub,:) = sample1(sub,:);
        else
        end
    end
    % stat:
    STATv=NaN(1,lb); % STATv: stat value
    ct=0;  % ct: count the lg;
    for j=1:lb
        ct=ct+1;
        t1=squeeze(new_1(:,j));
        t2=squeeze(new_2(:,j));
        t3=t1-t2;
        t3=t3(~isnan(t3));
        if isempty(t3) || length(t3)<2
            STATv(1,ct)=NaN;
        else
            [P,H,stats] = signrank(t1,t2);
            if  isfield(stats, 'zval')
                STATv(1,ct)=stats.zval;
            else
                STATv(1,ct)=NaN;
            end
        end
    end
    total(1,k) = nansum(STATv);
end
total95 = prctile(total,95);
% To compare:
tsum = nansum(abs(STATz(1,bins)));
