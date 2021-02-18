%% Dimensionality

%% 10 channel groups (clusters)
group =cell(1,10);
group{1} ={'E2','E3','E4','E5','E9','E10','E123','E124'}'; % right-Frontal --- 8 channels
group{2} ={'E115','E114','E113','E109','E108','E107','E102','E101','E100'}'; % right-Tempral --- 9 channels
group{3} ={'E103','E104','E105','E106','E110','E111','E112','E116','E117','E118'}'; % right-Central --- 10 channels
group{4} ={'E80','E87','E93','E98','E97','E79','E86','E92','E78'}'; % right-Parietal --- 9 channels
group{5} ={'E77','E85','E91','E96','E76','E84','E90','E95','E83'}'; % right-Occipital --- 9 channels
group{6} ={'E22','E18','E26','E23','E19','E12','E27','E24'}'; % left-Frontal --- 8 channels
group{7} ={'E49','E44','E39','E56','E45','E40','E57','E50','E46'}'; % left-Tempral --- 9 channels
group{8} ={'E34','E28','E20','E35','E29','E13','E41','E36','E30','E7'}'; % left-Central --- 10 channels
group{9} ={'E47','E42','E37','E31','E51','E52','E53','E54','E61'}'; %  left-Parietal --- 9 channels
group{10} ={'E58','E59','E60','E67','E64','E65','E66','E71','E70'}'; % left-Occipital --- 9 channels

% load data  % here e.g. dataipcD113_os_cl (Aha, stimulus-locked)
dmg1 = NaN(10,1001);
for gp=1:10
    chan=length(group{gp});
    [~,~,zz]= intersect(group{gp}, dataipcD113_os_cl.label);
    xi=0:0.1:100;
    % tr = trials for Aha or Ctrl
    dm=zeros(length(tr),length(xi));
    for j=1:length(tr)
        test=dataipcD113_os_cl.trial{tr(j)}(zz,501:1001); % time: -0.5 - 0 sec or 0 - 0.5 sec
        tm=mean(test,2);    ts=std(test,0,2);
        tm2=repmat(tm,[1 501]);   ts2=repmat(ts,[1 501]);
        tnew=(test-tm2)./ts2;
        r=tnew';
        dv=zeros(size(r,1),size(r,1),chan);
        for i = 1:chan
            dv(:,:,i) = bsxfun(@minus,r(:,i),r(:,i)');
        end
        temp=sum(dv.^2,3);
        dst = sqrt(temp);
        msk = triu(ones(size(dv,1)),1);
        msk(msk==0) = NaN;
        numbs = msk.*dst;
        numbs = numbs(:);
        mnmx = [min(numbs),max(numbs)];
        rng = linspace(mnmx(1),mnmx(2),1000);
        c = 1;
        for ri = rng
            s(c) = sum(numbs<ri); % count the number < ri
            c = c+1;
        end
        gd=gradient(log(s),log(rng));
        dm(j,:)=interp1(rng,gd,xi,'linear'); % size: length(tr)*length(xi)
    end
    dm1=nanmean(dm,1); % size: 1*length(xi)
    dmg1(gp,:) = dm1;
end
save dmg1 dmg1  % others same

%% calculate LI for Aha, Ctrl, 1st stage, 3rd stage
LI1=(dmg1(1:5,:)-dmg1(6:10,:))./(dmg1(1:5,:)+dmg1(6:10,:));
LI2=(dmg2(1:5,:)-dmg2(6:10,:))./(dmg2(1:5,:)+dmg2(6:10,:));
LIbl1=(dmblg1(1:5,:)-dmblg1(6:10,:))./(dmblg1(1:5,:)+dmblg1(6:10,:));
LIbl2=(dmblg2(1:5,:)-dmblg2(6:10,:))./(dmblg2(1:5,:)+dmblg2(6:10,:));
LIgroup1=LI1-LIbl1;
LIgroup2=LI2-LIbl2;
% save data

%% stat: signrank test for LI
PP=zeros(5,1001);
HH=zeros(5,1001);
for i=1:5
    for j=1:1001
        t1=squeeze(LIgroupall1(:,i,j)); % for all subjects
        t2=squeeze(LIgroupall2(:,i,j));
        t3=t1-t2;
        t3=t3(~isnan(t3));
        if isempty(t3) || length(t3)<2
            P=NaN;
            H=0;
        else
            [P,H] = signrank(t1,t2);
        end
        PP(i,j)=P;
        HH(i,j)=H;
    end
end
%% plot LI
gps={'F','T','C','P','O'}';
xi=0:0.1:100;
figure;
set(gcf,'outerposition',get(0,'screensize'));
set(0, 'defaultTextInterpreter', 'tex');
for gp=1:5
    subplot(3,5,gp)
    for j=1:1001
        if PP(gp,j)<0.05
            y=[-10 10];
            x=[xi(j) xi(j)];
            bg=plot(log(x),y);
            set(bg,'Color',[0.7 0.7 0.7],'LineWidth', 1.5);hold on;
        end
    end
    hold on
    plot(log(xi),LIgroupavg1(gp,:),'r',log(xi),LIgroupavg2(gp,:),'b','LineWidth',2);
    xlim([0 2]);
    ylim([-.05  .05]);
    box on
    set(gca,'layer','top');
    title(char(gps{gp}));
end

%% stat: signrank test for dim
PP=zeros(10,1001);
HH=zeros(10,1001);
for i=1:10
    for j=1:1001
        t1=squeeze(dmgroupall1(:,i,j));
        t2=squeeze(dmgroupall2(:,i,j));
        t3=t1-t2;
        t3=t3(~isnan(t3));
        if isempty(t3) || length(t3)<2
            P=NaN;
            H=0;
        else
            [P,H] = signrank(t1,t2);
        end
        PP(i,j)=P;
        HH(i,j)=H;
    end
end
%% curve plot
gps={'RF','RT','RC','RP','RO','LF','LT','LC','LP','LO'}';
xi=0:0.1:100;
xii=log(xi);
figure;
set(gcf,'outerposition',get(0,'screensize'));
set(0, 'defaultTextInterpreter', 'tex');
for gp=1:10
    subplot(3,5,gp);
    for j=1:1001
        if PP(gp,j)<0.05
            y=[-10 10];
            x=[xi(j) xi(j)];
            bg=plot(log(x),y);
            set(bg,'Color',[0.7 0.7 0.7],'LineWidth',1.5);hold on;
        end
    end
    hold on
    curve=plot(xii,dmgroupavg1(gp,:),'r',xii,dmgroupavg2(gp,:),'b');
    set(curve, 'LineWidth',1.5);
    xlim([0 2]);
    ylim([-.1 .1]);
    box on
    set(gca,'layer','top');
    legend(curve,'Aha','Ctrl',1);
    title(char(gps{gp}));
    hold off
end

%% Note
% Cluster statistics same as that in Coh_PLV analysis.
