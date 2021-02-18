%  Behavior
 
%% Collect correct-responsed trials
% filename={......}; % subjects
RT1=cell(1,25);
RT2=cell(1,25);
lost=NaN(1,25);
for sub=1:length(filename)
    Subject=filename{sub};
    loaddata = ['load ......, Subject, '_piclistreact_new' '.mat'];    eval(loaddata);
    tr1= piclistreact(:,2)==1 & piclistreact(:,8)==1;
    RT1{1,sub}=piclistreact(tr1,6)/10000;
    tr2= piclistreact(:,2)==2 & piclistreact(:,8)==1;
    RT2{1,sub}=piclistreact(tr2,6)/10000;
    lost(1,sub)=160-length(piclistreact);
end
%% Collect non-correct-response trials
% filename={......}; % subjects
RT1cuo=cell(1,25);
RT2cuo=cell(1,25);
lost=NaN(1,25);
for sub=1:length(filename)
    Subject=filename{sub};
    loaddata = ['load ......, Subject, '_piclistreact_new' '.mat'];    eval(loaddata);
    tr1= piclistreact(:,2)==1 & piclistreact(:,8)==0;
    RT1cuo{1,sub}=piclistreact(tr1,6)/10000;
    tr2= piclistreact(:,2)==2 & piclistreact(:,8)==0;
    RT2cuo{1,sub}=piclistreact(tr2,6)/10000;
end
%% Count
a1=cellfun(@length, RT1);
a0=cellfun(@length, RT1cuo);
c1=cellfun(@length, RT2);
c0=cellfun(@length, RT2cuo); 
%
All1=[];
All2=[];
for i=1:25
    All1=[All1;RT1{1,i}];
    All2=[All2;RT2{1,i}];
end
%
c1=0;c2=0;
for i=1:25
    c1=c1+length(RT1{1,i});
    c2=c2+length(RT2{1,i});
end
%
med1=NaN(1,25);
med2=NaN(1,25);
for i=1:25
    med1(1,i)=median(RT1{1,i});
    med2(1,i)=median(RT2{1,i});
end

%% Plot1
figure,
b1=bar(x,y1);  % Percentage for Aha
set(b1,'FaceColor','r','EdgeColor','w');
alpha(0.4);
hold on
b2=bar(x,y2);  % Percentage for Ctrl
set(b2,'FaceColor','b','EdgeColor','w')
alpha(0.4);
hold on
c1 = polyfit(x, y1, 7);
d1 = polyval(c1, x, 1);
c2 = polyfit(x, y2, 7);
d2 = polyval(c2, x, 1);
curve=plot(x,d1,'r',x,d2,'b','LineWidth',2); hold on

xlabel('Reaction Time/s','fontsize',12);
ylabel('Distribution/%','fontsize',12);
axis xy
h=legend(curve,'Aha','Ctrl',1,'FontName','Helvetica');
xlim([0 8]);
ylim([0 5]);
hold off
 
%% Plot2
figure,
boxplot([med1', med2'],'orientation','horizontal','color','rb'); 
axis ij
xlim([0 8]);
xlabel('Reaction Time/s','fontsize',12);

%% Normality test
lg=25;
PPn=NaN(1,lg);
HHn=NaN(1,lg);
for j=1:lg % subs
    [h,p] = lillietest(RT1{1,j},0.05); % same for RT2
    PPn(1,j)=p;
    HHn(1,j)=h;
end


