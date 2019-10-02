%% Plots the model output to generate figure 2
load('DataNorthKivu.mat'); % load the data the model was fit to
close all;
% Specify dimensions for the subplots
hei=0.36;
wid=0.42;
xxs=linspace(0.075,0.99-0.42,2);
yys=linspace(0.97-hei,0.145,2);

figure('units','normalized','outerposition',[0 0 1 1])

f=3; % Saturation function used 
load(['ModelFit']); % Load data to plot    

    startDateofSim = datenum('04-30-2018'); % Start date of the outbreak is April 30, 2018
    XTL=datestr([startDateofSim+7.*[0:5:59]],'mm/dd/yy'); % Index for x-tick in the figure
    endDate = startDateofSim+15*7; %% Doing the week vacciantion started and not the day datenum('08-12-2018');
    TVac=endDate-startDateofSim;   % Week of vaccinatioon
%% Incidence 
subplot('Position',[xxs(1) yys(1) wid hei]);

for ii=1:length(IData)+1
    p3=patch(ii-1+[-0.35 0.35 0.35 -0.35],[UBI(ii) UBI(ii) LBI(ii) LBI(ii)],[0 0 1],'FaceAlpha',.2); hold on % Plot upper and lower bounds
    set(p3,'LineStyle','none');     
    plot(ii-1+linspace(-0.35,0.35,101),mleI(ii).*ones(101,1),'color',[0 0 1],'LineWidth',2); % plot the mle
end
scatter([1:length(IData)],IData,40,'k','filled'); % plot the data
plot((TVac./7).*ones(1001,1),linspace(0,140,1001),'k-.','LineWidth',1.5); hold on; % plot the line for vaccinatoon
xlim([0 length(IData)+0.35]); % Set xlim
ylim([0 140]); % set ylim
box off
set(gca,'linewidth',2,'tickdir','out','Fontsize',16,'YTick',[0:20:140],'XTick',[1:5:length(IData)],'XMinorTick','on','XTickLabel',XTL,'YMinortick','on')
h=ylabel({'Number of','cases'},'Fontsize',20);
xtickangle(45)
posll=h.Extent(1);

 %% Avg Reff over the week
 subplot('Position',[xxs(2) yys(1) wid hei]);

for ii=1:length(IData)
    p3=patch(ii+[-0.35 0.35 0.35 -0.35],[UBRW(ii) UBRW(ii) LBRW(ii) LBRW(ii)],[0 0 1],'FaceAlpha',.2); hold on % plot upper/lower bound for R_E
    set(p3,'LineStyle','none');     
 
    p3=patch(ii+[-0.35 0.35 0.35 -0.35],[IQRRW(1,ii) IQRRW(1,ii) IQRRW(2,ii) IQRRW(2,ii)],[0 0 1],'FaceAlpha',.25); hold on % plot IQR for R_
    set(p3,'LineStyle','none');     
    plot(ii+linspace(-0.35,0.35,101),mleRW(ii).*ones(101,1),'color',[0 0 1],'LineWidth',2); % plot mle
end

plot(linspace(0,60,1001),ones(1,1001),'-.','color',[0.6 0.6 0.6],'LineWidth',1.5); %plot line for R_E=1

plot((TVac./7).*ones(1001,1),linspace(0,3,1001),'k-.','LineWidth',1.5); hold on; % plot line ofr vaccination
xlim([0 length(IData)+.35]);
ylim([0 3]);
box off
set(gca,'linewidth',2,'tickdir','out','Fontsize',16,'YTick',[0:0.5:5],'XTick',[1:5:length(IData)],'XMinorTick','on','XTickLabel',XTL,'YMinortick','on')
%xlabel('Weeks','Fontsize',20);
h=ylabel({'Reproductive number','({\it R}_E )'},'Fontsize',20);
xtickangle(45)
posrl=h.Extent(1);

%% Effectiveness of vaccination
 
subplot('Position',[xxs(2) yys(2) wid hei]);
for ii=16:length(IData)
    p3=patch(ii+[-0.35 0.35 0.35 -0.35],[UBEW(ii) UBEW(ii) LBEW(ii) LBEW(ii)],[0 0 1],'FaceAlpha',.2); hold on % plot upper/lower bound
    set(p3,'LineStyle','none');     

    p3=patch(ii+[-0.35 0.35 0.35 -0.35],[IQREW(1,ii) IQREW(1,ii) IQREW(2,ii) IQREW(2,ii)],[0 0 1],'FaceAlpha',.25); hold on  % plot IQR
    set(p3,'LineStyle','none');     
    plot(ii+linspace(-0.35,0.35,101),mleEW(ii).*ones(101,1),'color',[0 0 1],'LineWidth',2); % plot mle
end

plot((TVac./7).*ones(1001,1),linspace(0,1,1001),'k-.','LineWidth',1.5); hold on; % plot time of vaccination
xlim([0 length(IData)+.35]);
ylim([0 1]);
box off
set(gca,'linewidth',2,'tickdir','out','Fontsize',16,'YTick',[0:0.2:1],'XTick',[1:5:length(IData)],'XMinorTick','on','XTickLabel',XTL,'YMinortick','on')
xlabel('Week of illness onset','Fontsize',20);
h=ylabel({'Population level','effectiveness of vaccination'},'Fontsize',20);
xtickangle(45)
text(min(posrl,h.Extent(1)),max(ylim),char(64+4),'Fontsize',32,'FontWeight','bold');

subplot('Position',[xxs(2) yys(1) wid hei]);
text(min(posrl,h.Extent(1)),max(ylim),char(64+2),'Fontsize',32,'FontWeight','bold');                
          
%% Time to isolation
         
 subplot('Position',[xxs(1) yys(2) wid hei]);
        
for ii=1:length(IData)
    p3=patch(ii+[-0.35 0.35 0.35 -0.35],[UBGW(ii) UBGW(ii) LBGW(ii) LBGW(ii)],[0 0 1],'FaceAlpha',.2); hold on % plot upper/lower bound
    set(p3,'LineStyle','none');    

    p3=patch(ii+[-0.35 0.35 0.35 -0.35],[IQRGW(1,ii) IQRGW(1,ii) IQRGW(2,ii) IQRGW(2,ii)],[0 0 1],'FaceAlpha',.25); hold on % plot IQR
    set(p3,'LineStyle','none');  
    plot(ii+linspace(-0.35,0.35,101),mleGW(ii).*ones(101,1),'color',[0 0 1],'LineWidth',2);% plot mle
end

plot((TVac./7).*ones(1001,1),linspace(0,100,1001),'k-.','LineWidth',1.5); hold on; % plot time of vaccination
xlim([0 length(IData)+.35]);
ylim([0 12]);
box off
set(gca,'linewidth',2,'tickdir','out','Fontsize',16,'YTick',[0:2:16],'XTick',[1:5:length(IData)],'XMinorTick','on','XTickLabel',XTL,'YMinortick','on')
xlabel('Week of illness onset','Fontsize',20);
h=ylabel({'Time to isolation','(days)'},'Fontsize',20);
xtickangle(45)
text(min(posll,h.Extent(1)),max(ylim),char(64+3),'Fontsize',32,'FontWeight','bold');
