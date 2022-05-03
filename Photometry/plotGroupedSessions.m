function [] = plotGroupedSessions(GroupedGCAMP)
%% Figure Properties
set(0,'defaultfigurecolor',[1 1 1])
FontSize = 12;
LineWidth = 2;
LineWidth_Legend = 3;
size_w = 150;
size_h = 300;
heat_size_w = 250;
heat_size_h = 125;
blue = [0, 0.4470, 0.7410];
red = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];
purple = [0.4940, 0.1840, 0.5560];
black = [0.25, 0.25, 0.25];
smooth_win = 20;
smooth_win_limited = 10;
interp_smooth_win = 1;
smooth_type = 'sgolay';
y_onset = [-0.5 2.0];
y_offset = [-0.5 2.0];
y_interp = [-0.5 2.0];
y_headentry = [-0.5 3.0];
x_onset = [-2 .1];
x_offset = [-.1 5];
%% Plot Single Session Data from GCAMP .mat
%% LP Onset
%% Mean LP ON Delta F Plot

y1 = smoothdata(GroupedGCAMP.LPON_Met,2,smooth_type,smooth_win);
y2 = smoothdata(GroupedGCAMP.LPON_Fail,2,smooth_type,smooth_win);
n1 = size(y1,1);
n2 = size(y2,1);

figure('Name',[GroupedGCAMP.training_day(1:4) ' Z Criteria LPON'],'NumberTitle','off','rend','painters','pos',[100 100 size_w size_h])
hold on
%Plot mean with shaded standard error
s = shadedErrorBar(GroupedGCAMP.plot_time, y1, {@mean, @(x) std(x) / sqrt(n1)}, 'lineprops', '-b', 'transparent',1);
set(s.edge,'LineWidth',6,'LineStyle','none')
set(s.mainLine,'Color',blue)
s.mainLine.LineWidth = LineWidth;
s.patch.FaceColor = blue;

t = shadedErrorBar(GroupedGCAMP.plot_time, y2, {@mean, @(x) std(x) / sqrt(n2)}, 'lineprops', '-r', 'transparent',1);
set(t.edge,'LineWidth',6,'LineStyle','none')
set(t.mainLine,'Color',red)
t.mainLine.LineWidth = LineWidth;
t.patch.FaceColor = red;

%this plots where the two different traces significantly differ at y = 0
plot(GroupedGCAMP.plot_time,GroupedGCAMP.Stats.LPON.p_val_for_graph, 'k:','LineWidth',3)
% x0=10;
% y0=10;
% width=800;
% height=600
% set(gcf,'position',[x0,y0,width,height])


% title(['Training Days: ' GroupedGCAMP.training_day(1:4) ', n = ' num2str(n1)])
title(['Onset'])
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'-', 'Color', blue);
h(2) = plot(NaN,NaN,'-', 'Color', red);
axis tight
legend(h,{'Rewarded', 'Unrewarded'})
set(h,'LineWidth',LineWidth_Legend);
legend boxoff
xlabel('Time from Onset (S)')
ylabel('Z-scored \Delta F / F')
set(gca,'FontSize',FontSize)
set(gca, 'FontName', 'Arial')
set(gca, 'TickDir', 'out')
%ylim([-0.5 2.0])
%ylim([-0.5 1.0])
ylim(y_onset)
%ylim([-0.5 1.0])
% ylim([-2 1.0]) %prob
xlim(x_onset)
%% Mean LP OFF Delta F Plot

y1 = smoothdata(GroupedGCAMP.LPOFF_Met,2,smooth_type,smooth_win);
y2 = smoothdata(GroupedGCAMP.LPOFF_Fail,2,smooth_type,smooth_win);
n1 = size(y1,1);
n2 = size(y2,1);

figure('Name',[GroupedGCAMP.training_day(1:4) ' Z Criteria LPOFF'],'NumberTitle','off','rend','painters','pos',[100 100 size_w size_h])
hold on
%Plot mean with shaded standard error
s = shadedErrorBar(GroupedGCAMP.plot_time, y1, {@mean, @(x) std(x) / sqrt(n1)}, 'lineprops', '-b', 'transparent',1);
set(s.edge,'LineWidth',6,'LineStyle','none')
set(s.mainLine,'Color',blue)
s.mainLine.LineWidth = LineWidth;
s.patch.FaceColor = blue;

t = shadedErrorBar(GroupedGCAMP.plot_time, y2, {@mean, @(x) std(x) / sqrt(n2)}, 'lineprops', '-r', 'transparent',1);
set(t.edge,'LineWidth',6,'LineStyle','none')
set(t.mainLine,'Color',red)
t.mainLine.LineWidth = LineWidth;
t.patch.FaceColor = red;

%this plots where the two different traces significantly differ at y = 0
plot(GroupedGCAMP.plot_time,GroupedGCAMP.Stats.LPOFF.p_val_for_graph, 'k:','LineWidth',3)

%title(['Training Days: ' GroupedGCAMP.training_day(1:4) ', n = ' num2str(n1)])
title(['Offset'])
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'-', 'Color', blue);
h(2) = plot(NaN,NaN,'-', 'Color', red);
axis tight
% legend(h,{'Rewarded', 'Unrewarded'})
% set(h,'LineWidth',LineWidth_Legend);
% legend boxoff
xlabel('Time from Offset (S)')
ylabel('Z-scored \Delta F / F')
set(gca,'FontSize',FontSize)
set(gca, 'FontName', 'Arial')
set(gca, 'TickDir', 'out')
ylim(y_offset)
% ylim([-2.0 2.0]) %prob
xlim(x_offset)

%% Mean First Head Entry After Reward vs no Reward Delta F Plot

y1 = smoothdata(GroupedGCAMP.FirstHE_RE,2,smooth_type,smooth_win);
y2 = smoothdata(GroupedGCAMP.FirstHE_Fail,2,smooth_type,smooth_win);
n1 = size(y1,1);
n2 = size(y2,1);

figure('Name',[GroupedGCAMP.training_day(1:4) ' Z Criteria FirstHE'],'NumberTitle','off','rend','painters','pos',[100 100 size_w size_h])
hold on
%Plot mean with shaded standard error
s = shadedErrorBar(GroupedGCAMP.plot_time, y1, {@mean, @(x) std(x) / sqrt(n1)}, 'lineprops', '-b', 'transparent',1);
set(s.edge,'LineWidth',6,'LineStyle','none')
set(s.mainLine,'Color',blue)
s.mainLine.LineWidth = LineWidth;
s.patch.FaceColor = blue;

t = shadedErrorBar(GroupedGCAMP.plot_time, y2, {@mean, @(x) std(x) / sqrt(n2)}, 'lineprops', '-r', 'transparent',1);
set(t.edge,'LineWidth',6,'LineStyle','none')
set(t.mainLine,'Color',red)
t.mainLine.LineWidth = LineWidth;
t.patch.FaceColor = red;


%this plots where the two different traces significantly differ at y = 0
plot(GroupedGCAMP.plot_time,GroupedGCAMP.Stats.HE_Fail.p_val_for_graph, 'k:','LineWidth',3)

%title(['Training Days: ' GroupedGCAMP.training_day(1:4) ', n = ' num2str(n1)])
title(['Head Entry'])
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'-', 'Color', blue);
h(2) = plot(NaN,NaN,'-', 'Color', red);
axis tight
% legend(h,{'Rewarded', 'Unrewarded'})
% set(h,'LineWidth',LineWidth_Legend);
% legend boxoff
xlabel('Time from First Head Entry (S)')
ylabel('Z-scored \Delta F / F')
set(gca,'FontSize',FontSize)
set(gca, 'FontName', 'Arial')
set(gca, 'TickDir', 'out')
ylim(y_headentry)

%% Interpolated Duration Activity

% prior to plotting, remove any NaNs (which are LPs withouth at least 2
% samples)

Grouped_interp_All = rmmissing(GroupedGCAMP.LPInterp);
Grouped_interp_Met = rmmissing(GroupedGCAMP.LPInterp_Met);
Grouped_interp_Fail = rmmissing(GroupedGCAMP.LPInterp_Fail);



y1 = smoothdata(Grouped_interp_Met,2,smooth_type,interp_smooth_win);
y2 = smoothdata(Grouped_interp_Fail,2,smooth_type,interp_smooth_win);
y3 = smoothdata(Grouped_interp_All,2,smooth_type,interp_smooth_win);
n1 = size(y1,1);
n2 = size(y2,1);
n3 = size(y3,1);

percent_of_press = (100/GroupedGCAMP.interp_length)*(1:GroupedGCAMP.interp_length);

figure('Name',[GroupedGCAMP.training_day(1:4) ' Z Criteria Duration'],'NumberTitle','off','rend','painters','pos',[100 100 size_w size_h])
hold on
%Plot mean with shaded standard error
s = shadedErrorBar(percent_of_press, y1, {@mean, @(x) std(x) / sqrt(n1)}, 'lineprops', '-b', 'transparent',1);
set(s.edge,'LineWidth',6,'LineStyle','none')
set(s.mainLine,'Color',blue)
s.mainLine.LineWidth = LineWidth;
s.patch.FaceColor = blue;

t = shadedErrorBar(percent_of_press, y2, {@mean, @(x) std(x) / sqrt(n2)}, 'lineprops', '-r', 'transparent',1);
set(t.edge,'LineWidth',6,'LineStyle','none')
set(t.mainLine,'Color',red)
t.mainLine.LineWidth = LineWidth;
t.patch.FaceColor = red;


%this plots where the two different traces significantly differ at y = 0
plot(percent_of_press,GroupedGCAMP.Stats.LPInterp.p_val_for_graph, 'k:','LineWidth',3)

% p = shadedErrorBar(percent_of_press, y3, {@mean, @(x) std(x) / sqrt(n3)}, 'lineprops', '-r', 'transparent',1);
% set(p.edge,'LineWidth',6,'LineStyle','none')
% set(p.mainLine,'Color',black)
% p.mainLine.LineWidth = LineWidth;
% p.patch.FaceColor = black;

%title(['Training Days: ' GroupedGCAMP.training_day(1:4) ', n = ' num2str(n1)])
title(['Duration'])
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'-', 'Color', blue);
h(2) = plot(NaN,NaN,'-', 'Color', red);
axis tight
% legend(h,{'Rewarded', 'Unrewarded'})
% set(h,'LineWidth',LineWidth_Legend);
% legend boxoff
xlabel('Lever Press Duration (%)')
ylabel('Z-scored \Delta F / F')
set(gca,'FontSize',FontSize)
set(gca, 'FontName', 'Arial')
set(gca, 'TickDir', 'out')
ylim(y_interp)
%ylim([-1.5 2.0]) %prob
end


