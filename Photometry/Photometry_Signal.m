function [GCAMP] = Photometry_Signal(GCAMP, session_file_name_photometry)
[photometry_data] = photometry_extract(session_file_name_photometry);

%% Photometry Data Extraction and Normalization
raw_gcamptimestamps = photometry_data(:,1);
raw_gcampdata = photometry_data(:,2);
% if there are an odd # of rows after removing headers, delete the last row
% in both wave and wavetime
A = raw_gcampdata(1:2:end,:);
B = raw_gcampdata(2:2:end,:);

%truncate the data from the LED that has more datapoints to match data
%dimension
if length(A) > length(B)
  A = A(1:end-1,:);
elseif length(A) < length(B)
  B = B(1:end-1,:);
end

%the raw fluorescence values from LED that are consistently greater belong
%to 470 nm (i.e. LED1)
if sum(A>B) > (length(A)*.99) %% 0.99 seems arbitrary here?
  LED1 = A;
  LED2 = B;
else
  LED1 = B;
  LED2 = A;
end

% Timestamps for gcampdata
if LED1 == A
    k = 1;
elseif LED1 == B
    k = 2;
end
gcampdata_timestamps = raw_gcamptimestamps(k:2:end,:);

%% Denoising the signal 
LED1dn = medfilt1(LED1,5); %median filter gets rid of random/ salt and pepper noise 
%If you have rare outliers an order of 5 should be good. 
%Otherwise you can increase to 7, 9 or even 11 but not more because you`ll
%distort.
fc = 10; %cutoff frequency in Hz (anything lower than this passes)
fs = 20; %sampling rate in Hz
%below returns transfer function coefficients of an nth-order lowpass digital filter 
%Butterworth filter with normalized cutoff frequency Wn
%For digital filters, the cutoff frequencies must lie between 0 and 1
%where 1 corresponds to the Nyquist frequency
%rad/sample)
[b,a] = butter(2,0.99999999,'low');
LED1dn2 = filtfilt(b,a,LED1dn); %filtfilt is zero-phase so no temporal distortion

%% High-pass filter method of accounting for slow decay 
fc2 = 0.001; %cutoff frequency in Hz (anything higher than this f passes)
[b2,a2] = butter(2,fc2/(fs/2),'high'); %this accounts for slow time course changes on the order of ~16 minutes
LED1c = filtfilt(b2,a2,LED1dn2); %filtfilt is zero-phase
gcampdata = LED1c+50000; %add a constant so the trace doesn't go down to zero on y-axis
SR = 20;
time = (1:length(gcampdata))/SR; % in seconds
%% Obtain Baseline Trace
%estimate 10th percentile of F using 15 sec. sliding window 
%this is to generate baseline fluorescence estimate
win = SR * 15;
winright = 0;
p = 10;
%[ ys ] = running_percentile(gcampdata, win, p); %ege
[ ys ] = running_percentile_filter(gcampdata, win, winright, p); %drew
delta_F = gcampdata - ys;
delta_FoverF = delta_F ./ ys ;
%exclusion criterion: only sessions where 97.5th percentile of delta_FoverF
%exceeded 1% are included
delta_F_p = delta_FoverF .* 100; %percent change from baseline fluorescence
exc = prctile(delta_F_p, 97.5)

pass = (exc > 1) 
if pass == 0
    ['excluding ' session_file_name_photometry]
end

%% Plot raw GCAMP data
% Make some pretty colors for later plotting
% http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
red = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];
cyan = [0.3010, 0.7450, 0.9330];
yellow = [0.9290, 0.6940, 0.1250];
purple = [0.4940, 0.1840, 0.5560];
gray1 = [.7 .7 .7];
gray2 = [.8 .8 .8];
figure('Name',[session_file_name_photometry],'NumberTitle','off','Position',[100, 100, 800, 400])
hold on;
p1 = plot(time, LED1,'color',red,'LineWidth',1);
p2 = plot(time, gcampdata,'color',cyan,'LineWidth',1);
p3 = plot(time, ys,'color',green,'LineWidth',1);
%p3 = plot(time, Y_exp_fit_all(time),'color',green,'LineWidth',10);
p4 = plot(time,delta_F_p,'color',gray1,'LineWidth',1);
ylabel('F','fontsize',16);
title([session_file_name_photometry])
axis tight;
legend([p1 p2 p3 p4], {'Raw Signal', 'Denoised and Decay Corrected','Baseline','df/F'});

figure
p4 = plot(time,delta_F_p,'color',gray1,'LineWidth',1);
%% Save data to GCAMP structure

GCAMP.gcampdata = delta_F_p;
GCAMP.gcampdata_timestamps = gcampdata_timestamps;
GCAMP.gcampdata_delta_F = delta_F;
GCAMP.gcampdata_delta_FoverF = delta_FoverF;
GCAMP.LED1 = LED1;
GCAMP.LED2 = LED2;
GCAMP.raw_gcamptimestamps = raw_gcamptimestamps;
GCAMP.raw_gcampdata = raw_gcampdata;
GCAMP.SR = SR;

end

