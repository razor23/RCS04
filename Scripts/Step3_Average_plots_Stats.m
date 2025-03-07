%% Analyse Dystonia data %%

clear all; close all; clc

% Load up the data folder with all the files in %

InDir = '/Volumes/Vinith2/Starr/Dystonia/RawData';
Groups = {'RCS04L','RCS04R'};
Devices = {'DeviceNPC700418H','DeviceNPC700412H'};

%% Load file

for G=1%:length(Groups)
    load(sprintf('%s/%s/StimSess_gamma_Dec5.mat',InDir,Groups{G}));
    ONfiles  = StimSess(:,1);
    OFFfiles = StimSess(:,2);
    ONfiles =  ONfiles(~cellfun('isempty',ONfiles)); % there are some empty cells, geting rid of em
    OFFfiles = OFFfiles(~cellfun('isempty',OFFfiles)); % there are some empty cells, geting rid of em

   
    
    coherence1ON{G} = [];   coherence1OFF{G} = [];
    coherence2ON{G} = [];   coherence2OFF{G} = [];
    coherence3ON{G} = [];   coherence3OFF{G} = [];

    psdout1ON{G} = [];  psdout1OFF{G} = [];
    psdout2ON{G} = [];  psdout2OFF{G} = [];
    psdout3ON{G} = [];  psdout3OFF{G} = [];
    psdout4ON{G} = [];  psdout4OFF{G} = [];
    % ON files data
 
fprintf('\n Printing ON files \n');

    
    for f=[1:3 9 10 15]
        
        try
        FlName = ONfiles{f};
        fprintf('\n Loading file %s from %s',FlName,Groups{G});
        path = sprintf('%s/%s/%s/%s',InDir,Groups{G},FlName,Devices{G});
        Infile = sprintf('%s/%s/%s/%s/processed_data_gamma_Dec5.mat',InDir,Groups{G},FlName,Devices{G});
        load(Infile);

         psdout1ON{G} = cat(2, psdout1ON{G}(:,:), Data.psd0);
         psdout3ON{G} = cat(2, psdout3ON{G}(:,:), Data.psd2);
         psdout4ON{G} = cat(2, psdout4ON{G}(:,:), Data.psd3);
         
         coherence1ON{G} = cat(2,coherence1ON{G}(:,:),Data.cxy);
         coherence2ON{G} = cat(2,coherence2ON{G}(:,:),Data.cyz);
         coherence3ON{G} = cat(2,coherence3ON{G}(:,:),Data.czx);
         
        catch
            fprintf('\n error \n');
            continue
        end
    end


  fprintf('\n Printing OFF files \n');

    
    for f=1:length(OFFfiles)
         try
        FlName = OFFfiles{f};
        fprintf('\n Loading file %s from %s',FlName,Groups{G});
        path = sprintf('%s/%s/%s/%s',InDir,Groups{G},FlName,Devices{G});
        Infile = sprintf('%s/%s/%s/%s/processed_data_gamma_Dec5.mat',InDir,Groups{G},FlName,Devices{G});
        load(Infile);
        
         psdout1OFF{G} = cat(2, psdout1OFF{G}(:,:), Data.psd0);
         psdout3OFF{G} = cat(2, psdout3OFF{G}(:,:), Data.psd2);
         psdout4OFF{G} = cat(2, psdout4OFF{G}(:,:), Data.psd3);
         
         coherence1OFF{G} = cat(2,coherence1OFF{G}(:,:),Data.cxy);
         coherence2OFF{G} = cat(2,coherence2OFF{G}(:,:),Data.cyz);
         coherence3OFF{G} = cat(2,coherence3OFF{G}(:,:),Data.czx);        

         catch
            fprintf('\n error \n');
            continue
        end
    end   
    
%    Time-Frequency Spectograms
    
%     figure;pspectrum(Specout1.ON,250,'spectrogram','FrequencyLimits',[1 120]);shading('interp');caxis([-80 -40]); title(sprintf('Stim:ON Channel:1  recorded over %d hours in Group: %s ',fix(Totalminutes/60),Groups{G}));xlabel('time');
%     figure;pspectrum(Specout2.ON,250,'spectrogram','FrequencyLimits',[1 120]);shading('interp'); caxis([-80 -40]); title(sprintf('Stim:ON Channel:2 recorded over %d hours in Group: %s ',fix(Totalminutes/60),Groups{G}));xlabel('time');
%     figure;pspectrum(Specout3.ON,250,'spectrogram','FrequencyLimits',[1 120]);shading('interp'); caxis([-80 -40]); title(sprintf('Stim:ON Channel:3 recorded over %d hours in Group: %s ',fix(Totalminutes/60),Groups{G}));xlabel('time');
%     figure;pspectrum(Specout4.ON,250,'spectrogram','FrequencyLimits',[1 120]);shading('interp'); caxis([-80 -40]); title(sprintf('Stim:ON Channel:4 recorded over %d hours in Group: %s ',fix(Totalminutes/60),Groups{G}));xlabel('time');

end

%% Coherence plotting

figure; plot(mean(coherence1ON{G}(:,:),2),'linewidth',2);xlim([0 90]); ylim([0 .1]); 
hold on;box off;
plot(mean(coherence1OFF{G}(:,:),2),'linewidth',2);xlim([0 90]); ylim([0 .1]);

figure; plot(mean(coherence2ON{G}(:,:),2),'linewidth',2);xlim([0 90]);ylim([0 .1]); 
hold on; box off;
plot(mean(coherence2OFF{G}(:,:),2),'linewidth',2);xlim([0 90]);ylim([0 .1]); 

figure; plot(mean(coherence3ON{G}(:,:),2),'linewidth',2);xlim([0 90]);ylim([0 .1]); 
hold on; box off;
plot(mean(coherence3OFF{G}(:,:),2),'linewidth',2);xlim([0 90]);ylim([0 .1]); 


% title(sprintf('Coherence between STN & C1 Group: %s',Groups{G}));
%% PSD Plotting

options=struct;  % Use once is for initial plot
options.handle=figure('Position',[800 400 800 400]); %Use once is for initial plot

c =5;  % change according to plot

switch c   
    case 1
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
    case 2
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
    case 3
options.color_area = rgb('DarkSalmon'); %red
options.color_line = rgb('Red');
    case 4
options.color_area = rgb('YellowGreen');    % yellow
options.color_line = rgb('Green');
    case 5
options.color_area = rgb('white');    % FOR OFF VALUES
options.color_line = rgb('Black');
end

options.line_width = 1.5;
options.x_axis= Data.h1(:,1);
options.error='sem';
options.alpha= .5;

plot_areaerrorbar(log10(psdout4OFF{G}(:,:))',options)
xlim([0 90]);
set(gca,'fontweight','bold','fontsize',24);
hold on
box off
% ylabel('Power (log_1_0\muV^2/Hz)');
%title(sprintf('PSD Stim:ON Group: %s',Groups{G}));   %change
clear psdout1ON psdout1OFF psdout4ON psdout3ON   psdout3OFF psdout4OFF 

% legend('','OFF','','.1A','','.3-.5A','','.6A-.7A','','.8A-1.1A') % LEFT
% legend('','.1A','','.3-.5A','','.6A-.8A','','OFF') % RIGHT

%% Transparent Graph plotting

figure;
for i=1: size(psdout1ON{1,G},2)
p=plot(log10(psdout1ON{G}(:,i)),'b-');
p.Color(4) = .1;
hold on
end

figure;
for i=1: size(psdout3ON{1,G},2)
p=plot(log10(psdout3ON{G}(:,i)),'b-');
p.Color(4) = .1;
hold on
end

figure;
for i=1: size(psdout4ON{1,G},2)
p=plot(log10(psdout4ON{G}(:,i)),'b-');
p.Color(4) = .1;
hold on
end

%% Figure 2 All OFF plots

figure('Position',[800 400 800 400]); 
plot(mean(log10(psdout4OFF{G}(:,:))')','b','LineWidth',1.5);
xlim([0 90]);ylim([-8 -4]);
%hline = findobj(gcf, 'type', 'line');set(hline(1),'LineStyle','.-') %change LineStyle here accordingly
set(gca,'fontweight','bold','fontsize',24);
hold on
box off



figure('Position',[800 400 800 400]) ; 
plot(mean(coherence2OFF{G}(:,:),2),'m','LineWidth',1.5);xlim([0 90]);ylim([0 .1]);
%hline = findobj(gcf, 'type', 'line');set(hline(1),'LineStyle','-.') %change LineStyle here accordingly
set(gca,'fontweight','bold','fontsize',24);
hold on; box off;


%% Figure 4

figure('Position',[800 400 800 400]) ; 
plot(mean(coherence2ON{G}(:,:),2),'b','LineWidth',1.5);xlim([0 90]);ylim([0 .1]);
% hline = findobj(gcf, 'type', 'line');set(hline(1),'LineStyle','-.') %change LineStyle here accordingly
set(gca,'fontweight','bold','fontsize',24);
hold on; box off;






 %title(sprintf('Power Spectrum Density Channel 1 Stim:ON Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Power (log_1_0\muV^2/Hz)');
%figure; plot(log10(mean(psdout1ON{G}(:,:),2)'));

%% Whole Session PLotting

for G=1:length(Groups)
figure; plot(Data.h1,log10(mean(fftout1.ON{G}(:,:),1)'),'linewidth',2);xlim([0 100]); title(sprintf('Power Spectrum Density Channel 1 Stim:ON Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Power (log_1_0\muV^2/Hz)');
figure; plot(Data.h1,log10(mean(fftout2.ON{G}(:,:),1)'),'linewidth',2);xlim([0 100]); title(sprintf('Power Spectrum Density Channel 2 Stim:ON Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Power (log_1_0\muV^2/Hz)');
figure; plot(Data.h1,log10(mean(fftout3.ON{G}(:,:),1)'),'linewidth',2);xlim([0 100]); title(sprintf('Power Spectrum Density Channel 3 Stim:ON Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Power (log_1_0\muV^2/Hz)');
figure; plot(Data.h1,log10(mean(fftout4.ON{G}(:,:),1)'),'linewidth',2);xlim([0 100]); title(sprintf('Power Spectrum Density Channel 4 Stim:ON Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Power (log_1_0\muV^2/Hz)');

figure; plot(Data.h1,log10(mean(fftout1.OFF{G}(:,:),1)'),'r');xlim([0 100]); title(sprintf('Power Spectrum Density Channel 1 Stim:OFF Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Power (log_1_0\muV^2/Hz)');
figure; plot(Data.h1,log10(mean(fftout2.OFF{G}(:,:),1)'),'r');xlim([0 100]); title(sprintf('Power Spectrum Density Channel 2 Stim:OFF Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Power (log_1_0\muV^2/Hz)');
figure; plot(Data.h1,log10(mean(fftout3.OFF{G}(:,:),1)'),'r');xlim([0 100]); title(sprintf('Power Spectrum Density Channel 3 Stim:OFF Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Power (log_1_0\muV^2/Hz)');
figure; plot(Data.h1,log10(mean(fftout4.OFF{G}(:,:),1)'),'r');xlim([0 100]); title(sprintf('Power Spectrum Density Channel 4 Stim:OFF Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Power (log_1_0\muV^2/Hz)');
end

%% Ttests

% Figure 4

% Coherence Stats

for s=1:90
[h1 p1 c1 stats1] = ttest2(coherence1ON{G}(s,:),coherence1OFF{G}(s,:)); 
pval1(s,1)=p1*90;
tval1(s,1)=stats1.tstat;
[h2 p2 c2 stats2] = ttest2(coherence2ON{G}(s,:),coherence2OFF{G}(s,:)); 
pval2(s,1)=p2*90;
tval2(s,1)=stats2.tstat;
end

% Figure 3

% stats

for s=1:90
[h7 p7 c7 stats7] = ttest2(log10(psdout1ON{G}(s,:)),(log10(psdout1OFF{G}(s,:))));
pval7(s,1)=p7*90;
tval7(s,1)=stats7.tstat;
[h8 p8 c8 stats8] = ttest2(log10(psdout3ON{G}(s,:)),(log10(psdout3OFF{G}(s,:))));
pval8(s,1)=p8*90;
tval8(s,1)=stats8.tstat;
[h9 p9 c9 stats9] = ttest2(log10(psdout4ON{G}(s,:)),(log10(psdout4OFF{G}(s,:))));
pval9(s,1)=p9*90;
tval9(s,1)=stats9.tstat;
end

% Figure 2

% stats

for s=1:90
[h10 p10 c10 stats10] = ttest2(log10(psdout1OFF{G}(s,:)),(log10(psdout3OFF{G}(s,:))));
pval10(s,1)=p10*90;
tval10(s,1)=stats10.tstat;
[h11 p11 c11 stats11] = ttest2(log10(psdout1OFF{G}(s,:)),(log10(psdout4OFF{G}(s,:))));
pval11(s,1)=p11*90;
tval11(s,1)=stats11.tstat;
[h12 p12 c12 stats12] = ttest2(log10(psdout3OFF{G}(s,:)),(log10(psdout4OFF{G}(s,:))));
pval12(s,1)=p12*90;
tval12(s,1)=stats12.tstat;
end

% Coherence Stats

for s=1:90
[h13 p13 c13 stats13] = ttest2(coherence1OFF{G}(s,:),coherence2OFF{G}(s,:)); 
pval13(s,1)=p13*90;
tval13(s,1)=stats13.tstat;
end
