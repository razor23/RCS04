%% Analyse Dystonia data %%

clear all; close all; clc

% Load up the data folder with all the files in %

InDir = '/Volumes/Vinith2/Starr/Dystonia/RawData';
Groups = {'RCS04L','RCS04R'};
Devices = {'DeviceNPC700418H','DeviceNPC700412H'};

%% Load file

for G=1:length(Groups)
    load(sprintf('%s/%s/StimSess_final.mat',InDir,Groups{G}));
    ONfiles  = StimSess(:,1);
    OFFfiles = StimSess(:,2);
    ONfiles =  ONfiles(~cellfun('isempty',ONfiles)); % there are some empty cells, geting rid of em
    OFFfiles = OFFfiles(~cellfun('isempty',OFFfiles)); % there are some empty cells, geting rid of em

    psdout1ON{G} = [];  psdout1OFF{G} = [];
    psdout2ON{G} = [];  psdout2OFF{G} = [];
    psdout3ON{G} = [];  psdout3OFF{G} = [];
    psdout4ON{G} = [];  psdout4OFF{G} = [];
    
    coherence1ON{G} = [];   coherence1OFF{G} = [];
    coherence2ON{G} = [];   coherence2OFF{G} = [];
    coherence3ON{G} = [];   coherence3OFF{G} = [];


    % ON files data
 
fprintf('\n Printing ON files \n');

    
    for f=1:19
        
        try
        FlName = ONfiles{f};
        fprintf('\n Loading file %s from %s',FlName,Groups{G});
        path = sprintf('%s/%s/%s/%s',InDir,Groups{G},FlName,Devices{G});
        Infile = sprintf('%s/%s/%s/%s/processed_data_final.mat',InDir,Groups{G},FlName,Devices{G});
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
        Infile = sprintf('%s/%s/%s/%s/processed_data_final.mat',InDir,Groups{G},FlName,Devices{G});
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

figure; plot(mean(coherence1ON{G}(:,:),2),'linewidth',2);xlim([0 60]); ylim([0 .15]); title(sprintf('Coherence between STN & M1 Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Coherence');
hold on
plot(mean(coherence1OFF{G}(:,:),2),'linewidth',2);xlim([0 60]); ylim([0 .15]);title(sprintf('Coherence between STN & M1 Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Coherence'); legend('on','off')

figure; plot(mean(coherence2ON{G}(:,:),2),'linewidth',2);xlim([0 60]);ylim([0 .15]); title(sprintf('Coherence between STN & M2 Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Coherence');
hold on
plot(mean(coherence2OFF{G}(:,:),2),'linewidth',2);xlim([0 60]);ylim([0 .15]); title(sprintf('Coherence between STN & M2 Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Coherence'); legend('on','off')

figure; plot(mean(coherence3ON{G}(:,:),2),'linewidth',2);xlim([0 60]);ylim([0 .15]); title(sprintf('Coherence between M1 & M2 Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Coherence');
hold on
plot(mean(coherence3OFF{G}(:,:),2),'linewidth',2);xlim([0 60]);ylim([0 .15]); title(sprintf('Coherence between M1 & M2 Group: %s',Groups{G}));xlabel('Frequency (Hz)'); ylabel('Coherence'); legend('on','off')

%% PSD Plotting

options=struct;  %Just once is fine
options.handle=figure('Position',[800 400 800 400]); %Just once is fine


% options.color_area = [128 193 219]./255;    % Blue theme
% options.color_line = [ 52 148 186]./255;
% options.color_area = [243 169 114]./255;    % Orange theme
% % options.color_line = [236 112  22]./255;
% options.color_area = rgb('DarkSalmon'); %red
% options.color_line = rgb('Red');
options.color_area = rgb('YellowGreen');    % yellow
options.color_line = rgb('YellowGreen');

options.line_width = 1;
options.x_axis= Data.h1(:,1);
options.error='sem';
options.alpha= 0.5;

plot_areaerrorbar(log10(psdout1ON{G}(:,:))',options)
xlabel('Frequency(hz)','Fontsize',10,'Fontname','Arial');
ylabel('Power (log_1_0\muV^2/Hz)');
title(sprintf('PSD Stim:ON Group: %s',Groups{G}));   %change
hold on
box off



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

%LEFT
[h1 p1] = ttest2(mean(coherence1ON{1}([4:7],:),1),mean(coherence1OFF{1}([4:7],:),1)); % STN M1
[h2 p2] = ttest2(mean(coherence2ON{1}([4:7],:),1),mean(coherence2OFF{1}([4:7],:),1)); % STN M2
[h3 p3] = ttest2(mean(coherence3ON{1}([4:7],:),1),mean(coherence3OFF{1}([4:7],:),1)); % M1(+9-8) M2 (+11-10)

%RIGHT
[h4 p4] = ttest(mean(coherence1ON{2}([4:7],:),2),mean(coherence1OFF{2}([4:7],:),2)); %STN M1
[h5 p5] = ttest(mean(coherence2ON{2}([4:7],:),2),mean(coherence2OFF{2}([4:7],:),2)); %STN M2
[h6 p6] = ttest(mean(coherence3ON{2}([4:7],:),2),mean(coherence3OFF{2}([4:7],:),2)); % M1(+9-8) M2 (+11-10)


[h7 p7] = ttest(mean(log10(psdout1ON{G}([4:7],:)),2),mean(log10(psdout1ON{G}([4:7],:)),2));
[h8 p8] = ttest(mean(log10(psdout1ON{G}([4:7],:)),2),mean(log10(psdout1ON{G}([4:7],:)),2));
[h9 p9] = ttest(mean(log10(psdout1ON{G}([4:7],:)),2),mean(log10(psdout1ON{G}([4:7],:)),2));

