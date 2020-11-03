%% Analyse Dystonia data %%

clear all; close all; clc

% Load up the data folder with all the files in %

InDir = '/Volumes/Vinith2/Starr/Dystonia/RawData';
Groups = {'RCS04L','RCS04R'};
Devices = {'DeviceNPC700418H','DeviceNPC700412H'};
Date1 = datetime('13-Jul-2019'); %set the day you want to start the analysis from
Date4 = datetime('30-Jul-2019'); %visit dates
Date5 = datetime('29-Aug-2019'); %visit dates

ch1 = '+2-0 lpf1-100Hz lpf2-1700Hz sr-250Hz';
ch2 = '+3-1 lpf1-450Hz lpf2-1700Hz sr-250Hz';
ch3 = '+9-8 lpf1-450Hz lpf2-1700Hz sr-250Hz';
ch4 = '+11-10 lpf1-450Hz lpf2-1700Hz sr-250Hz';


High =1;
filterOrder = 3*fix(250/High); %250 being SR
fs =250;

%% Load file

for G=1:length(Groups)
    load(sprintf('%s/%s/database_from_device_settings.mat',InDir,Groups{G}));

    for f = 1:length(masterTableUse.session)
        FlName = masterTableUse.session{f};
        fprintf('\n \n Loading file %s from %s',FlName,Groups{G});
        path = sprintf('%s/%s/%s/%s',InDir,Groups{G},FlName,Devices{G});
        Infile = sprintf('%s/%s/%s/%s/RawDataTD.mat',InDir,Groups{G},FlName,Devices{G});
        try            
            load(Infile)
            if masterTableUse.senseSettings{f,1}.samplingRate == 250 

% Trimmed data

temp1outTr=temp1raw(fs*10:fs*2000);
temp2outTr=temp2raw(fs*10:fs*2000);
temp3outTr=temp3raw(fs*10:fs*2000);

% Here I was just looking at beta and seeing if the beta amplitude
% correlates across the two cortical sites.
[b,a] = butter(4, [13,30]/(fs/2));
temp2out=filtfilt(b,a,temp2outTr);
temp3out=filtfilt(b,a,temp3outTr);

% So there is a correlation between the beta amplitudes across the cortical
% contacts - but it is very low.
figure; scatter(abs(hilbert(temp2outTr)), abs(hilbert(temp3outTr)));
[rho,p]=corr(abs(hilbert(temp2outTr)), abs(hilbert(temp3outTr)))

%%

% So things to do next would be:
% 1. Run the coherence on the artifact rejected data BUT - make sure that
% when you artifact reject in any channel - you do a matched rejection in
% the other channels - so that all the channels are always aligned and
% always of the same length.
% 2. See how it looks on the larger dataset - maybe it is low here - just
% because this is a weird piece of data - seems unlikely - but possible.
% 3. HOW ABOUT TRYING IF YOU JUST LOOK AT PERIODS WHEN THE CORTICAL THETA IS
% HIGH if the coherence goes up. So that would look something like this:

%%
% Plot the theta range signal
[b,a] = butter(4, [3,7]/(fs/2));
temp1outTheta=filtfilt(b,a,temp1outTr);
temp2outTheta=filtfilt(b,a,temp2outTr);
temp3outTheta=filtfilt(b,a,temp3outTr);

t1Amp=abs(hilbert(temp1outTheta));
t2Amp=abs(hilbert(temp2outTheta));
t3Amp=abs(hilbert(temp3outTheta));

% Plot the histogram of the theta amplitudes %
figure;
subplot(3,1,1);
histogram(t2Amp)
xlim([0 0.02])
title('Cortex 1');
subplot(3,1,2);
histogram(t3Amp)
xlim([0 0.02])
title('Cortex 2');
subplot(3,1,3);
histogram(t1Amp)
title('Subcortical');
xlim([0 0.02])

% Having eyeballs the histogra -
% Use a cortical theta threshold of 0.01 to find just when the theta is
% a bit higher
whr=t3Amp>0.02;
% Smooth out a bit to find longer periods increased theta
whrSm=smooth(whr,250);
figure; plot(whrSm)

whrTheta=find(whrSm>0.3);

T1ThetaEnriched=temp1outTr(whrTheta);
T2ThetaEnriched=temp2outTr(whrTheta);
T3ThetaEnriched=temp3outTr(whrTheta);

[cxyThetaEnrich1,f] = mscohere(T1ThetaEnriched,T2ThetaEnriched,fs*5,round(fs/2),[1:90],fs);
[cxyThetaEnrich2,f] = mscohere(T1ThetaEnriched,T3ThetaEnriched,fs*5,round(fs/2),[1:90],fs);

figure;
plot(f,cxyThetaEnrich1);
hold on;
plot(f,cxyThetaEnrich2);

% SO THIS SEEMS TO SHOW - THAT IF YOU ENRICH THE DATA - AND JUST LOOK WHEN THETA IS HIGHER - YOU DO SEE
% COHERENCE PEAKS THAT LOOK A BIT MORE EXPECTED.

% I WOULD RUN SOMETHING LIKE THIS ON THE ALIGNED ARTIFACT REJECTED DATA -
% WITH AND WITHOUT THETA ENRICHMENT AND SEE WHAT IT SHOWS

%% THETA CROSS CORRELATION AND BURST ANALYSIS

% Here are our theta filtered data
% Get the amplitudes

t1amp=abs(hilbert(temp1outTheta));
t2amp=abs(hilbert(temp2outTheta));
t3amp=abs(hilbert(temp3outTheta));

% Now before you do the cross correlation - always mean subtract:
t1ampMs=t1amp-mean(t1amp);
t2ampMs=t2amp-mean(t2amp);
t3ampMs=t3amp-mean(t3amp);

[r,lags]=xcorr(t1ampMs,t3ampMs,'coeff');

% Time (get it into ms, 4ms per sample at FS 250Hz);
lagms=lags*4;
figure; plot(lagms,r);

%% NEXT - GET THE BURST FEATURES

% THIS IS WHERE YOU SET THE THRESHOLD.
% We'll set at the 75% percentile...
% BASELINE (SEE TINKHAUSER, LITTLE AND BROWN PAPER IN BRAIN).

% Starting with the subcortical version first:
D_amp1=t1amp';
% Use these instead for the cortical  versions
% D_amp2=t1amp
% D_amp2=t1amp


% Really important - that when you run on the ON stim state here - you use
% the threshold for the OFF stim data (so that the threshold is the same
% across all the different recordings. You might want to write a little
% script and to work out what the 75th percentile is across the whole OFF
% stim dataset and use that for both the OFF and ON stim analyses

% Here just working it out with the current session data - but make sure
% the threshold stays constant across all data analyses or the results will
% be seriously confounded.
thresh=prctile(D_amp1,75);

% Set initial variables %
pkAs=[];
zA=zeros(size(D_amp1,1),size(D_amp1,2));
tsLocs=[];
tsAmp=[];
psst=[];
BrL=[];
BrH=[];
% Set first value as 0 
D_amp1(1)=0;

% Find start of each burst
tmpS=D_amp1;
idxl = tmpS>=thresh;
idxl(1) = 0;
idx = find(idxl);
yest = tmpS(idx-1)<thresh;
st1k=idx(yest);



% Set some more variables
Mx=[];

ed1k=[];
st2k=[];
ed2k=[];
InaS=[];

% Analyse each burst in turn %
for k=1:length(st1k)
    st=st1k(k);
    
    if ~isempty(find(D_amp1(st:end)<thresh))
        edt=find(D_amp1(st:end)<thresh);
        % -2 fence post problem and also wanting point
        % above threshold (inside)
        ed=edt(1)+st-2;
    else
        ed=size(D_amp1,2);
    end
    ed1k(k)=ed;
    
    % Get max peak height and location
    amp=D_amp1(st:ed);
    [BurstHeight(k),In]=max(amp);
    Ina=In+st-1;
    BurstPeakPosition(k)=Ina;
    BurstLength(k)=ed1k(k)-st1k(k);
    
end

% Sanity check - plot the beginning of each burst %
figure; 
plot(D_amp1); 
hold on; 
% Put on the start of each burst
zA2=zA;
zA2(st1k)=1*thresh;
plot(zA2,'g');
% Put on the end of each burst
zA3=zA;
zA3(ed1k)=1*thresh;
plot(zA3,'k');
ThreshShow=thresh*ones(size(D_amp1,1),size(D_amp1,2));
plot(ThreshShow,'r')



            end
        end
    end
end






