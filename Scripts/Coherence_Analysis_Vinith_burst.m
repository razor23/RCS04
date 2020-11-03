clear all
close all

%% Load data

%% Change this to link to your data location %%
load('C:\Users\slittle\Box\Research\RCS04_data\Vinith_data\Coherence_sample_data.mat');

figure; plot(temp1raw); hold on; plot(temp2raw); title('Raw - data - see massive artifact at the end');

% Plot the PSD and Coherence for the cortical and STN data
fs=250;
% So if you run the analyis on the earlier data that doesn't have an
% artifact:
[cxy,f] = mscohere(temp1raw(fs*10:fs*2000),temp2raw(fs*10:fs*2000),fs*5,round(fs/2),[1:90],fs)
[cxy2,f] = mscohere(temp1raw(fs*10:fs*2000),temp3raw(fs*10:fs*2000),fs*5,round(fs/2),[1:90],fs)
[pxx1,f] = pwelch(temp1raw(fs*10:fs*2000),fs*5,fs,[1:90],fs)
[pxx2,f] = pwelch(temp2raw(fs*10:fs*2000),fs*5,fs,[1:90],fs)
[pxx3,f] = pwelch(temp3raw(fs*10:fs*2000),fs*5,fs,[1:90],fs)

figure; 
subplot(2,1,1); 
plot(f,pxx1); hold on; plot(f,pxx2); plot(f,pxx3);
title('PSDs - no artifacts');
subplot(2,1,2);
plot(f,cxy); 
hold on; plot(f,cxy2);
title('Coherence - no artifacts');

% Gets rid of that high frequency increased coherence (which must be due to
% that artifact). 

% Also you can show that this goes crazy if you don't remove the
% artifact....
fs=250;
[cxy,f] = mscohere(temp1raw,temp2raw,fs*5,round(fs/2),[1:90],fs)
[pxx1,f] = pwelch(temp1raw,fs*5,fs,[1:90],fs)
[pxx2,f] = pwelch(temp2raw,fs*5,fs,[1:90],fs)

figure; 
subplot(2,1,1); 
plot(f,pxx1); hold on; plot(f,pxx2); title('PSDs - with artifacts');
subplot(2,1,2);
plot(f,cxy); title('Coherence - without artifacts');

% SO THE HIGH FREQUENCY INCREASED COHERENCE IS INDEED DUE TO THE
% ARTIFACTS

% BUT THERE IS VERY LITTLE LOW FREQUENCY COHERENCE HERE %

% First sanity check: Check the algorithm
% Run the coherence on the same data and check that you get high coherence:
[cxycheck,f] = mscohere(temp1raw,temp1raw,fs*5,round(fs/2),[1:90],fs)
figure; plot(f,cxycheck);

% Comes out at 1 - so that means that the coherence algorithm is working...

%% SANITY CHECK 2 - what about looking at the coherence of the two cortical channels which are close to each other.
% We would expect to see some coherence between these two sites.
% Plot the cortical data %%

[cxyCor,f] = mscohere(temp2raw(fs*10:fs*2000),temp3raw(fs*10:fs*2000),fs*5,round(fs/2),[1:90],fs)
[pxx5,f] = pwelch(temp2raw(fs*10:fs*2000),fs*5,fs,[1:90],fs)
[pxx6,f] = pwelch(temp3raw(fs*10:fs*2000),fs*5,fs,[1:90],fs)

figure; 
subplot(2,1,1);
plot(f,pxx5); hold on; plot(f,pxx6)
subplot(2,1,2);
plot(f,cxyCor);

% This is surprisingly low! I would have thought that this would have been
% much higher....

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
whr=t1Amp>0.01;
% Smooth out a bit to find longer periods increased theta
whrSm=smooth(whr,250);
figure; plot(whrSm)

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

