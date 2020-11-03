%% Beta Burst %%

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
    OFFfiles = OFFfiles(~cellfun('isempty',OFFfiles));
    thrsh1{G}=[];thrsh2{G}=[];thrsh3{G}=[];
    
    fprintf('\n Printing OFF files \n');
    
    for f=1:length(OFFfiles)
        try
            FlName = OFFfiles{f};
            fprintf('\n Loading file %s from %s',FlName,Groups{G});
            path = sprintf('%s/%s/%s/%s',InDir,Groups{G},FlName,Devices{G});
            Infile = sprintf('%s/%s/%s/%s/processed_data_final.mat',InDir,Groups{G},FlName,Devices{G});
            load(Infile);
            
            %% Plot the theta range signal
            
            [b,a] = butter(4,[3,7]/(250/2));
            temp1outTheta=filtfilt(b,a,Data.temp0);
            temp2outTheta=filtfilt(b,a,Data.temp2);
            temp3outTheta=filtfilt(b,a,Data.temp3);
            
            t1amp=abs(hilbert(temp1outTheta));
            t2amp=abs(hilbert(temp2outTheta));
            t3amp=abs(hilbert(temp3outTheta));
            
            Thresh1=prctile(t1amp,75);
            Thresh2=prctile(t2amp,75);
            Thresh3=prctile(t3amp,75);
            
            thrsh1{G} = cat(1,thrsh1{G},Thresh1);
            thrsh2{G} = cat(1,thrsh2{G},Thresh2);
            thrsh3{G} = cat(1,thrsh3{G},Thresh3);
            
             clear t1amp t2amp t3amp
        catch
            fprintf('\n error \n');
             clear t1amp t2amp t3amp
           continue
        end
    end
end

clear Data
%% Get the threshold values

T1=mean(thrsh1{1}(:,:)); %75th Percentile value OFF files LEFT STN1
T2=mean(thrsh2{1}(:,:)); %75th Percentile value OFF files LEFT M1
T3=mean(thrsh3{1}(:,:)); %75th Percentile value OFF files LEFT M2
T4=mean(thrsh1{2}(:,:)); %75th Percentile value OFF files RIGHT STN1
T5=mean(thrsh2{2}(:,:)); %75th Percentile value OFF files RIGHT M1
T6=mean(thrsh3{2}(:,:)); %75th Percentile value OFF files RIGH  M2












