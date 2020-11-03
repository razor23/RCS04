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

%% Run loops to get 
 
    for f=1:length(OFFfiles)     
        try
            FlName = OFFfiles{f};
            fprintf('\n Loading file %s from %s',FlName,Groups{G});
            path = sprintf('%s/%s/%s/%s',InDir,Groups{G},FlName,Devices{G});
            Infile = sprintf('%s/%s/%s/%s/processed_data_final.mat',InDir,Groups{G},FlName,Devices{G});
            load(Infile);
            
            [b,a] = butter(4,[3,7]/(250/2));
            s=1;
            
            for u=1:floor(size(Data.temp0,1)/(60*250))
                
            l= 250*60*u;
                                          
            temp1outTheta=filtfilt(b,a,Data.temp0(s:l));
            temp2outTheta=filtfilt(b,a,Data.temp2(s:l));
            temp3outTheta=filtfilt(b,a,Data.temp3(s:l));
            
            s=l+1;
            
            t1amp=abs(hilbert(temp1outTheta));
            t2amp=abs(hilbert(temp2outTheta));
            t3amp=abs(hilbert(temp3outTheta));
            
            D_amp1 = t1amp;
            D_amp2 = t2amp;
            D_amp3 = t3amp;
            
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
            idxl = tmpS>=T1;
            idxl(1) = 0;
            idx = find(idxl);
            yest = tmpS(idx-1)<T1;
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
                
                if ~isempty(find(D_amp1(st:end)<T1))
                    edt=find(D_amp1(st:end)<T1);
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
                BurstLength(k)=ed1k(k)-st1k(k);                           
            end       
            
           numburst = length(st1k);
           datasegmentsecs = length(D_amp1)/250;
           Burstrate = numburst/datasegmentsecs;
            
            % Sanity check - plot the beginning of each burst %
            figure;
            plot(D_amp1);
            hold on;
            % Put on the start of each burst
            zA2=zA;
            zA2(st1k)=1*T1;
            plot(zA2,'g');
            % Put on the end of each burst
            zA3=zA;
            zA3(ed1k)=1*T1;
            plot(zA3,'k');
            ThreshShow=T1*ones(size(D_amp1,1),size(D_amp1,2));
            plot(ThreshShow,'r');
            
            end
        end
    end













