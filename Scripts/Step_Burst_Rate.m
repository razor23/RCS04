%% Burst Rate

clear all; close all; clc

% Load up the data folder with all the files in %

InDir = '/Volumes/Vinith2/Starr/Dystonia/RawData';
Groups = {'RCS04L','RCS04R'};
Devices = {'DeviceNPC700418H','DeviceNPC700412H'};

T1L=.0102;
T2L=.0096;
T3L=.0143;

T1R=.0113;
T2R=.0140;
T3R=.0122;
totminON=0;
totminOFF=0;

%% Load files

%% OFF files

for G=1%:length(Groups)  % Change depenfind on LEFT or RIGHT
    load(sprintf('%s/%s/StimSess_final.mat',InDir,Groups{G}));
    ONfiles  = StimSess(:,1);
    OFFfiles = StimSess(:,2);
    ONfiles =  ONfiles(~cellfun('isempty',ONfiles)); % there are some empty cells, geting rid of em
    OFFfiles = OFFfiles(~cellfun('isempty',OFFfiles));
    BrateOFF{G}=[]; BrateON{G}=[];
    BlenOFF{G} = []; BlenON{G}=[];
    BhtOFF{G}=[]; BhtON{G}=[];
    
                fprintf('\n ***** OFF Files *****');

    
    for f=1:length(OFFfiles)
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
            
            %% SubCortical One Electrode (+2-0)
            
            % Set first value as 0
            D_amp1(1)=0;
            
            % Find start of each burst
            tmpS=D_amp1;
            idxl = tmpS>=T1L;         % CHANGE THE VALUE DEPENDING ON L OR R
            idxl(1) = 0;
            idx = find(idxl);
            yest = tmpS(idx-1)<T1L;    % CHANGE THE VALUE DEPENDING ON L OR R
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
                
                if ~isempty(find(D_amp1(st:end)<T1L))   % CHANGE THE VALUE DEPENDING ON L OR R
                    edt=find(D_amp1(st:end)<T1L); % CHANGE THE VALUE DEPENDING ON L OR R
                    % -2 fence post problem and also wanting point
                    % above threshold (inside)
                    ed=edt(1)+st-2;
                else
                    ed=size(D_amp1,2);
                end
                ed1k(k)=ed;
                
                  if st>ed
                    continue;
                else
                
                % Get max peak height and location
                amp=D_amp1(st:ed);
                [BurstHeight1(k),In]=max(amp);
                Ina=In+st-1;
                BurstLength1(k)=abs(ed1k(k)-st1k(k)); %sometimes negative value? Hence the Abs
                  end
            end
            
            numburst(u,1) = length(st1k);
            datasegmentsecs(u,1) = length(D_amp1)/250;
            Burstrate(u,1) = numburst(u,1)/datasegmentsecs(u,1);
            Blen(u,1) = mean(BurstLength1); %mean length for all bursts in a minute
            Bht(u,1) = mean(BurstHeight1);
            
            clear tmpS k In Ina st1k ed edt yest ed1k
            %% Cortical 1 Electrode (+8-9)
            
            % Set first value as 0
            D_amp2(1)=0;
            
            % Find start of each burst
            tmpS=D_amp2;
            idxl = tmpS>=T2L;  % CHANGE THE VALUE DEPENDING ON L OR R
            idxl(1) = 0;
            idx = find(idxl);
            yest = tmpS(idx-1)<T2L;  % CHANGE THE VALUE DEPENDING ON L OR R
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
                
                if ~isempty(find(D_amp2(st:end)<T2L))  % CHANGE THE VALUE DEPENDING ON L OR R
                    edt=find(D_amp2(st:end)<T2L);   % CHANGE THE VALUE DEPENDING ON L OR R
                    % -2 fence post problem and also wanting point
                    % above threshold (inside)
                    ed=edt(1)+st-2;
                else
                    ed=size(D_amp2,2);
                end
                ed1k(k)=ed;
                
                 if st>ed
                    continue;
                else
                
                % Get max peak height and location
                amp=D_amp2(st:ed);
               [BurstHeight2(k),In]=max(amp);
               Ina=In+st-1;
                BurstLength2(k)=abs(ed1k(k)-st1k(k)); %sometimes negative value? Hence the Abs

                 end
            end
            
            numburst(u,2) = length(st1k);
            datasegmentsecs(u,2) = length(D_amp2)/250;
            Burstrate(u,2) = (numburst(u,2)/datasegmentsecs(u,2));
            Blen(u,2) = mean(BurstLength2); %mean length for all bursts in a minut
            Bht(u,2) = mean(BurstHeight2);


            clear tmpS k In Ina st1k ed edt yest ed1k
            
            %% Cortical 2 Electrode (+10-11)
            
            % Set first value as 0
            D_amp3(1)=0;
            
            % Find start of each burst
            tmpS=D_amp3;
            idxl = tmpS>=T3L;  % CHANGE THE VALUE DEPENDING ON L OR R
            idxl(1) = 0;
            idx = find(idxl);
            yest = tmpS(idx-1)<T3L;  % CHANGE THE VALUE DEPENDING ON L OR R
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
                
                if ~isempty(find(D_amp3(st:end)<T3L))  % CHANGE THE VALUE DEPENDING ON L OR R
                    edt=find(D_amp3(st:end)<T3L);  % CHANGE THE VALUE DEPENDING ON L OR R
                    % -2 fence post problem and also wanting point
                    % above threshold (inside)
                    ed=edt(1)+st-2;
                else
                    ed=size(D_amp3,2);
                end
                ed1k(k)=ed;
                
                
                if st>ed
                    continue;
                else
                 
                % Get max peak height and location                
                amp=D_amp3(st:ed);                  
                [BurstHeight3(k),In]=max(amp);
                Ina=In+st-1;
                BurstLength3(k)=abs(ed1k(k)-st1k(k));
                end
            end
            
            numburst(u,3) = length(st1k);
            datasegmentsecs(u,3) = length(D_amp3)/250;
            Burstrate(u,3) = numburst(u,3)/datasegmentsecs(u,3);
            Blen(u,3) = mean(BurstLength3); %mean length for all bursts in a minut
            Bht(u,3) = mean(BurstHeight3);        
            clear tmpS k Ina In st1k ed edt yest ed1k
        end %minute
        clear Data
        totminOFF = totminOFF+u;
       BrateOFF{G}= cat(1,BrateOFF{G}(:,:),Burstrate);
       BlenOFF{G}=cat(1,BlenOFF{G}(:,:),Blen);
       BhtOFF{G}=cat(1,BhtOFF{G}(:,:),Bht);
    end % file
    %BR = 
end% Group

clear Blen Burstrate 
%% ON files

for G=1%:length(Groups)
    load(sprintf('%s/%s/StimSess_final.mat',InDir,Groups{G}));
    ONfiles  = StimSess(:,1);
    OFFfiles = StimSess(:,2);
    ONfiles =  ONfiles(~cellfun('isempty',ONfiles)); % there are some empty cells, geting rid of em
    OFFfiles = OFFfiles(~cellfun('isempty',OFFfiles));
    BrateON{G}=[];
    
            fprintf('\n ***** ON Files *****');
            
          
                
                for f=[1:3 9 10 15]
                    
                      try
                    FlName = ONfiles{f};
                    fprintf('\n Loading file %s from %s',FlName,Groups{G});
                    path = sprintf('%s/%s/%s/%s',InDir,Groups{G},FlName,Devices{G});
                    Infile = sprintf('%s/%s/%s/%s/processed_data_final.mat',InDir,Groups{G},FlName,Devices{G});
                    load(Infile);
                    
                    catch
                        fprintf('Error opening file %s',FlName);
                        continue;
             
                      end
        
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
            
            %% SubCortical One Electrode (+2-0)
            
            % Set first value as 0
            D_amp1(1)=0;
            
            % Find start of each burst
            tmpS=D_amp1;
            idxl = tmpS>=T1L;         % CHANGE THE VALUE DEPENDING ON L OR R
            idxl(1) = 0;
            idx = find(idxl);
            yest = tmpS(idx-1)<T1L;    % CHANGE THE VALUE DEPENDING ON L OR R
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
                
                if ~isempty(find(D_amp1(st:end)<T1L))   % CHANGE THE VALUE DEPENDING ON L OR R
                    edt=find(D_amp1(st:end)<T1L); % CHANGE THE VALUE DEPENDING ON L OR R
                    % -2 fence post problem and also wanting point
                    % above threshold (inside)
                    ed=edt(1)+st-2;
                else
                    ed=size(D_amp1,2);
                end
                ed1k(k)=ed;
                
                     
                if st>ed
                    continue;
                else
                
                % Get max peak height and location
                amp=D_amp1(st:ed);
                [BurstHeight1(k),In]=max(amp);
                Ina=In+st-1;
                BurstLength1(k)=abs(ed1k(k)-st1k(k));
                end
            end
            
            numburst(u,1) = length(st1k);
            datasegmentsecs(u,1) = length(D_amp1)/250;
            Burstrate(u,1) = numburst(u,1)/datasegmentsecs(u,1);
            Blen(u,1) = mean(BurstLength1); %mean length for all bursts in a minute
            Bht(u,1) = mean(BurstHeight1);

            
            clear tmpS k In Ina st1k ed edt yest ed1k
            %% Cortical 1 Electrode (+8-9)
            
            % Set first value as 0
            D_amp2(1)=0;
            
            % Find start of each burst
            tmpS=D_amp2;
            idxl = tmpS>=T2L;  % CHANGE THE VALUE DEPENDING ON L OR R
            idxl(1) = 0;
            idx = find(idxl);
            yest = tmpS(idx-1)<T2L;  % CHANGE THE VALUE DEPENDING ON L OR R
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
                
                if ~isempty(find(D_amp2(st:end)<T2L))  % CHANGE THE VALUE DEPENDING ON L OR R
                    edt=find(D_amp2(st:end)<T2L);   % CHANGE THE VALUE DEPENDING ON L OR R
                    % -2 fence post problem and also wanting point
                    % above threshold (inside)
                    ed=edt(1)+st-2;
                else
                    ed=size(D_amp2,2);
                end
                ed1k(k)=ed;
                
                % Get max peak height and location
                     
                if st>ed
                    continue;
                else
                amp=D_amp2(st:ed);
                [BurstHeight2(k),In]=max(amp);
                Ina=In+st-1;
                BurstLength2(k)=abs(ed1k(k)-st1k(k)); %negatice value sometimes
                end
            end
            
            Blen(u,2) = mean(BurstLength2); %mean length for all bursts in a minute

            numburst(u,2) = length(st1k);
            datasegmentsecs(u,2) = length(D_amp2)/250;
            Burstrate(u,2) = (numburst(u,2)/datasegmentsecs(u,2));
            Bht(u,2) = mean(BurstHeight2);

            clear tmpS k In Ina st1k ed edt yest ed1k
            
            %% Cortical 2 Electrode (+10-11)
            
            % Set first value as 0
            D_amp3(1)=0;
            
            % Find start of each burst
            tmpS=D_amp3;
            idxl = tmpS>=T3L;  % CHANGE THE VALUE DEPENDING ON L OR R
            idxl(1) = 0;
            idx = find(idxl);
            yest = tmpS(idx-1)<T3L;  % CHANGE THE VALUE DEPENDING ON L OR R
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
                
                if ~isempty(find(D_amp3(st:end)<T3L))  % CHANGE THE VALUE DEPENDING ON L OR R
                    edt=find(D_amp3(st:end)<T3L);  % CHANGE THE VALUE DEPENDING ON L OR R
                    % -2 fence post problem and also wanting point
                    % above threshold (inside)
                    ed=edt(1)+st-2;
                else
                    ed=size(D_amp3,2);
                end
                ed1k(k)=ed;
                
                % Get max peak height and location
                     
                if st>ed
                    continue;
                else
                amp=D_amp3(st:ed);
                [BurstHeight3(k),In]=max(amp);
                Ina=In+st-1;
                BurstLength3(k)=abs(ed1k(k)-st1k(k));
                end
            end
            
            numburst(u,3) = length(st1k);
            datasegmentsecs(u,3) = length(D_amp3)/250;
            Burstrate(u,3) = numburst(u,3)/datasegmentsecs(u,3);
            Blen(u,3) = mean(BurstLength3);
            Bht(u,3) = mean(BurstHeight3);

            clear tmpS k Ina In st1k ed edt yest ed1k
        end %minute        
        clear Data
        totminON = totminON+u;
       BrateON{G}= cat(1,BrateON{G}(:,:),Burstrate);
       BlenON{G}=cat(1,BlenON{G}(:,:),Blen);
       BhtON{G}=cat(1,BhtON{G}(:,:),Bht);

       
    end % file
    %BR = 
end% Group

%% Stats

% Figure 4

[h11 p11 c1 stats1] = ttest2((BrateOFF{1,1}(:,1)),(BrateON{1,1}(:,1))); %LEFT
[h21 p21 c2 stats2] = ttest2((BrateOFF{1,1}(:,2)),(BrateON{1,1}(:,2))); %LEFT
[h31 p31 c3 stats3 ] = ttest2((BrateOFF{1,1}(:,3)),(BrateON{1,1}(:,3))); %LEFT

[h41 p41 c4 stats4] = ttest2((BrateOFF{1,2}(:,1)),(BrateON{1,2}(:,1))); %RIGHT
[h51 p51 c5 stats5] = ttest2((BrateOFF{1,2}(:,2)),(BrateON{1,2}(:,2))); %RIGHT
[h61 p61 c6 stats6] = ttest2((BrateOFF{1,2}(:,3)),(BrateON{1,2}(:,3))); %RIGHT






