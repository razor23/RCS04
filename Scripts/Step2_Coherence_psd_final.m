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
    z=1;x=1;k=1;
    daterror=0;
    lengtherror=0;
    channelerror=0;
    sleeperror=0;
    diffSR=0;
    filesprocessed=0;
    throwingerror=0;
    TotalONminutes=0;
    TotalOFFminutes=0;
    for f = 1:length(masterTableUse.session)
        FlName = masterTableUse.session{f};
        fprintf('\n \n Loading file %s from %s',FlName,Groups{G});
        path = sprintf('%s/%s/%s/%s',InDir,Groups{G},FlName,Devices{G});
        Infile = sprintf('%s/%s/%s/%s/RawDataTD.mat',InDir,Groups{G},FlName,Devices{G});
        try            
            load(Infile)
            if masterTableUse.senseSettings{f,1}.samplingRate == 250 % onlylooking at 250 data
                
                % Date check
                
                d=split(datestr(masterTableUse.senseSettings{f,1}.time));
                Date2 = datetime(d(1));
                Date = Date2-Date1;
                
                if Date<0
                    fprintf('\n This file date is not valid');
                    daterror=daterror+1;
                    continue;
                end
               
                if Date2 == Date4
                    fprintf('\n This file date is not valid');
                    daterror=daterror+1;
                    continue;
                end
               
                if Date2 == Date5
                    fprintf('\n This file date is not valid');
                    daterror=daterror+1;
                    continue;
                end
               
                %time check (only awake times)
                
                timecheck1 = masterTableUse.timeStart(f,1); %start time
                timecheck2 = masterTableUse.timeEnd(f,1); %endtime
                
                if ismember(hour(timecheck1),[23 0:7])  %start time
                    fprintf('\n patient was sleeping');
                    sleeperror=sleeperror+1; %sleep count
                    continue;
                end
                
                if ismember(hour(timecheck2),[0:7]) %endtime
                    fprintf('\n patient was sleeping');
                    sleeperror=sleeperror+1; %sleep count
                    continue;
                end
                
                % length check (at least one minute recd i.e >15000 samples)
                
                if size(outdatcomplete,1)<15000
                    fprintf('\t This file size is tooo small \n ');
                    lengtherror=lengtherror+1; %count
                    continue;
                end
                
                % channel check (making sure all channels are the same)
                
               % if strcmp(ch1,masterTableUse.senseSettings{f,1}.chan1{1,1}) && strcmp(ch2, masterTableUse.senseSettings{f,1}.chan2{1,1})
                    if strcmp(ch3,masterTableUse.senseSettings{f,1}.chan3{1,1}) && strcmp(ch4, masterTableUse.senseSettings{f,1}.chan4{1,1})
                        fprintf('\n channel check done');
                    else
                        fprintf('\n channel error');
                        channelerror = channelerror+1;
                        cerror{channelerror,1}=f;
                        continue;
                    end
%                         
%                 else
%                         fprintf('\t channel error \n');
%                         channelerror = channelerror+1;
%                         cerror{channelerror,1}=f;
%                         continue;
%                 end
                
                % Gather Stim on/off info
                

                onoff = masterTableUse.stimStatus{f,1}.stimulation_on;   
                
                if onoff == 1
                    StimSess{z,1} = FlName;
                    stimtype{z,1} =  masterTableUse.stimStatus{f,1}.amplitude_mA;
                    channelinfo{z,1} = masterTableUse.stimStatus{f,1}.electrodes;
                    z=z+1;
      
                    
                else
                    StimSess{x,2} = FlName;
                    x=x+1;
                end
                

                
            else
                fprintf('\n File %s from %s has a different SR',FlName,Groups{G});  % SR check
                diffSR = diffSR+1;
                continue;
            end
            

            %% Mean subtract to remove DC Offset
            
            temp0 = outdatcomplete.key0 - mean(outdatcomplete.key0);
            temp1 = outdatcomplete.key1 - mean(outdatcomplete.key1);
            temp2 = outdatcomplete.key2 - mean(outdatcomplete.key2);
            temp3 = outdatcomplete.key3 - mean(outdatcomplete.key3);
            
            
            %% Keeping a copy of the DC offset removed raw data for the next step
            
            Data.temp0 = temp0;
            Data.temp1 = temp1;
            Data.temp2 = temp2;
            Data.temp3 = temp3;

            %% Artifact Rejection minute by minute
            
           fprintf('\n Artifact Rejecting...');
           
           stdtemp0 = 3*std(temp0);
           stdtemp1 = 3*std(temp1);
           stdtemp2 = 3*std(temp2);
           stdtemp3 = 3*std(temp3);
           
           
           index0 = union(find(temp0>mean(temp0)+stdtemp0),find(temp0<mean(temp0)-(stdtemp0)));
           index1 = union(find(temp1>mean(temp1)+stdtemp1),find(temp1<mean(temp1)-(stdtemp1)));
           index2 = union(find(temp2>mean(temp2)+stdtemp2),find(temp2<mean(temp2)-(stdtemp2)));
           index3 = union(find(temp3>mean(temp3)+stdtemp3),find(temp3<mean(temp3)-(stdtemp3)));
           
           temp0(union(index0,union(index2,index3))) = NaN;
           temp1(union(index0,union(index2,index3))) = NaN;
           temp2(union(index0,union(index2,index3))) = NaN;
           temp3(union(index0,union(index2,index3))) = NaN;
            
            % Remove missing NaN values
            
            temp0 = rmmissing(temp0);
            temp1 = rmmissing(temp1);
            temp2 = rmmissing(temp2);
            temp3 = rmmissing(temp3);
            
            %% Filter Data
            
            fprintf('\n Highpass filtering...');
            
            Nf = fs/2; %nyquist frequency
            B = fir1(filterOrder,[High/Nf],'high');
            ftemp0 = filtfilt(B,1,temp0);
            ftemp1 = filtfilt(B,1,temp1);           
            ftemp2 = filtfilt(B,1,temp2);
            ftemp3 = filtfilt(B,1,temp3);
            
            %% Coherence minute by minute
            
            fprintf('\n Computing minute by minute Coherence...');

            s=1;
            
            for i=1:(floor(size(ftemp0,1)/(250*60)))
            
            l = 250*i*60;
            
            [Data.cxy(:,i),Data.fone(:,i)]=mscohere(ftemp0(s:l),ftemp2(s:l),fs,fs/2,[1:60],fs);
            [Data.cyz(:,i),Data.fone(:,i)]=mscohere(ftemp0(s:l),ftemp3(s:l),fs,fs/2,[1:60],fs);
            [Data.czx(:,i),Data.fone(:,i)]=mscohere(ftemp2(s:l),ftemp3(s:l),fs,fs/2,[1:60],fs);
            
            s = l+1;
            
            end
            
            clear i
                       
             %% PSDs of the 1 minute chunks
            
             fprintf('\n Computing minute by minute PSDs');

             s=1;
             for i=1:(floor(size(ftemp0,1)/(250*60)))
             l = 250*i*60;  
             [Data.psd0(:,i),Data.h1(:,i)] = pwelch(temp0(s:l),fs*5,fs,[1:60],fs);
             [Data.psd2(:,i),Data.h2(:,i)] = pwelch(temp2(s:l),fs*5,fs,[1:60],fs);
             [Data.psd3(:,i),Data.h3(:,i)] = pwelch(temp3(s:l),fs*5,fs,[1:60],fs);
             s = l+1;
             end

            if onoff==1
            TotalONminutes = i+TotalONminutes;
            else
            TotalOFFminutes = i+TotalOFFminutes;
            end
            
            clear i
            %% SAVING
                     
            filesprocessed = filesprocessed+1;    %file count
            fprintf('\n This file was processed fully');
            fprintf('\n SAVING');
            save(sprintf('%s/processed_data_final.mat',path),'Data','-v7.3')

            
             clear Data
             clear outdatcomplete
             clear timecheck1 timecheck2 Date Date2
        catch           
           fprintf('\n file throwing error is %s and findex is %d',FlName,f);
           throwerror{k,1} = f;
           k=k+1;
        end % main try loop
    end % f for loop   
    save(sprintf('%s/%s/StimSess_final.mat',InDir,Groups{G}),'stimtype','channelinfo','StimSess','TotalONminutes','TotalOFFminutes','throwingerror','lengtherror','diffSR','daterror','filesprocessed','sleeperror','channelerror','stimtype','-v7.3')
end % Groups


%% plot check sample code for psd & coherence

%plot(Data.h1(:,1),log10((Data.psd0(:,1)')))
%plot(Data.cxy(:,1))