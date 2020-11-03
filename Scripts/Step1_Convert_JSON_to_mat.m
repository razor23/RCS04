%% Convert JSON to .MAT files

clear all; close all; clc

% Load up the data folder with all the files in %

InDir = '/Users/laptop-0008/Desktop/Starr/Dystonia/RawData';
Groups = {'RCS04L','RCS04R'};
Devices = {'DeviceNPC700418H','DeviceNPC700412H'};

for G=1:length(Groups)
    Files = dir(sprintf('%s/%s/Session*',InDir,Groups{G}));
    js=1;
    for f = 1:length(Files)
        FlName = Files(f).name;
        path = sprintf('%s/%s/%s/%s',InDir,Groups{G},FlName,Devices{G});
        folderfiles = dir(path);
        index = find(ismember({folderfiles.name},'RawDataTD.mat'));
        
  boxdir = '/Users/laptop-0008/Desktop/Starr/Dystonia/RawData/RCS04L'  
  create_database_from_device_settings_files(boxdir);
       
  try
        % create patient database
    load(fullfile(boxdir,'database_from_device_settings.mat'));
    idxTable = cellfun(@(x) istable(x), masterTableOut.stimState);
    idxNonZero = masterTableOut.duration > seconds(5);
    idxUse = idxTable & idxNonZero;
    masterTableUse = masterTableOut(idxUse,:);
    allDeviceSettingsOut = allDeviceSettingsOut(idxUse);
    for s = 1:size(masterTableUse,1)
        [pn,fn] = fileparts(allDeviceSettingsOut{s});
        timeStart = report_start_end_time_td_file_rcs(fullfile(pn,'RawDataTD.json'));
        isValidTime = ~isempty(timeStart.duration); 
        if isValidTime
            timeStart.startTime.TimeZone             = 'America/Los_Angeles';
            timeStart.endTime.TimeZone               = 'America/Los_Angeles';
            masterTableUse.idxkeep(s) = 1;
            masterTableUse.timeStart(s) = timeStart.startTime;
            masterTableUse.timeEnd(s) = timeStart.endTime;
            masterTableUse.duration(s) = timeStart.duration;
        else
            masterTableUse.idxkeep(s) = 0;
            masterTableUse.timeStart(s) = NaT;
            masterTableUse.timeEnd(s) = NaT;
            masterTableUse.duration(s) = seconds(0);
        end
    end
    
    masterTableUse = masterTableUse(logical(masterTableUse.idxkeep),:);
    masterTableUse.duration.Format = 'hh:mm:ss';
        
        
       catch
           fprinft ('that didn work');
    
       end
       
       
        if index
            continue;
        else
            
            try
                fprintf('JSON file %s from %s is missing \n',FlName,Groups{G});
                Infile = sprintf('%s/%s/%s/%s/RawDataTD.json',InDir,Groups{G},FlName,Devices{G});
                i=i+1;
                jsonobj = deserializeJSON(Infile);
                [outdat, srates] = unravelData(jsonobj);
                outdatcomplete = outdat;
                SR = getSampleRate(srates);
                unqsrates = unique(SR);
                save(sprintf('%s/RawDataTD.mat',path),'outdatcomplete','unqsrates');
            catch
                fprintf('Loading file %s from %s is throwing an error \n',FlName,Groups{G});
                JSONERRORS{js,G} = FlName;
                js= js+1;
            end
        end % conversion if
    end % f for loop
    outpath = sprintf('%s/%s/%s/%s',InDir,Groups{G});
    save(sprintf('%s/Errorfiles.mat',outpath),'JSONERRORS');
end % Groups