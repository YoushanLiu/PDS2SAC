function pds2sac24

% pad2sac24 - convert files in 24-bit PDS format into SAC format
% It can cut data into segments
% Author: Youshan Liu
% Affiliation: Institute of Geology and Geophysics, Chinese Academy of Sciences
%
clear all;
clc;
fclose all;

% pool = gcp;
% if (~isempty(pool))
%     delete(pool);
% else
%     pool = parpool();
% end


downsampling_rate = 20.0;
Seconds_segment = 3600;


inpath = './JD.Group4.Raw';
outpath = './SAC.Group4';
inpath = './JD.Group5.PDS.Raw';
outpath = './SAC.Group5.PDS';
% inpath = './ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½';
% outpath = './Group2';
% inpath = './JD2';
% outpath = './Group2';


if (~exist(outpath, 'dir'))
    mkdir(outpath);
end

network = 'NCISP9';
stainfo = read_stainfo('STATIONS_JD.dat');


stage_path = [inpath, '/'];
stage_folders_list = dir([stage_path, 'DATA_*']);
nstage_folders_list = length(stage_folders_list);

for istage = 1:1:nstage_folders_list

    station_path = [stage_path, stage_folders_list(istage).name, '/'];

    %station_folders_list = dir([station_path, '????']);
    station_folders_list = dir([station_path, '5*']);
    nstation_folders_list = length(station_folders_list);

    parfor istation = 1:1:nstation_folders_list
%     for istation = 1:1:6

%         istation
        daily_path = [station_path, station_folders_list(istation).name, '/'];

        daily_files_list = dir([daily_path, '*_*']);
        nfiles = length(daily_files_list);

        for ifile = 1:nfiles

            filename = [daily_path, daily_files_list(ifile).name];
%             filename = './èƒ¶ä¸œï¿?ï¿½ï¿½/å±±ä¸œï¿?ï¿½ï¿½æ•°æ®/2017.3.24/1714251_.1140'

% filename = './ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½/ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?2017.03.22/1710100_.1268';
% filename = './èƒ¶ä¸œï¿?ï¿½ï¿½/å±±ä¸œï¿?ï¿½ï¿½æ•°æ®/2017.3.27/2009195_.1051'
%2017.03.22/1709452_.1266
%2017.03.22/1709562_.1265
% filename = './èƒ¶ä¸œå››ç»„/èƒ¶ä¸œå››ç»„æ•°æ®/2017.03.22/1710100_.1268'
% filename = './èƒ¶ä¸œå››ç»„/èƒ¶ä¸œå››ç»„æ•°æ®/2017.03.22/1709562_.1265'
% filename = './èƒ¶ä¸œå››ç»„/èƒ¶ä¸œå››ç»„æ•°æ®/2017.03.25/1810105_.1145'

            fprintf('\n\nConverting %s ...\n', filename);
            readpds(filename, network, stainfo, downsampling_rate, Seconds_segment, outpath);

        end
        
    end

end

% delete(pool);



function readpds(filename, network, stainfo, downsampling_rate, Seconds_segment, outpath)


fidin = fopen(filename, 'r', 'n');

fseek(fidin, 0, 'eof');
nbytes = ftell(fidin);
npackets = ceil(nbytes/512) - 3;

frewind(fidin);

% read control block
packet = fread(fidin, 512, 'uchar');


% get header terms
% pdshead.head = strcat(packet(1:4).');
% pdshead.version = strcat(packet(5:8).');
% pdshead.dasid = strcat(packet(36:-1:33).');
% pdshead.station = strcat(packet(40:-1:37).');
% pdshead.location = location_codes(packet(41));
% pdshead.filter_type = packet(43);
% pdshead.gain = 2^packet(45);
% pdshead.sampling_rate = 100*2^packet(47);
% pdshead.delta = 1.0 / pdshead.sampling_rate;
station = strcat(packet(40:-1:37).'); % This indicates big_endian
sps_type = packet(47);
sampling_rate = 100*2^sps_type;
decimate_rate = fix(sampling_rate/downsampling_rate);
dt = 1.0 / sampling_rate;
gain_inv = 1.0 / (2^packet(45));
filter_type = packet(43);
filter_delay = 0;
if (1 == filter_type)     % linear phase filter
    if (0 == sps_type)
        filter_delay = 0.2300;
    elseif(1 == sps_type)
        filter_delay = 0.1150;
    elseif(2 == sps_type)
        filter_delay = 0.0575;
    end
elseif(2 == filter_type)  % minimum phase filter
    if (0 == sps_type)
        filter_delay = 0.4500; % [0 - 40  Hz]
    elseif(1 == sps_type)
        filter_delay = 0.2250; % [0 - 80  Hz]
    elseif(2 == sps_type)
        filter_delay = 0.1125; % [0 - 160 Hz]
    end
    % correction is not applied
    % it may be exp(-1i*2*pi*f*tau(f)), where tau = max(linspace(0, fmax,
    % nf), fmax), fmax is the corresponding cut-off frequency for minimum
    % phase filter
end




% read first clock block
fseek(fidin, 512, 0);
% clock_block1 = fread(fidin, 512, 'uchar');
% clock_block1 = reshape(clock_block1, 8, 64);
% ptr_clock1 = clock_block1(1, 64);
% clock_block1 = zeros(8, 64);
% for k = 1:1:64
%     clock_block1(:,k) = fread(fidin, 8, 'uchar');
% end
% ptr_clock1 = clock_block1(1,64);
% packet = fread(fidin, 512, 'uchar');
% packet = reshape(packet, 8, 64);
% nclock_diff = hex2dec(packet(1,63));
% for k = 1:1:nclock_diff
%     GPS_time.day = packet(1,k);
%     GPS_time.hour = packet(2,k);
%     GPS_time.minute = packet(3,k);
%     GPS_time.second = packet(4,k);
%     DAS_time.minute = packet(5,k);
%     DAS_time.second = packet(6,k);
%     DAS_time.microsecond = (packet(7,k)*10 + packet(8,k));
% end



% read second clock block
fseek(fidin, 512, 0);
% clock_block2 = fread(fidin, 512, 'uchar');
% clock_block2 = reshape(clock_block2, 8, 64);
% ptr_clock2 = clock_block2(1, 64);
% clock_block2 = zeros(8, 64);
% for k = 1:1:64
%     clock_block2(:,k) = fread(fidin, 8, 'uchar');
% end
% ptr_clock2 = clock_block2(1,64);


% stafields = {'starttime'; 'stla'; 'stlo'; 'stel'; 'num_gps'; ...
%              'gps_status'; 'gps_lock_num'; 'gps_lock_time'};
% stahead = cell(size(stafields,1), npackets);
% stahead(1,:) = {' '};
% stahead(2:7,:) = {nan};
% stahead(8,:) = {' '};
% stahead = cell2struct(stahead, stafields, 1);





% minimum GPS stars
% if the number of GPS stars less than it, I consider its GPS is unlocked.
nGPS_min = 4;




offset = 2^24;


half_dt = 0.5*dt;
% Seconds_segment = 3600;
Seconds_packet = 49*dt;
% Seconds_half_segment = 0.5*Seconds_segment;
npts_segment = nearest(Seconds_segment/dt);


[b, a] = butter(2, 2*dt*0.499*downsampling_rate, 'low');
sacn = zeros(npts_segment, 1);
sace = zeros(npts_segment, 1);
sacz = zeros(npts_segment, 1);

% wfn = zeros(50,1);
% wfe = zeros(50,1);
% wfz = zeros(50,1);

starttime = datetime();





if (strcmp(filename, './ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½/ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½/2017403/2714492_.1131'))
    station = '2127';
end
% if (strcmp(filename, './JD2/JD2/2017403/2714492_.1131'))
%     station = '2127';
% end
if (strcmp(filename, './ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½/ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½/2017.4.05/2908084_.0586'))
    station = '3161';
end
% if (strcmp(filename, './JD3/JD3/2017.4.05/2908084_.0586'))
%     station = '3161';
% end
% ç¬¬äºŒï¿?2017403/2714492_.1131
if (~isempty(findstr(filename, '2714492_.1131')))
    station = '2127';
    fprintf('Error: head keyword station is wrongle set in filename %s\n', filename);
end
% ç¬¬ä¸‰ï¿?2017.4.05/2908084_.0586
if (~isempty(findstr(filename, '2908084_.0586')))
    station = '3161';
    fprintf('Error: head keyword station is wrongle set in filename %s\n', filename);
end

% station
if (strcmp('4339', station))
    station = '4239';
    fprintf('Error: station is wrongle set in filename %s\n', filename);
end
if (strcmp('0802', station))
    station = '4266';
    fprintf('Error: station is wrongle set in filename %s\n', filename);
end


path_splitted = regexp(filename, '[\\\/]', 'split');
path_splitted = regexp(path_splitted{end-1}, '_', 'split');
station = path_splitted{1};


% filter_delay
% station
% return





% %ista = str2double(station(2:end));
% ista = strmatch(station, stainfo.stnm);
% stla = stainfo.stla(ista);
% stlo = stainfo.stlo(ista);
% stel = stainfo.stel(ista);

ista = strmatch(station, stainfo.stnm);
% ista = find(1 == contains(stainfo.stnm, station));
% ista = find(1 == startsWith(stainfo.stnm, station));
stla = str2double(cell2mat(stainfo.stla(ista)));
stlo = str2double(cell2mat(stainfo.stlo(ista)));
stel = str2double(cell2mat(stainfo.stel(ista)));


% station
% stla
% return


npts = 0;
nskip = 0;
for k = 1:1:npackets

    wfn = zeros(50,1);
    wfe = zeros(50,1);
    wfz = zeros(50,1);

    packet = fread(fidin, 500, 'uchar');
    packet = reshape(packet, 10, 50);
    fseek(fidin, 12, 0);


    % datetime
    Year = packet(10,1);
    if (Year < 70)
        Year = Year + 2000;
    else
        Year = Year + 1900;
    end
    Month = packet(10,2);
    Day = packet(10,3);
    Hour = packet(10,4);
    % convert to UTC time
    Hour = Hour - 8;
    Minute = packet(10,5);
    Second = packet(10,6);
    MicroSecond = packet(10,7)*10 + packet(10,8);


    starttime_prev = starttime;
    starttime = datetime(Year, Month, Day, Hour, Minute, Second, MicroSecond, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');
    % endtime = datetime(Year, Month, Day, Hour, Minute, Second, MicroSecond + 1e3*(49*dt), 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');
    if (1 == filter_type)
        starttime = starttime + seconds(filter_delay);
    end
    % endtime = starttime;
    %endtime.Second = endtime.Second + 49*dt;
    endtime = starttime + seconds(Seconds_packet);



    % if the number of GPS in current packet is less than the nGPS_min
    % (i.e. 4), I consider this packet is GPS unlocked. Then I discard it.
    if (packet(10,33) < nGPS_min)
        if (npts > 0)
            % write into a new file because of discontinuous clock
            % write sac file
            save_sacfile(sacz, sacn, sace, dt, npts, decimate_rate, b, a, network, ...
                    station, stla, stlo, stel, starttime_segment, outpath, nskip);
        end
        npts = 0;
        nskip = 0;
        continue;
    end

 


    is_continuous = true;
    if (k > 1)
        % if (abs(etime(datevec(starttime), datevec(starttime_prev))) > 50.00001*dt)
        if (round(etime(datevec(starttime), datevec(starttime_prev))/dt) ~= 50)
            is_continuous = false;
        end
    end



    if (0 == npts)
        %starttime_segment = starttime;
        %starttime_segment.Second = starttime_segment.Second + nskip*dt;

        % midtime_segment = starttime_segment;
        % midtime_segment.Second = midtime_segment.Second + Seconds_half_segment;
        % endtime_segment = datetime(midtime_segment.Year, midtime_segment.Month, midtime_segment.Day, midtime_segment.Hour, ...
        %                                                    0, Seconds_segment, 0, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');
        % endtime_segment.Second = endtime_segment.Second - dt;

        % midtime_segment = starttime_segment;
        % midtime_segment.Second = midtime_segment.Second + Seconds_half_segment;
        % endtime_segment = datetime(midtime_segment.Year, midtime_segment.Month, midtime_segment.Day, midtime_segment.Hour+1, ...
        %                                                                  0, 0, 0, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');
        %timestr = sprintf('%4.4d-%3.3dT%2.2d:%2.2d:%2.2d.%3.3d', starttime_segment.Year, dayofyear(starttime_segment)+1, 0, 0, 0, 0);
        %endtime_segment = datetime(timestr, 'format', 'uuuu-DDD''T''HH:mm:ss.SSS');
        % endtime_segment.Second = endtime_segment.Second - dt;

        %endtime_segment = starttime_segment + hours(1);
        %endtime_segment.Minute = 59; endtime_segment.Second = 60 - 0.5*dt;
 
        starttime_segment = starttime + seconds(nskip*dt);
        endtime_segment = datetime(starttime_segment.Year, starttime_segment.Month, starttime_segment.Day, starttime_segment.Hour, ...
                                                                            59, 59, 999, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');
    end
% starttime_segment
% endtime_segment
% return



    % get eatra header information
    % latitude
    % stla_key = char(packet(10,9));
    % stla = (str2double(char(packet(10,10)))*10 + str2double(char(packet(10,11)))) + ...
    %       ((str2double(char(packet(10,12)))*10 + str2double(char(packet(10,13)))) + ...
    %        1.e-4*(str2double(char(packet(10,14)))*1e3 + str2double(char(packet(10,15)))*1e2 + ...
    %        str2double(char(packet(10,16)))*1e1 + str2double(char(packet(10,17))))) / 60.0;
    % if (('S' == stla_key) || ('s' == stla_key))
    %     stla = -stla;
    % end
    % stahead(k).stla = stla;

    % longitude
    % stlo_key = char(packet(10,18));
    % stlo = (str2double(char(packet(10,19)))*100 + str2double(char(packet(10,20)))*10 + str2double(char(packet(10,21)))) + ...
    %       ((str2double(char(packet(10,22)))*10 + str2double(char(packet(10,23)))) + ...
    %        1.e-4*(str2double(char(packet(10,24)))*1e3 + str2double(char(packet(10,25)))*1e2 + ...
    %        str2double(char(packet(10,26)))*1e1 + str2double(char(packet(10,27))))) / 60.0;
    % if (('W' == stlo_key) || ('w' == stlo_key))
    %    stlo = -stlo;
    % end
    % stahead(k).stlo = stlo;

    % elevtion
    % stel = str2double(strcat([packet(10,28), packet(10,29), packet(10,30), packet(10,31), packet(10,32)]));
    % stahead(k).stel = stel;
% stla
% stlo
% stel


    % %ista = str2double(station(2:end));
    % ista = strmatch(station, stainfo.stnm);
    % stla = stainfo.stla(ista);
    % stlo = stainfo.stlo(ista);
    % stel = stainfo.stel(ista);



    % GPS
    % number of GPS
    % stahead(k).num_gps = packet(10,33);
    % GPS status
    % stahead(k).gps_status = packet(10,34);
    % number of GPS locked
    % stahead(k).gps_lock_num = packet(10,35);
    % time of last locked GPS
    % Day = packet(10,39);
    % Hour = packet(10,40);
    % convert to UTC time
    % Hour = Hour - 8;
    % Minute = packet(10,41);
    % Second = packet(10,42);
    % MicroSecond = packet(10,43)*10 + packet(10,44);
    % stahead(k).gps_lock_time = datetime(Year, Month, Day, Hour, Minute, Second, MicroSecond, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');



    % convert data
    wfn(:) = (packet(1,:)*256 + packet(2,:))*256 + packet(3,:);
    indx = find(packet(1,:) >= 128);
    wfn(indx) = wfn(indx) - offset;

    wfe(:) = (packet(4,:)*256 + packet(5,:))*256 + packet(6,:);
    indx = find(packet(4,:) >= 128);
    wfe(indx) = wfe(indx) - offset;

    wfz(:) = (packet(7,:)*256 + packet(8,:))*256 + packet(9,:);
    indx = find(packet(7,:) >= 128);
    wfz(indx) = wfz(indx) - offset;
    
    % minus is introduced by the south-positive
    wfn = -wfn*gain_inv;
    wfe = +wfe*gain_inv;
    wfz = +wfz*gain_inv;



% if(k >= 4166)
% k
% starttime
% endtime
% starttime_segment
% endtime_segment
% end


    if (~is_continuous)
        % write into a new file because of discontinuous clock
        % write sac file
        save_sacfile(sacz, sacn, sace, dt, npts, decimate_rate, b, a, network, ...
                station, stla, stlo, stel, starttime_segment, outpath, nskip);

        npts = 0;
        nskip = 0;

        % starttime_segment = starttime;

        % midtime_segment = starttime_segment;
        % midtime_segment.Second = midtime_segment.Second + Seconds_half_segment;
        % endtime_segment = datetime(midtime_segment.Year, midtime_segment.Month, midtime_segment.Day, midtime_segment.Hour, ...
        %                                                   0, Seconds_segment, 0, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');
        % endtime_segment.Second = endtime_segment.Second - dt;

        % midtime_segment = starttime_segment;
        % midtime_segment.Second = midtime_segment.Second + Seconds_half_segment;
        % endtime_segment = datetime(midtime_segment.Year, midtime_segment.Month, midtime_segment.Day, midtime_segment.Hour+1, ...
        %                                                                  0, 0, 0, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');
        %timestr = sprintf('%4.4d-%3.3dT%2.2d:%2.2d:%2.2d.%3.3d', starttime_segment.Year, dayofyear(starttime_segment)+1, 0, 0, 0, 0);
        %endtime_segment = datetime(timestr, 'format', 'uuuu-DDD''T''HH:mm:ss.SSS');
        % endtime_segment.Second = endtime_segment.Second - dt;

        %endtime_segment = starttime_segment + hours(1);
        %endtime_segment.Minute = 59; endtime_segment.Second = 60 - 0.5*dt;


        starttime_segment = starttime;
        %endtime_segment = starttime_segment + seconds(Seconds_segment - half_dt);
        endtime_segment = datetime(starttime_segment.Year, starttime_segment.Month, starttime_segment.Day, starttime_segment.Hour, ...
                                                                            59, 59, 999, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');

        % cat waveform into sac arrays
        sacz(npts+1:npts+50) = wfz(1:50);
        sacn(npts+1:npts+50) = wfn(1:50);
        sace(npts+1:npts+50) = wfe(1:50);
        npts = npts + 50;
        continue;
    end



    if (is_continuous)

        if (endtime_segment >= endtime)
            %
            % starttime_segment                       endtime_segment
            %         |                                      |
            %         v                                      v                      
            %         |______________________________________|
            %                                  |_________|
            %                              starttime  endtime
            %
            % cat waveform into sac arrays
            sacz(npts+1:npts+50) = wfz(1:50);
            sacn(npts+1:npts+50) = wfn(1:50);
            sace(npts+1:npts+50) = wfe(1:50);
            npts = npts + 50;
            if (npackets == k)
                % At the end of the file, save the last segment into files
                % write sac files
                save_sacfile(sacz, sacn, sace, dt, npts, decimate_rate, b, a, network, ...
                        station, stla, stlo, stel, starttime_segment, outpath, nskip);
                break;
            end
        else
            %
            % starttime_segment                       endtime_segment
            %         |                                      |
            %         v                                      v                      
            %         |______________________________________|
            %                                           |_________|
            %                                       starttime  endtime
            %
            % extract waveform into sac arrays
            time_diff = etime(datevec(endtime_segment), datevec(starttime));

            %if (time_diff > 0.5*dt)
            if (time_diff >= 0.0)
                %
                % starttime_segment                       endtime_segment
                %         |                                      |
                %         v                                      v                      
                %         |______________________________________|
                %                                    <------     |_________|
                %                                            starttime  endtime
                %
                % npts_part = min(nearest(time_diff/dt), 50);
                npts_part = min(ceil(time_diff/dt), 50);
                sacz(npts+1:npts+npts_part) = wfz(1:npts_part);
                sacn(npts+1:npts+npts_part) = wfn(1:npts_part);
                sace(npts+1:npts+npts_part) = wfe(1:npts_part);
                npts = npts + npts_part;

% npts
% k
% starttime
% starttime_segment
% endtime
% endtime_segment
% etime(datevec(endtime_segment), datevec(starttime))
% starttime_segment
% nskip
% npts_out = length(nskip+1:5:npts)




% k
% starttime_segment
% endtime_segment

                % write sac files
                nskip = save_sacfile(sacz, sacn, sace, dt, npts, decimate_rate, b, a, network, ...
                                station, stla, stlo, stel, starttime_segment, outpath, nskip);



% starttime_segment
% time0 = starttime_segment;
% starttime_segment = starttime;
% starttime_segment.Second = starttime_segment.Second + (npts_part + nskip)*dt
% 
% nskip
% npts_part
% npts
% 
% time1 = time0;
% time1.Second = time0.Second + (npts_out-1)*0.05
% time2 = time0;
% time2.Second = time0.Second + (371531-1)*dt
% 
% time4 = time1;
% time4.Second = time4.Second + 0.05
% % time0.Second + (npts_out-1)*0.05
% % time0.Second + (371531-1)*dt
% % (npts_out-1)*0.05
% % (371531-1)*dt
% 
% % time3 = time0;
% % time3.Second = time3.Second + 3.7153e+03
% % starttime
% wfz(npts_part+4)
% k
% % sacz(1)
% % return
% if (14631 == k)
%     return
% end



                npts = 0;

                % starttime_segment = starttime;
                % starttime_segment.Second = starttime_segment.Second + (npts_part - 1 + nskip + 1)*dt;
                % starttime_segment.Second = starttime_segment.Second + (npts_part + nskip)*dt;
                % midtime_segment = starttime_segment;
                % midtime_segment.Second = midtime_segment.Second + Seconds_half_segment;
                % endtime_segment = datetime(midtime_segment.Year, midtime_segment.Month, midtime_segment.Day, midtime_segment.Hour, ...
                %                                                   0, Seconds_segment, 0, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');

                %timestr = sprintf('%4.4d-%3.3dT%2.2d:%2.2d:%2.2d.%3.3d', starttime_segment.Year, dayofyear(starttime_segment)+1, 0, 0, 0, 0);
                %endtime_segment = datetime(timestr, 'format', 'uuuu-DDD''T''HH:mm:ss.SSS');
                % endtime_segment.Second = endtime_segment.Second - dt;

                %endtime_segment = starttime_segment + hours(1);
                %endtime_segment.Minute = 59; endtime_segment.Second = 60 - 0.5*dt;

                % midtime_segment = starttime_segment;
                % midtime_segment.Second = midtime_segment.Second + Seconds_half_segment;
                % endtime_segment = datetime(midtime_segment.Year, midtime_segment.Month, midtime_segment.Day, midtime_segment.Hour+1, ...
                %                                                                   0, 0, 0, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');
                % endtime_segment.Second = endtime_segment.Second - dt;


                starttime_segment = starttime + seconds((npts_part + nskip)*dt);
                % starttime_segment = starttime + seconds((npts_part - 1 + nskip + 1)*dt);
                endtime_segment = starttime_segment + seconds(Seconds_segment - half_dt);


                n = 50 - npts_part;

% starttime_segment
% endtime_segment
% if (strcmp('1610560_.1172', filename(end-12:end)))
%     dt
%     endtime_segment
%     starttime
%     endtime
%     etime(datevec(endtime_segment), datevec(starttime))/dt
%     nearest(etime(datevec(endtime_segment), datevec(starttime))/dt)
%     npts_part
%     n
% end
% n
% npts
% npts_part
                sacz(npts+1:npts+n) = wfz(npts_part+1:50);
                sacn(npts+1:npts+n) = wfn(npts_part+1:50);
                sace(npts+1:npts+n) = wfe(npts_part+1:50);
                npts = npts + n;
            else
                %
                % starttime_segment                      endtime_segment
                %         |                                      |
                %         v                                      v                      
                %         |______________________________________|
                %                                                |_________|  ------>
                %                                            starttime  endtime
                %
                % write sac files
                nskip = save_sacfile(sacz, sacn, sace, dt, npts, decimate_rate, b, a, network, ...
                                station, stla, stlo, stel, starttime_segment, outpath, nskip);


                npts = 0;

                % starttime_segment = starttime;

                % starttime_segment.Second = starttime_segment.Second + (nskip + 1)*dt;
                % midtime_segment = starttime_segment;
                % midtime_segment.Second = midtime_segment.Second + Seconds_half_segment;
                % endtime_segment = datetime(midtime_segment.Year, midtime_segment.Month, midtime_segment.Day, midtime_segment.Hour, ...
                %                                                   0, Seconds_segment, 0, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');

                %timestr = sprintf('%4.4d-%3.3dT%2.2d:%2.2d:%2.2d.%3.3d', starttime_segment.Year, dayofyear(starttime_segment)+1, 0, 0, 0, 0);
                %endtime_segment = datetime(timestr, 'format', 'uuuu-DDD''T''HH:mm:ss.SSS');
                % endtime_segment.Second = endtime_segment.Second - dt;
                %endtime_segment = starttime_segment + hours(1);
                %endtime_segment.Minute = 59; endtime_segment.Second = 60 - 0.5*dt;

                % midtime_segment = starttime_segment;
                % midtime_segment.Second = midtime_segment.Second + Seconds_half_segment;
                % endtime_segment = datetime(midtime_segment.Year, midtime_segment.Month, midtime_segment.Day, midtime_segment.Hour+1, ...
                %                                                                   0, 0, 0, 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSS');
                % endtime_segment.Second = endtime_segment.Second - dt;


                starttime_segment = starttime + seconds((nskip+1)*dt);
                endtime_segment = starttime_segment + seconds(Seconds_segment - half_dt);


                % cat waveform into sac arrays
                sacz(npts+1:npts+50) = wfz(1:50);
                sacn(npts+1:npts+50) = wfn(1:50);
                sace(npts+1:npts+50) = wfe(1:50);
                npts = npts + 50;
            end % end of if (time_diff > 0.5*dt)

        end % end of if (endtime_segment >= endtime)

    end % end of if (is_continuous)


% if (3502+7200*2 == k)
% % if (3502 == k)
%     k
%     break
% end


end


clear b a;
clear packet;
clear wfn wfe wfz;
clear sacn sace sacz;


fclose(fidin);



function nskip = save_sacfile(sacz, sacn, sace, dt, npts, decimate_rate, b, a, network, ...
                         station, stla, stlo, stel, starttime_segment, outpath, nskip)


% preprocessing
% downsampling
[sacz_out, sacn_out, sace_out, dt_out, npts_out, nskip] = downsampling(sacz, sacn, sace, dt, npts, ...
                                                                      decimate_rate, b, a, nskip);
% end of preprocessing


% t1 = (0:1:npts-1)'*dt;
% t2 = (0:1:npts_out-1)'*dt_out;
% 
% dt
% dt_out
% 
% max(t1)
% max(t2)
% 
% figure(4);
% subplot(3,1,1);
% hold off;
% % plot(wf(:,1), 'r');
% plot(t1, sacz(1:npts), 'r');
% hold on;
% plot(t2, sacz_out(1:npts_out), '--b');
% hold off;
% subplot(3,1,2);
% hold off;
% % plot(wf(:,2), 'r');
% plot(t1, sacn(1:npts), 'r');
% hold on;
% plot(t2, sacn_out(1:npts_out), '--b');
% hold off;
% subplot(3,1,3);
% hold off;
% % plot(wf(:,3), 'r');
% plot(t1, sace(1:npts), 'r');
% hold on;
% plot(t2, sace_out(1:npts_out), '--b');
% hold off;
% return


% remove mean
sacn_out = detrend(sacn_out, 'constant');
sace_out = detrend(sace_out, 'constant');
sacz_out = detrend(sacz_out, 'constant');
% remove linear trend
sacn_out = detrend(sacn_out, 'linear');
sace_out = detrend(sace_out, 'linear');
sacz_out = detrend(sacz_out, 'linear');


%sacn_out = sacn_out - mean(sacn_out);
%sace_out = sace_out - mean(sace_out);
%sacz_out = sacz_out - mean(sacz_out);
%if (npts_out > 6)
%    sacn_out = detrend(sacn_out);
%    sace_out = detrend(sace_out);
%    sacz_out = detrend(sacz_out);
%end


Year = starttime_segment.Year;
Hour = starttime_segment.Hour;
Minute = starttime_segment.Minute;
julday = day(starttime_segment, 'dayofyear');
% create output path
datestr = sprintf('%4.4d%3.3d', Year, julday);
station_daily_path = [outpath, '/', station, '/', datestr, '/'];
if (~exist(station_daily_path, 'dir'))
    mkdir(station_daily_path);
end


sec = starttime_segment.Second;
nzsec = fix(sec);
nzmsec = fix((sec - nzsec)*1000);
datestr = sprintf('%4.4d.%3.3d', Year, julday);
timestr = sprintf('%2.2d.%2.2d.%2.2d.%3.3d', Hour, Minute, nzsec, nzmsec);




prefix = [station_daily_path, datestr, '.', timestr, '.', network, '.', station, '..'];
% output sac file
cmp = 'BHN';
%outfile = [prefix, 'BHN.SAC'];
SAC = initi_sacheader([prefix, cmp, '.SAC'], starttime_segment, julday, nzsec, nzmsec, ...
                           stla, stlo, stel, network, station, npts_out, dt_out, cmp);
SAC.DATA1 = sacn_out;
writesac(SAC);
clear SAC;

cmp = 'BHE';
%outfile = [prefix, 'BHE.SAC'];
SAC = initi_sacheader([prefix, cmp, '.SAC'], starttime_segment, julday, nzsec, nzmsec, ...
                            stla, stlo, stel, network, station, npts_out, dt_out, cmp);
SAC.DATA1 = sace_out;
writesac(SAC);
clear SAC;

cmp = 'BHZ';
%outfile = [prefix, 'BHZ.SAC'];
SAC = initi_sacheader([prefix, cmp, '.SAC'], starttime_segment, julday, nzsec, nzmsec, ...
                           stla, stlo, stel, network, station, npts_out, dt_out, cmp);
SAC.DATA1 = sacz_out;
writesac(SAC);
clear SAC;


fprintf('%s is done ... \n', [datestr, '.', timestr, '.', network, '.', station, '..', 'BH*.SAC']);




function [sacz_out, sacn_out, sace_out, dt_out, npts_out, nskip] = downsampling(sacz, sacn, sace, dt, npts, ...
                                                                                decimate_rate, b, a, nskip)


sacz_segment = sacz(1:npts);
sacn_segment = sacn(1:npts);
sace_segment = sace(1:npts);
if ((1 ~= decimate_rate) && (npts > 12))
	% downsampling
	sacz_segment = filtfilt(b, a, sacz_segment);
	sacn_segment = filtfilt(b, a, sacn_segment);
	sace_segment = filtfilt(b, a, sace_segment);
end
indx = (nskip+1:decimate_rate:npts);
sacz_out = sacz_segment(indx);
sacn_out = sacn_segment(indx);
sace_out = sace_segment(indx);

npts_out = length(indx);
dt_out = decimate_rate*dt;

clear indx sacz_segment sacn_segment sace_segment;


% t1 = nskip*dt + (0:1:npts-1)'*dt;
% t2 = nskip*dt + (0:1:npts_out-1)'*dt_out;
% figure(2);
% n = 501;
% m = length((nskip+1:decimate_rate:n));
% hold off;
% % plot((nskip+1:n), sacz(nskip+1:n), 'r');
% % hold on;
% % plot((nskip+1:decimate_rate:n), sacz_out(1:m), '--b');
% plot((nskip+1:npts), sacz(nskip+1:npts), 'r');
% hold on;
% plot((nskip+1:decimate_rate:npts), sacz_out(1:npts_out), '--b');
% hold off;
% 
% 
% dt
% dt_out
% npts
% npts_out
% decimate_rate
% 
% max(t1)
% max(t2)
% 
% figure(1);
% subplot(3,1,1);
% hold off;
% % plot(wf(:,1), 'r');
% plot(t1, sacz(1:npts), 'r');
% hold on;
% plot(t2, sacz_out(1:npts_out), '--b');
% hold off;
% subplot(3,1,2);
% hold off;
% % plot(wf(:,2), 'r');
% plot(t1, sacn(1:npts), 'r');
% hold on;
% plot(t2, sacn_out(1:npts_out), '--b');
% hold off;
% subplot(3,1,3);
% hold off;
% % plot(wf(:,3), 'r');
% plot(t1, sace(1:npts), 'r');
% hold on;
% plot(t2, sace_out(1:npts_out), '--b');
% hold off;
% % return


% nskip = decimate_rate - (npts - ((nskip+1) + (npts_out-1)*decimate_rate)) - 1;
nskip = nskip + decimate_rate*npts_out - npts;




function SAC = initi_sacheader(outfile, starttime, julday, nzsec, nzmsec, stla, stlo, stel, network, station, npts, dt, cmp)

SAC = sacstruct(1);
SAC.FILENAME = outfile;
SAC.DELTA = dt;
SAC.SCALE = 1;
SAC.NPTS = npts;
SAC.B = 0.0;
SAC.O = 0.0;
SAC.E = (npts-1)*dt;
SAC.NZYEAR = starttime.Year;
%SAC.NZJDAY = day(starttime, 'dayofyear');
SAC.NZJDAY = julday;
SAC.NZHOUR = starttime.Hour;
SAC.NZMIN = starttime.Minute;
%SAC.NZSEC = fix(starttime.Second);
%SAC.NZMSEC = fix((starttime.Second - fix(starttime.Second))*1000);
SAC.NZSEC = nzsec;
SAC.NZMSEC = nzmsec;
SAC.NVHDR = 6;
SAC.STLA = stla;
SAC.STLO = stlo;
SAC.STEL = stel;
% SAC.MAG = 0;
SAC.EVLA = 0.0;
SAC.EVLO = 0.0;
SAC.EVDP = 0.0;
% SAC.IFTYPE = 1;
% SAC.IZTYPE = 9;
SAC.IFTYPE = 'ITIME';
SAC.IZTYPE = 'IB';
SAC.LEVEN = 1;
SAC.KSTNM = station;
SAC.KNETWK = network;
CMP = cmp(end:end);
if (('Z' == CMP) || ('z' == CMP))         % Z
    SAC.CMPAZ  = 0.0;
    SAC.CMPINC = 0.0;
    SAC.KCMPNM = 'BHZ';
elseif (('N' == CMP) || ('n' == CMP))     % N
    SAC.CMPAZ  = 0.0;
    SAC.CMPINC = 90.0;
    SAC.KCMPNM = 'BHN';
elseif (('E' == CMP) || ('e' == CMP))     % E
    SAC.CMPAZ  = 90.0;
    SAC.CMPINC = 90.0;
    SAC.KCMPNM = 'BHE';
end



function sta = read_stainfo(infile)

fid = fopen(infile, 'r');
    dat = [];
    while(1)
       line = fgetl(fid);
       if (-1 == line)
           break;
       end
       line = strtrim(line);
       line = regexp(line, '[\t ]*', 'split');
       %dat = [dat; [line{1}, line{2}, line{3}, line{4}]]
       dat = [dat; line];
   end
fclose(fid);
dat = dat';


stafields = {'stnm'; 'stla'; 'stlo'; 'stel'};
sta = cell(size(stafields,1), 1);
sta(1,:) = {dat(1,:)};
sta(2,:) = {dat(2,:)};
sta(3,:) = {dat(3,:)};
sta(4,:) = {dat(4,:)};
sta = cell2struct(sta, stafields, 1);

% dat = load(infile, '-ascii');
% 
% stafields = {'stnm'; 'stla'; 'stlo'; 'stel'};
% sta = cell(size(stafields,1), 1);
% sta(1,:) = {num2str(dat(:,1))};
% sta(2,:) = {dat(:,2)};
% sta(3,:) = {dat(:,3)};
% sta(4,:) = {dat(:,4)};
% sta = cell2struct(sta, stafields, 1);



function s = sacstruct(n)
% S = SACSTRUCT(N); returns an Nx1 SAC structure with all fields
% initiated to their undefined (default) valueSAC. N defaults to 1. 
%

%	Xiaoning Yang 2008 (modified from subfunction init_sac of readsac.m
%   in the MATSEIS package)

% check argumentSAC.
if (nargin < 1)
    n = 1;
end

% initialize output structure.
sacfields = {'FILENAME'; 'DELTA'; 'DEPMIN'; 'DEPMAX'; 'SCALE'; ...
    'ODELTA'; 'B'; 'E'; 'O'; 'A'; 'INTERNAL1'; 'T0'; 'T1'; 'T2'; 'T3'; ...
    'T4'; 'T5'; 'T6'; 'T7'; 'T8'; 'T9'; 'F'; 'RESP0'; 'RESP1'; 'RESP2'; ...
    'RESP3'; 'RESP4'; 'RESP5'; 'RESP6'; 'RESP7'; 'RESP8'; 'RESP9'; ...
    'STLA'; 'STLO'; 'STEL'; 'STDP'; 'EVLA'; 'EVLO'; 'EVEL'; 'EVDP'; ...
    'MAG'; 'USER0'; 'USER1'; 'USER2'; 'USER3'; 'USER4'; 'USER5'; ...
    'USER6'; 'USER7'; 'USER8'; 'USER9'; 'DIST'; 'AZ'; 'BAZ'; 'GCARC'; ...
    'INTERNAL2'; 'INTERNAL3'; 'DEPMEN'; 'CMPAZ'; 'CMPINC'; 'XMINIMUM'; ...
    'XMAXIMUM'; 'YMINIMUM'; 'YMAXIMUM'; 'UNUSED1'; 'UNUSED2'; ...
    'UNUSED3'; 'UNUSED4'; 'UNUSED5'; 'UNUSED6'; 'UNUSED7'; 'NZYear'; ...
    'NZJDay'; 'NZHour'; 'NZMIN'; 'NZSEC'; 'NZMSEC'; 'NVHDR'; 'NORID'; ...
    'NEVID'; 'NPTS'; 'INTERNAL4'; 'NWFID'; 'NXSIZE'; 'NYSIZE'; ...
    'UNUSED8'; 'IFTYPE'; 'IDEP'; 'IZTYPE'; 'UNUSED9'; 'IINST'; ...
    'ISTREG'; 'IEVREG'; 'IEVTYP'; 'IQUAL'; 'ISYNTH'; 'IMAGTYP'; ...
    'IMAGSRC'; 'UNUSED10'; 'UNUSED11'; 'UNUSED12'; 'UNUSED13'; ...
    'UNUSED14'; 'UNUSED15'; 'UNUSED16'; 'UNUSED17'; 'LEVEN'; 'LPSPOL'; ...
    'LOVROK'; 'LCALDA'; 'UNUSED18'; 'KSTNM'; 'KEVNM'; 'KHOLE'; 'KO'; ...
    'KA'; 'KT0'; 'KT1'; 'KT2'; 'KT3'; 'KT4'; 'KT5'; 'KT6'; 'KT7'; ...
    'KT8'; 'KT9'; 'KF'; 'KUSER0'; 'KUSER1'; 'KUSER2'; 'KCMPNM'; ...
    'KNETWK'; 'KDATRD'; 'KINST'; 'DATA1';};

cl = cell(size(sacfields,1),n);
cl(2:111, :) = {-12345};
cl(112:end-1, :) = {' '};
s = cell2struct(cl, sacfields, 1);
