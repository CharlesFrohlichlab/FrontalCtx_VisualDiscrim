%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%
% SpikeSort.m
%
%
%
% Sorts spikes from MEA and MCrack. Based on previous code developed over the years by
% Flavio Frohlich. Spikes realigned based on negative peak.
%d
% (c) Frohlich Lab 2012, 2013
%
% Author: Flavio Frohlich, Davis Bennett, Steve Schmidt
%
% History
%
% 1/21/2014 (ZZ) Cleaned up code to accomodate data mule with struct data
% 11/04/2013(SS) Updated to work with channel numbers as they appear for
% the ketamine project (11 to 106)
% 09/25/2013(SS) Updated to work with new spike/trial extraction script
% 11/15/2012(SS) Added Realign code from Davis
% 10/18/2012(SS) Updated code to only prompt sorting for spikes in area
%                of interest as determined by ElectrodeLayers2.csv
% 7/26/2012 (FF) User selects channels, checks if already sorted, nicer
%                plots and other minor upgrades
% 7/25/2012 (FF) Plot all raw traces, remove excess samples before neg
%                phase
% 7/24/2012 (FF) Initial version
%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% Begin: Prepare workspace
clear all % delete all variables
close all % close all windows
clc % clear command window
tic % start clock

if matlabpool('size') > 0
    matlabpool close
end

%matlabpool open 4;
% End: Prepare workspace
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$





% TEMP


%used for mapping and determining eletrodes in area
%channelMap = [47,48,46,45,38,37,28,36,27,17,26,16,35,25,15,14,24,34,13,23,12,22,33,21,32,31,44,43,41,42,52,51,53,54,61,62,71,63,72,82,73,83,64,74,84,85,75,65,86,76,87,77,66,78,67,68,55,56,58,57;];
channelMap = [1:16];
% Plot parameters
perfig = 8; % Number of panels per figure
my_colors = ['b';'r';'g';'k';'b';'r';'g';'k';'b';'r';'g';'k';'b';'r';'g';'k';'b';'r';'g';'k';'b';'r';'g';'k';'b';'r';'g';'k';'b';'r';'g';'k';]; % Color sequence



% Begin: User selects file to process
[file_name, path_name] = uigetfile('*_Spikes.mat','Select MU to sort','Z:\Ferret Data\VisualDiscrim_Ephys_Analysis\');
save_subdir = 'Spike Sorting\';
%Depth Data CSV
inCSV = importdata([path_name 'Electrodes.csv']);
temp_ind = strfind(file_name,'_Spikes');
file_id = file_name(1:temp_ind-1);
tr_value = str2num(file_name(end-4));
numFile = file_name(end-4);
% allInd = load([path_name file_name(1:end-11) 'Indices' numFile '.mat'], ['allInd' numFile]);
% allInd = eval(['allInd.allInd' numFile]);
% End: User selects file to process

% %determine electrodes from file name
% 
% % fileSN = str2num(file_id(7:11));
% % inds = find(inCSV(:,1) == fileSN);
% % sites = inCSV.data(inds,2);
% % electrodes = zeros(size(sites));
% % for iElec = 1:length(sites)
% %     electrodes(iElec) = find(channelMap == sites(iElec));
% % end
% % clear inds sites

% Begin: Load files
spike_file_name = file_name;
disp(['Begin: Load in ' spike_file_name '.'])
load([path_name spike_file_name]);

%Extract experiment parameters
Fs = spine.spec.Fs;
num_channels =  length(spine.spec.chnames);
traces4 = spikes.spike_traces_4;

all_spikes_waveform = cell(num_channels,1);

% cut spike waveforms to 50 sample windows from sample 32 to 82
for i_chn = 1:num_channels
    %all_spikes_waveform{i_chn} = traces4{i_chn}(:,32:end-20).*1e3;
    all_spikes_waveform{i_chn} = traces4{i_chn}(:,32:end-20).*1e3;
    %all_spikes_waveform{i_chn} = traces4{i_chn}(:,20:end-20).*1e3;
end
clear traces4;

disp(['End: Load in ' spike_file_name '.'])


% End: Load files

% Begin: Reduce number of spikes used for sorting since sorting makes
% template which then gets matched to all spikes
%% get random spikes and align peaks
time_ind = -inf;
minpeak_search_window = 7:15;
tic
% Loop through channels to collate spikes and align by negative peak
for i_chn = 1:num_channels
    num_spikes_sort = 3000; % number of spikes used for sorting
    num_spikes =  size(all_spikes_waveform{i_chn},1); % number of spikes in channel
    num_spikes_sort = min(num_spikes, num_spikes_sort);
    shuffle_index_temp = randperm(num_spikes); % shuffles a list of spikes in current channel
    shuffle_index_temp = shuffle_index_temp(1:num_spikes_sort); % chooses first num_spikes_sort (3000) spikes
    sort_spikes_waveform{i_chn} = all_spikes_waveform{i_chn}(shuffle_index_temp,:); % get all spike waveforms of the selected spikes
    % Supersample data 4x with splines
    % align to peaks
    alignedSpks = zeros(num_spikes_sort,size(sort_spikes_waveform{i_chn},2));
    for i = 1:length(sort_spikes_waveform{i_chn})
% % %         spkwave = sort_spikes_waveform{i_chn}(i,:);
% % %         spl = spline(1:length(spkwave), spkwave);
% % %         xx = linspace(1,length(spkwave),4*length(spkwave));
% % %         spkspline(i,:) = ppval(spl,xx);
% % %         localMin = find(spkspline(i,:) == min(spkspline(i,80:120)));
% % %         temp = 130 - localMin;
% % %         alignedSpks(i, temp + 30 + (1:size(spkspline,2))) = spkspline(i,:);
        localMin = find(sort_spikes_waveform{i_chn}(i,:) == min(sort_spikes_waveform{i_chn}(i,minpeak_search_window))); % find the times of min peaks
        temp = 16 - localMin;
        alignedSpks(i, temp + (1:size(sort_spikes_waveform{i_chn},2))) = sort_spikes_waveform{i_chn}(i,:);
    end
    my_spikes = alignedSpks(:, :);
    % Change number of spikes to sort (not used here)
    my_spikes = sort_spikes_waveform{i_chn};
    my_num_spikes = size(my_spikes,1);
    %spikes_to_sort = my_spikes(1:my_num_spikes,:)';
    sort_spikes_waveform{i_chn} = my_spikes(1:my_num_spikes,:);
    %spike_times_to_sort = my_spike_time(:,1:my_num_spikes);
    clear my_spikes
    %time_ind = max(time_ind, size(all_spikes_waveform{i_chn},2));
    time_ind = max(time_ind, size(sort_spikes_waveform{i_chn},2));
    %num_spikes_sort = 3000;
end
toc
%%
% End: Reduce number of spikes used for sorting since sorting makes
% template which then gets matched to all spikes

% for iChn = 1:60
%     sort_spikes_waveform{iChn} = sort_spikes_waveform{iChn}';
% end
% Begin: Show channels and choose which ones to sort
chn_to_sort_suggested = [];
time_msec = (1:time_ind)./Fs*1000;
time_window = [0 3];
amp_window = [-0.3 0.3];
peak_threshold = 0.04;

figure(1)
figure_name{1} = [file_id ': Raw Spike Waveforms'];
set(gcf,'Name',figure_name{1})
clf

% Creates figure of waveforms within each channel
for i_chn = 1:num_channels
    subplot(num_channels/8,8,i_chn,'align')
    if ~isempty(sort_spikes_waveform{i_chn})
        neg_peaks(i_chn) = min(min(sort_spikes_waveform{i_chn}'));
        pos_peaks(i_chn) = max(max(sort_spikes_waveform{i_chn}'));
        if neg_peaks(i_chn) < -peak_threshold || pos_peaks(i_chn) > peak_threshold
            spike_color = 'r';
            chn_to_sort_suggested = [chn_to_sort_suggested i_chn];
        else
            spike_color = 'k';
        end
        if length(sort_spikes_waveform{i_chn}(:,1)) > 1000
            x=50;
        else
            x=1;
        end
        plot(time_msec, sort_spikes_waveform{i_chn}(1:x:end,:)',spike_color)
        %no idea why this next line doesn't resolve the plotting issue
        %plot(sort_spikes_waveform{1,i_chn},spike_color);
    end
    title(num2str(i_chn))
    xlim(time_window)
    ylim(amp_window)
    % axis off
end
s = suptitle(figure_name{1});
set(s,'Interpreter','none');
set(s,'FontSize',10)

set(gcf,'Position',[680 49 1159 948])

chn_to_sort_suggested2 = chn_to_sort_suggested;

prompt={'Select Channels. No commas.'};
name='Enter as M vector';
numlines=1;

% prompts user to choose which channels to sort
answer=inputdlg(prompt,name,numlines);
chn_to_sort = str2num(answer{1}); %don't use commas
if isempty(chn_to_sort)
    disp('User aborted.')
    return;
end

% End: Show channels and choose which ones to sort



%%

Ts = 1/Fs;
spikes_to_sort_window = 5:40;

for i_chn = 1:length(chn_to_sort)
    
    sort_file_output = [path_name save_subdir file_id '_Sorted_' num2str(chn_to_sort(i_chn)) '.mat'];
    
    if exist(sort_file_output,'file')
        next_step = questdlg('Resort?', ...
            'Channel Sorted', ...
            'Resort','Skip','Skip');
        if strcmp(next_step,'Skip') == 1 %
            continue;
        end
    end
    
    % Prepare vector with waveforms to sort
    spikes_to_sort =  (sort_spikes_waveform{chn_to_sort(i_chn)}(:,spikes_to_sort_window))';
    
    % Set-up variables
    spikes_sorted = [];
    next_free = 1;
    
    
    while 1 % Loop until everything sorted
        
        % Begin: Query number of subclusters
        prompt = {'# subclusters:'}; dlg_title = 'Prepare for overclustering';
        num_lines = [1 50]; def = {'16'};     options.WindowStyle='normal';
        answer = inputdlg(prompt,dlg_title,num_lines,def,options);
        num_c = str2num(answer{1});
        % End: Query number of subclusters
        
        % Make clusters
        [spikes_temp  n_units] = SpikeSort_OverCluster(spikes_to_sort, num_c);
        if isempty(spikes_temp)
            break;
        end
       size(spikes_to_sort)
       size(spikes_temp)
       
        % Begin: Collect sorted units (last "unit" gets discarded)
        for i = 1:n_units - 1 % Last "unit" contains all the remaining spikes
            spikes_sorted{i + next_free - 1} = spikes_temp{i};
        end
        next_free = next_free + n_units -1;
        n_units_sorted  = next_free - 1;
        % End: Collect sorted units
        
        
        % Begin: Show all sorted units
        numfig = ceil(n_units_sorted/8); % Number of figures
        for i_fig = 1:1:numfig
            
            % Prepare figure
            figure(i_fig+10)
            set(gcf,'PaperPositionMode','auto')
            set(gcf,'Name','Sorted units')
            clf
            
            for i_sub = 1:perfig
                i_unit = (i_fig-1)*perfig + i_sub;
                if i_unit <= n_units_sorted
                    subplot(min(perfig, n_units_sorted),1,i_sub,'align')
                    cla
                    plot(spikes_sorted{i_unit},'Color',my_colors(i_unit))
                    axis tight
                    title(['Unit # ' num2str(i_unit) ' Spikes: ' num2str(size(spikes_sorted{i_unit},2)) ])
                    ylim([-0.10 0.100])
                    axis off
                end
            end
            
            set(gcf,'Position',[347+(175*i_fig) 100 175 110*n_units_sorted])
            s = suptitle(['Chn: ' num2str(chn_to_sort(i_chn)) '  Sorted Units ' num2str(i_fig) '/' num2str(numfig)]);
            set(s,'FontSize',10);
            set(gcf,'MenuBar','none')
            set(gcf,'ToolBar','none')
        end
        
        choice = questdlg('Next step?', ...
            'SpikeSort_Overcluster', ...
            'Cluster unassigned waveforms','Done','Done');
        
        % Handle response
        if strcmp(choice,'Done') == 1 % Build cluster with remaining waveforms
            break
        else
            spikes_to_sort = spikes_temp{n_units};
            close all
        end
        
    end % Loop until everything sorted
    
    
    % Now that all the waveforms are clustered, we move on to a final
    % review which allows the user to eliminate and merge units.
    
    while 1
        choice = questdlg('Combine Clusters Here?', ...
            'SpikeSort_Overcluster', 'Yes','No','Yes');
        
        if strcmp(choice,'No') == 1 % Build cluster with remaining waveforms
            break
        else
            %prompt for which clusters to combine
            prompt = {'Clusters to Combine:'}; dlg_title = 'Prepare for overclustering';
            num_lines = [1 50]; def = {'0'};     options.WindowStyle='normal';
            answer = inputdlg(prompt,dlg_title,num_lines,def,options);
            num_c = str2num(answer{1});
            if num_c ~= 0
                figure
                set(gcf,'position',[47    79   416   702])
                hold on
                for iClus = 1:length(num_c)
                    plot(spikes_sorted{num_c(iClus)},my_colors(num_c(iClus)));
                end
                choice = questdlg('Combine?', ...
                    'SpikeSort_Overcluster', 'Yes','No','Yes');
                if strcmp(choice,'Yes') == 1 % Combine groups
%                     oldSorted = spikes_sorted;
%                     spikes_sorted = oldSorted;
                    %clear spikes_sorted;
                    %                 if min(num_c) > 1
                    %                     spikes_sorted{1:min(num_c)-1} = oldSorted{1:min(num_c)};
                    %                 end
                    %                 spikes_sorted{min(num_c)} = [];
                    sort(num_c);
                    for iClus = 2:length(num_c)
                        spikes_sorted{min(num_c)} = [spikes_sorted{min(num_c)} spikes_sorted{num_c(iClus)}];
                    end
                    for iClus = 2:length(num_c)
                        spikes_sorted(num_c(iClus)) = [];
                    end
                    n_units_sorted = n_units_sorted - (length(num_c)-1);
                    n_units = n_units - (length(num_c)-1);
                    %clear oldSorted;
                    close all
                    % Begin: Show all sorted units
                    numfig = ceil(n_units_sorted/8); % Number of figures
                    for i_fig = 1:1:numfig
                        
                        % Prepare figure
                        figure(i_fig+10)
                        set(gcf,'PaperPositionMode','auto')
                        set(gcf,'Name','Sorted units')
                        clf
                        
                        for i_sub = 1:perfig
                            i_unit = (i_fig-1)*perfig + i_sub;
                            if i_unit <= n_units_sorted
                                subplot(min(perfig, n_units_sorted),1,i_sub,'align')
                                cla
                                plot(spikes_sorted{i_unit},'Color',my_colors(i_unit))
                                axis tight
                                title(['Unit # ' num2str(i_unit) ' Spikes: ' num2str(size(spikes_sorted{i_unit},2)) ])
                                ylim([-0.10 0.100])
                                axis off
                            end
                        end
                        
                        set(gcf,'Position',[347+(175*i_fig) 100 175 110*n_units_sorted])
                        s = suptitle(['Chn: ' num2str(chn_to_sort(i_chn)) '  Sorted Units ' num2str(i_fig) '/' num2str(numfig)]);
                        set(s,'FontSize',10);
                        set(gcf,'MenuBar','none')
                        set(gcf,'ToolBar','none')
                    end
                end
            end
        end
    end
        
    % Begin: Plot all units
    close all
    numfig = ceil(n_units_sorted/4); % Number of figures
    for i_fig = 1:1:numfig
        % Prepare figure
        figure(i_fig+500)
        clf
        
        for i_sub = 1:perfig
            i_unit = (i_fig-1)*perfig + i_sub;
            if i_unit <= n_units_sorted
                subplot(min(perfig, n_units_sorted),1,i_sub)
                cla
                plot(spikes_sorted{i_unit},'Color',my_colors(i_unit))
                axis tight
                title(['Unit # ' num2str(i_unit) '  Spikes: ' num2str(size(spikes_sorted{i_unit},2)) ])
                ylim([-0.20 0.20])
                axis off
            end
        end
        set(gcf,'Position',[347+(175*i_fig) 100 175 110*n_units_sorted])
        
        set(gcf,'PaperPositionMode','auto')
        
        s = suptitle(['Review units Chn: ' num2str(chn_to_sort(i_chn))]);
        set(s,'FontSize',10)
        set(gcf,'MenuBar','none')
        set(gcf,'ToolBar','none')
    end
    % End: Plot all units
    
    % Begin: Query user about fate of units
    for i=1:n_units_sorted
        prompt{i} = ['Unit ' num2str(i) ':']; def{i} = num2str(i);
    end
    dlg_title = 'Review units'; num_lines = [1 30]; options.WindowStyle='normal';
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    % End: Query user about fate of units
    
    
    % Begin: Process user input to form final unit clusters
    index = 1;
    spikes_sorted_final = cell(n_units_sorted,1);
    %indsSortedFinal = ;
    for i=1:n_units_sorted
        if str2num(answer{i}) ~= 0 % check if unit is eliminated by user
            temp_final_cl = str2num(answer{i});
            % Merge units as requested by user
            for i_unit = temp_final_cl
                spikes_sorted_final{index} = [spikes_sorted_final{index} spikes_sorted{i_unit}];
            end
            index = index + 1;
        end
    end
    spikes_sorted_final(index:end) = [];
    num_sorted_units = max((index-1),1);
    % End: Process user input to form final unit clusters
    
    % Begin: Plot final units
    close all
    numfig = ceil(num_sorted_units/4);
    
    for i_fig = 1:1:numfig
        figure(i_fig+100)
        clf
        for i_sub = 1:perfig
            i_unit = (i_fig-1)*perfig + i_sub;
            
            if i_unit <= num_sorted_units && length(spikes_sorted_final)>0
                subplot(min(perfig, num_sorted_units),1,i_sub)
                cla
                plot(spikes_sorted_final{i_unit},'Color',my_colors(i_unit))
                axis tight
                title(['Unit # ' num2str(i_unit) '  Spikes: ' num2str(size(spikes_sorted_final{i_unit},2)) ])
                if i_unit == num_sorted_units
                    title(['Unit # ' num2str(i_unit) '  Spikes: ' num2str(size(spikes_sorted_final{i_unit},2)) ])
                end
                ylim([-0.20 0.200])
                axis off
            end
        end
        set(gcf,'MenuBar','none')
        set(gcf,'ToolBar','none')
        set(gcf,'Position',[347 100 273 110*num_sorted_units])
        set(gcf,'PaperPositionMode','auto')
        set(gcf,'Name',['Saved Units: ' num2str(i_fig) '/' num2str(numfig)])
        s = suptitle(['Chn: ' num2str(chn_to_sort(i_chn))]);
        set(s,'FontSize',10)
    end 
    
    
    save(sort_file_output,'spikes_sorted_final',  'num_sorted_units')
    disp('Hit any key to continue...')
    pause
    
    close all;
    
end

display('Complete!')
