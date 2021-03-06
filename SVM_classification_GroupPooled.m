%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SVM_classfication_xxx
% This script uses linear support vector machine method to examine if
% neural population activity (su+mu) could decode the touch location or
% perceptual difficulty
%  
% 
%using window (400ms) to check the encoding accuracy around stim onset  or touch using sliding windows of 100ms
% input file: 
%
%output: 
%
% %(c) Frohlich Lab 2016
% Author: Charles Zhou
% Code Check: Yuhui Li 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc
% close all
% clear all
clear accuracy_group_chance accuracy_group decode_accuracy
SessionsSavePath = 'C:\Users\Frohlich Lab\Desktop\VisualDiscrimProject\SciRepRevisionPlots\SVM\Sessions';

%% Declare variables

minTrials = 40; % minimum number of trials needed to include a session in the analysis

fontsize_label=8;
fontsize_tick = 8;
fontsize_title = 10;

%%%%data in behavioral code%%%
Col_TrialNumber = 1;
Col_StimulusWindow = 2;
Col_TrialOnsetToInit = 4; % How long the spout light is on prior to trial initiation
Col_XCoordinateTouch = 5;
Col_YCoordinateTouch = 6;
Col_TouchTimeStampHr = 7;
Col_TouchTimeStampMin = 8;
Col_TouchTimeStampSec = 9;
Col_TrialOnsetToTouch = 10;
Col_HitMiss = 11;   % 0 for miss or 1 for touch; 2 for premature touch
Col_DrinkTimeStampHr = 12;
Col_DrinkTimeStampMin = 13;
Col_DrinkTimeStampSec = 14;
Col_TrialOnsetToDrink = 15;
Col_NumTimeOutTouch = 16; %number of touches during time out period
Col_NumITITouch=17;  %number of touches during ITI period after time out
Col_NumDelayTouch = 18;  %number of touches during waiting period, permature reponses

% define the range and bin size for the FR PSTH. Bin width is determined by
% win_ends(1)-win_starts(1). Binning calculation calculates spikes
% incrementing each win_step

% define epoch and task variable to analyze (for visual discrim)
EpochLabel = {'StimOn', 'Touch'};
TaskVarLabel = {'Diff', 'Loc'};
iEpoch = 2;
TaskVar = 1;
numConds = 2;

if iEpoch == 1 % for stim on centered epoch
    winStep=0.1;
    winStarts=-0.5:winStep:1.7;  
    winEnds=-0.1:winStep:2.1;
    
    CenterTime = 'Stim On';
elseif iEpoch == 2 % for touch
    winStep=0.1;
    winStarts=-2:winStep:2.6;  
    winEnds=-1.6:winStep:3.0;
    CenterTime = 'Touch';
end

numSteps=length(winStarts);
winWidth=winEnds(1)-winStarts(1);

t=(winStarts+winEnds)./2;
tmin=min(t);
tmax=max(t);


TrialsPerSesh = cellfun('length', ValidTrialsPerSession);% calculate min number of trials across sessions
ValidSessions = find( TrialsPerSesh > minTrials ); % take sessions with greater than x number of trials

% get number of trials in each condition
TrialIDs = {};
CondID = {};
NumCondID = [];
for iFile = ValidSessions
    % load in trial identifiers
    if TaskVar == 1
        TrialIDs{iFile} = EasyHardTrial{iFile}+1;
    elseif TaskVar == 2
        TrialIDs{iFile} = LeftRightTrial{iFile}+1;
    end
    
    % find trial indices and number of trials for each condition
    for iID = 1:numConds
       
        CondID{iFile,iID} = find( TrialIDs{iFile} == iID );
        NumCondID(iFile,iID) = length(CondID{iFile,iID});
        
    end
    
end

% find minimum num trials across all sessions and conditions
NumCondID( NumCondID == 0 ) = [];
MinBothCond = min(NumCondID);
numTrialsHit = MinBothCond*2;

% Prep trial IDs to load (dependent on min number of trials in either cond
TrialIDLoad = [ repmat(1,1,MinBothCond) repmat(2,1,MinBothCond) ];
%% Begin: Loop through (load) files
% get random sample of the minimum num of trials for each condition. Orient
% FR data and trial ID data such that even across SU and sessions, first
% trials are of condition1 and second half of trials are cond2.

% to do this, new SU and trial reference numbers are used to load in data
% to trial_FR

decodeAccuracy = [];
decodeAccuracyChance=[];

for Iter = 1:15
    
    clear TrialIDs
    
    trialFR = {};
    SUCount = 1; % used to reference new SU order
    for iFile = ValidSessions
        
        %%% compute histogram  and spike times within each trial
        suSpikesHit={};
        
        % tabulate trial identifiers and calculate FR PSTH
        for iSU=CombinedValidUnits{iFile} % used to reference old SU ID
            
            TrialCount = 1; % used to reference new trial order
            
            for iID = 1:numConds
                
                % get random trials from session for given condition
                RandTrials = randsample( length(CondID{iFile,iID}),MinBothCond );
                TrialsOfInt = CondID{iFile,iID}(RandTrials);
                
                for iTrial=TrialsOfInt
                    
                    % section spike times into SU and trial
                    suSpikesHit{iSU,TrialCount}= SpikeTimes_Trial_Allepoch{iEpoch, iFile, iTrial, iSU}; % SpikeTimes_Trial_Allepoch dim are epoch, file, trial, SU. Inside are spike times in seconds % su_spike_times{i_su_unit,1}(1,:);
                    
                    % calculate FR PSTH
                    for  iStep=1:numSteps
                        Twindow1=winStarts(1)+(iStep-1)*winStep;
                        Twindow2=winEnds(1)+(iStep-1)*winStep;
                        
                        % find number of spikes in this bin and normalize by bin
                        % width
                        trialFR{iStep}(SUCount,TrialCount,1) = length(find( suSpikesHit{iSU,TrialCount}>Twindow1 & suSpikesHit{iSU,TrialCount}<Twindow2 ))./winWidth;
                        trialFR{iStep}(SUCount,TrialCount,2) = TrialIDLoad(TrialCount); % define trial ID.
                        % 0's are left and easy; 1's are rights and hards
                        
                    end
                    TrialCount = TrialCount + 1;
                end
                
                clear RandTrials TrialsOfInt
                
            end
            SUCount = SUCount + 1;
        end
    end
    
    
    %% perform SVM analysis
    performance=[];
    performanceChance=[];
    
    
    predictors={};
    for i=1:SUCount-1
        predictors=[predictors, ['x' num2str(i)] ];
    end
    
    
    for iStep=1:numSteps % cycle through PSTH bins
        
        for iTrial=1:numTrialsHit
            
            disp(['working on bin ' num2str(iStep) ' trial ' num2str(iTrial) ])
            
            % get all trials other than target trial
            tempTrials= 1:numTrialsHit;
            newTrials=tempTrials(tempTrials~=iTrial);
            
            % get FR and trial ID for each bin
            x=trialFR{iStep}(:,newTrials,1)';
            y=TrialIDLoad(newTrials)'; % grab class labels
            
            % convert FR data to table; add FR and trial ID data into
            % struct
            svm_train=array2table(x);
            svm_train.Group=y;
            
            [trainedClassifier, validationAccuracy] = trainClassifierGroup(svm_train,predictors);
            
            clear x
            clear y
            
            % use this trial's data to predict trial condition using the
            % trained classifier
            x=trialFR{iStep}(:,iTrial,1)';
            y=TrialIDLoad(iTrial);
            svm_test=array2table(x);
            yfit = trainedClassifier.predictFcn(svm_test);
            
            % test if prediction is accurate
            if yfit==y
                accuracy=1;
            else
                accuracy=0;
            end
            
            performance(iTrial)=accuracy;
            
        end
        
        decodeAccuracy(Iter,iStep)=length(find(performance==1))./length(performance);
        
        clear x
        clear y
        %%% chance level; shuffle the class labels
        for iTrialCon=1:numTrialsHit
            
            % find all trials that aren't the target trial
            tempTrials= 1:numTrialsHit;
            newTrials=tempTrials(tempTrials~=iTrialCon);
            
            % grab class labels and training data
            x=trialFR{iStep}(:,newTrials,1)';
            y=TrialIDLoad(newTrials);
            y=y(randperm(length(y)))'; %shuffle labels
            
            svmTrainCon=array2table(x);
            svmTrainCon.Group=y;
            
            [trainedClassifier, validationAccuracy] = trainClassifierGroup(svmTrainCon,predictors);
            
            clear x
            clear y
            
            % use this trial's data to predict trial condition using the
            % trained classifier
            x=trialFR{iStep}(:,iTrialCon,1)';
            y=TrialIDLoad(iTrialCon);
            svmTestCon=array2table(x);
            yfit = trainedClassifier.predictFcn(svmTestCon);
            
            if yfit==y
                accuracyChance=1;
            else
                accuracyChance=0;
            end
            
            % populate prediction value for each trial
            performanceChance(iTrialCon)=accuracyChance;
            
        end
        decodeAccuracyChance(Iter,iStep)=length(find(performanceChance==1))./length(performanceChance);
        
    end
end

figure
clf
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[465 521 300 300])
plot(t, decodeAccuracyChance)
hold on
plot(t,decodeAccuracy,'r', 'LineWidth',2)
ylim([0 1])
xlim([tmin tmax])
xlabel( ['Time from ' CenterTime])
ylabel('Classifier Accuracy')
%title(['session' num2str(iFile+10) ' ' num2str(NumSUTaskPref) '/' num2str(NumSU(iFile)) 'With Preference'])
line([tmin tmax],[0.5 0.5],'Color',[0 0 0],'LineStyle','--')
line([0 0],[0 1],'Color',[0 0 0],'LineStyle','--')

% clear decode_accuracy_chance decode_accuracy performance_chance trial_FR


saveas(100+iFile,['S' num2str(iFile) '_' EpochLabel{iEpoch} TaskVarLabel{TaskVar}])

%     figure_file_name = [path_name_figure filesep file_id{iFile} 'align_touch_4windows' ];
%     print(gcf,'-depsc2','-tiff',figure_file_name)
%     print(gcf,'-djpeg100',figure_file_name)



%% Stats for plotting

% define period of time to perform stats on
Epoch = [ 16 26 ]; % FOR CENTER ON TOUCH ONLY NOW
EpochBins = Epoch(1):Epoch(2);

alpha = 0.05./size(decodeAccuracy,2);

% test significant difference between chance and actual
[h_SVM,p_SVM] = ttest( decodeAccuracyChance, decodeAccuracy,'alpha', alpha);

SVMSigBins = find(h_SVM == 1);
SigWinStart = winStarts(SVMSigBins);
SigWinEnd = winEnds(SVMSigBins);
%% plot mean across sessions

figure
clf
 set(gcf,'PaperPositionMode','auto');
    set(gcf,'Position',[465 521 300 300])
%num_file=10;
plot(t, nanmean(decodeAccuracyChance),'LineWidth',2)
errorbar(t, nanmean(decodeAccuracyChance),nanstd(decodeAccuracyChance)./sqrt(NumFiles))
hold on
plot(t, nanmean(decodeAccuracy),'r','LineWidth',2)
errorbar(t, nanmean(decodeAccuracy),nanstd(decodeAccuracy)./sqrt(NumFiles))
% plot significant stats
for iCount = 1:length(SVMSigBins)
    hold on
    line( [ SigWinStart(iCount) SigWinStart(iCount)+winStep ], [0.35 0.35], 'LineWidth',3)
end
title('SVM Group Average')
ylim([0 1])
xlim([tmin tmax])
xlabel( ['Time from ' CenterTime])
ylabel('Classifier Accuracy')
line([tmin tmax],[0.5 0.5],'Color',[0 0 0],'LineStyle','--')
line([0 0],[0 1],'Color',[0 0 0],'LineStyle','--')
legend('Chance','Actual')
% Save
% path_name_toSave = fullfile(path_name,'Spiking_EachSession')
% figure_file_name = [path_name_toSave,filesep num2str(animalID) '_S11-20_align_touch_4winds_AVG']
% print(gcf,'-depsc2','-tiff',figure_file_name)
% print(gcf,'-djpeg100',figure_file_name)

%% prep for stats

% mean and conf int.
[muhat,sigmahat,muci,sigmaci] = normfit(reshape(decodeAccuracy(:,21:31),[],1));  % gives 95% conf int.

% for comparing between SU and population

figure 
hold on
plot(t,  nanmean(decodeAccuracyChance,1), 'b')
plot(t,  nanmean(decodeAccuracy,1), 'r')
ylim([0 1])
xlim([tmin tmax])
% for comparison between SU and population

% get performance difference from chance 
for iFile = 1:NumFiles
    
    PerformanceDiff(iFile,:) = decodeAccuracy(iFile,EpochBins) - decodeAccuracyChance(iFile,EpochBins);
    
end

%%

SavePath = 'C:\Users\FrohlichLab\Desktop\CZSVM\AvgData\';
SaveFile =[SavePath 'Pooled_' CenterTime '_' TaskVarLabel{TaskVar} '_15Iter'];
save ( SaveFile,'decode_accuracy_chance', 'decode_accuracy', 'PerformanceDiff', 'win_step','win_starts','win_ends', 'iEpoch', 'TaskVar', 'h_SVM', 'p_SVM')   

%% 

figure (1401)
clf
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[465 121 1386 857])
for  iFile = 1:NumFiles
    
%          if i_file >=7
%            i_file_plot=i_file+1
%          else
%             i_file_plot=i_file
%          end
%         subplot(3,4,i_file_plot)
   subplot(3,4,iFile)
    plot(t, decodeAccuracyChance(iFile,:))
    hold on
    plot(t,decodeAccuracy(iFile,:),'r')
    ylim([0 1])
    xlim([tmin tmax])
    xlabel('Time from touch')
    ylabel('classifieraccuracy')
   % title(['ch', num2str(i_file_plot)])
   title(['ch', num2str(iFile+10)])
    line([tmin tmax],[0.5 0.5],'Color',[0 0 0],'LineStyle','--')
    line([0 0],[0 1],'Color',[0 0 0],'LineStyle','--')
    
    
%     path_name_toSave = fullfile(path_name,'Spiking_EachSession')
%     figure_file_name = [path_name_toSave,filesep num2str(animalID) '_S11-20_align_touch_4winds_individule']
%     print(gcf,'-depsc2','-tiff',figure_file_name)
%     print(gcf,'-djpeg100',figure_file_name)
end
