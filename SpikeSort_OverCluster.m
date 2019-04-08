function [sorted_spikes, num_units] = MEA_SpikeSort_OverCluster(my_spikes, num_c)

disp('Creating subclusters using kmeans')
pause(0.1)
tic
[clusters C]= kmeans(my_spikes',num_c, 'MaxIter',100);
toc
sc_spikes = cell(num_c,1);


% Collect spikes for each subcluster
for j=1:num_c

    sc_spikes{j} = [my_spikes(:,clusters==j)];
end


figure(300)
set(gcf,'Name','Linkage Analysis')
clf
Y = pdist(C);
dendrogram(linkage(Y),0);
set(gcf,'Position',[17 679 457 284])
set(gcf,'PaperPositionMode','auto')



% Plot parameters
perfig = 16; % Number of panels per figure
numfig = ceil(num_c/16); % Number of figures
for k = 1:1:numfig
    figure(k+1000)
        set(gcf,'Position',[747 26 681 906])
        set(gcf,'PaperPositionMode','auto')
        
    for tt = 1:perfig
        j = (k-1)*perfig + tt;
        if tt == 1
            clf
        end
        if j <= num_c
            subplot(perfig/2,2,tt,'align')
            if size(sc_spikes{j}) < 100
                plot(sc_spikes{j},'k')
            else
                plot(sc_spikes{j}(:,1:5:end),'k')
            end
            axis tight
            text(10,0.06,['Cluster # ' num2str(j) '  Spikes: ' num2str(size(sc_spikes{j},2))], 'FontSize',8)
            ylim([-0.100 0.100])
            axis off
            set(gcf,'Position',[50+(k*300)    71   300   906]);
            set(gcf,'PaperPositionMode','auto')
        end
    end
end

prompt = {'Unit 1:','Unit 2:', 'Unit 3:', 'Unit 4:','Unit 5:','Unit 6:','Unit 7:'};
dlg_title = 'Combine Subclusters';
num_lines = 1;
def = {'0','0','0','0','0','0','0'};
options.WindowStyle='normal';


while 1==1
    
    disp('Combine subclusters')
    
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    if isempty(answer)
        sorted_spikes = [];
        n_units = [];
        break;
    end
    
    ref_period = 0.001;
    allused_sc = [];
    my_clusters_temp = [];
    
    
    def = {'0','0','0','0','0','0','0'};
    for i=1:length(answer)
        
        if  str2num(answer{i}) > 0
            
            my_clusters_temp{i} = str2num(answer{i});
            allused_sc = [allused_sc my_clusters_temp{i}];
            def{i} = num2str(my_clusters_temp{i});
        end
        
    end
    
    num_units = length(my_clusters_temp);
    my_clusters_temp{num_units+1} = setdiff([1:1:num_c],allused_sc );
    def{num_units+1} = num2str(my_clusters_temp{num_units+1});
    num_units = num_units + 1;
    
    sorted_spikes = cell(num_units,1);

    sorted_isi = cell(num_units,1);
    
   
    for j=1:num_units % includes unassigned waveforms
        single_cluster = my_clusters_temp{j}  ;          
        % Loop through subclusters
        for i=1:length(single_cluster)

            sorted_spikes{j} = [sorted_spikes{j} my_spikes(:,clusters==single_cluster(i))];
        end
    end
    
    
    
   
%     figure(300)
%     clf
%     Y = pdist(C);
%     dendrogram(linkage(Y),0);
%     set(gcf,'Position',[4 34 1394 902])
%     set(gcf,'PaperPositionMode','auto')
    
    perfig = 4; % Number of panels per figure
    numfig = ceil(num_units/4); % Number of figures
    my_colors = ['b';'r';'g';'k';'b';'r';'g';'k'];
    for k = 1:1:numfig
        figure(k+200)
        clf
        set(gcf,'Position',[420+(400*k) 34 400 910])
        set(gcf,'PaperPositionMode','auto')
        
        for tt = 1:perfig
            j = (k-1)*perfig + tt;
            
            if j <= num_units
                subplot(perfig,1,tt)
                cla
                plot(sorted_spikes{j},'Color',my_colors(j))
                axis tight
                % text(20,0.02,['# ' num2str(j) '  Spikes: ' num2str(size(sorted_spikes{j},2)) '  RefViol: ' num2str(ref_viol(j))])
                %  title(['# ' num2str(j) '  Spikes: ' num2str(size(sorted_spikes{j},2)) '  RefViol: ' num2str(ref_viol(j))])
                title(['Unit # ' num2str(j) '  Spikes: ' num2str(size(sorted_spikes{j},2)) ])
                if j == num_units
                    title(['Unit # ' num2str(j) ' (Unassigned Waveformes) Spikes: ' num2str(size(sorted_spikes{j},2)) ])
                end
                ylim([-0.10 0.100])
                % axis off
            end
        end
    end
    
    % Construct a questdlg with three options
    
    
    
    
    choice = questdlg('Happy with Units?', ...
        'SpikeSort_Overcluster', ...
        'Yes','No','No');
    % Handle response
    if strcmp(choice,'Yes') == 1
        return
    end
    
end

disp('...done.')
end


