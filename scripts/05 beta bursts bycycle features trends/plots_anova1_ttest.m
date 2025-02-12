close all;clear all;clc;
%%
PATH_DATA='Z:\DBS';

DATE=datestr(now,'yyyymmdd');
format long
% from previous analysis:
A=2; % median distance stim1 onset - syl 1 onset
B=1.26; % median distance syl1 onset - syl3 offset
IT=1.93; % median inter trial distance (syl3 offset - stim 1 onset)
LEFT=1.5;
RIGHT=1.5; 

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

%% 
ii=1:numel(SUBJECTS);

for i=ii
    %open a subject dir
    if i==6;continue;end
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
    
    PATH_ANNOT=strcat(... 'Preprocessed data\Sync\annot');
    PATH_SIGNAL=strcat(... 'Preprocessed data\FieldTrip\');

    % take sessions with DBS LFP data
    session=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_session'));
    id_session=session.id(strcmpi(session.type, 'LEAD'));

    % upload useful annot tables
    coding=bml_annot_read(strcat(...,'_coding')); coding=coding(coding.session_id==id_session,:);
    cue=bml_annot_read(strcat(...,'_cue_precise')); cue=cue(cue.session_id==id_session, :);
    prod_triplet=bml_annot_read(strcat(...,'_produced_triplet'));prod_triplet=prod_triplet(prod_triplet.session_id==id_session,:);
    electrode=bml_annot_read(strcat(...,'_electrode'));
    electrode=electrode(:, {'id', 'starts','ends','duration','electrode','connector','port','HCPMMP1_label_1','HCPMMP1_weight_1'});
    
    load(fullfile(...+"_clean_syl.mat"));

    cfg=[];
    cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                % threshold_subj or threshold_weight
    electrode=bml_getEcogArea(cfg,electrode);

    tokeepCod=(~isnan(coding.syl1_onset));
    tokeepCue=(~isnan(cue.stim1_starts));
    switch string(SUBJECTS(i))
        case {'3003', '3010', '3015', '3020', '3025', '3027'};tokeepCue=tokeepCod;
        case {'3012','3014'};tokeepCod=tokeepCue;
        case {'3024'};tokeepCod=tokeepCod & tokeepCue;tokeepCue=tokeepCod;
        case {'3021'};idx=find(~tokeepCue)+1;idx=idx(1:2);tokeepCod(idx)=0;
    end
    center=coding.syl1_onset(tokeepCod);

    x=center - cue.stim1_starts(tokeepCue); %TRIAL SPECIFIC
    tf_baseline_t0=[-(x+0.9) -(x+0.1)];

    x=mean(coding.syl3_offset(tokeepCod) - center,'omitnan'); %MEDIATING
    tf_rebound_t0=[x+0.1 x+0.9];

    x_stim1=median(cue.stim1_starts(tokeepCue)-center,'omitnan');
    x_stim3_end=median(cue.stim3_ends(tokeepCue)-center,'omitnan');
    x_syl1=median(coding.syl1_onset(tokeepCod)-center,'omitnan');
    x_syl3_end=median(coding.syl3_offset(tokeepCod)-center,'omitnan');

    windows={'baseline','stimulus','prespeech','speech','rebound'};
    tag='window50';
    timing_windows=struct();
    for w=1:numel(windows)
        window=windows{w};
        switch window
            case 'baseline'
                timing_windows.(window).start= mean(tf_baseline_t0(:,1),'omitnan');
                timing_windows.(window).end= mean(tf_baseline_t0(:,2),'omitnan');
            case 'stimulus'
                timing_windows.(window).start= x_stim1;
                timing_windows.(window).end= x_stim3_end;
            case 'prespeech'
                timing_windows.(window).start= x_stim3_end;
                timing_windows.(window).end= x_syl1;
            case 'speech'
                timing_windows.(window).start= x_syl1;
                timing_windows.(window).end= x_syl3_end;
            case 'rebound'
                timing_windows.(window).start= tf_rebound_t0(1);
                timing_windows.(window).end= tf_rebound_t0(2);
        end
    end

    tab_cycle=readtable(strcat('annot/general CTAR/',SUBJECT,'_cycle_bursts.txt'));
    % filter out bursts <=100ms 
    tab_cycle=tab_cycle(tab_cycle.duration>0.10001,:); 

    %% channels selection and referencing
    % ecog
    choice='ecog';
    switch choice
        case 'area'
            channels=electrode.HCPMMP1_area;
            channels=channels(~(cellfun('isempty',channels)));
            idx=find(cellfun(@(x) strcmp(x,'AAC'),channels));
            % "PMC","POC","IFC","SMC","TPOJ","LTC","DLPFC","IFOC","AAC","IPC"
            if sum(idx)==0
                warning('No channels in this area');
            else
                cfg=[];
                cfg.channel     =chan_selected;
                ECOG=ft_selectdata(cfg,E);
            end
        case 'ecog'
            cfg=[];
            cfg.channel     ={'ecog*'};
            ECOG=ft_selectdata(cfg,E);
    end
    cfg=[];
    cfg.label   = electrode.electrode;
    cfg.group   = electrode.connector;
    cfg.method  = 'CTAR';   % using trimmed average referencing
    cfg.percent = 50;       % percentage of 'extreme' channels in group to trim
    ECOG=bml_rereference(cfg,ECOG);
    

    % dbs
    cfg=[];
    cfg.channel     ={'dbs_L*'};
    DBS=ft_selectdata(cfg,E);
    if strcmp(SUBJECT,"DBS3028")
        cfg=[];
        cfg.channel=['all', strcat('-', {'dbs_LS','dbs_L10'})];
        DBS=ft_selectdata(cfg,DBS);
    end
    DBS=DBS_referencing(DBS);

    % unify set of data
    cfg=[];
    E=ft_appenddata(cfg,ECOG,DBS);

    chan_selected=E.label;
    electrode=electrode(ismember(electrode.electrode,ECOG.label),:);

    % remasking NaNs with zeros
    cfg=[];
    cfg.value           =0;
    cfg.remask_nan      =true;
    cfg.complete_trial  =true;
    E=bml_mask(cfg,E);

    %% burst probability BYCYCLE
    measures={'frequency','probability','duration','volt_amp','time_rdsym','time_ptsym'};%power

    %time = linspace(-(A+LEFT), B+RIGHT, 250 );
    time=linspace(timing_windows.baseline.start-0.05,timing_windows.rebound.end+0.05,126);
    electrodes_analysis= unique(tab_cycle.chan_id);
    for e=1:numel(electrodes_analysis)
        if startsWith(electrodes_analysis{e},'dbs_R');continue;end
        tab = tab_cycle( strcmp(tab_cycle.chan_id,electrodes_analysis(e)),:);
        trial_ids=unique(tab.trial_id);figure
        
        % select single channel 
        cfg=[];
        cfg.channel     ={electrodes_analysis{e}}; 
        Echannel=ft_selectdata(cfg,E);
        fig1=figure('units','normalized','outerposition',[0.01 0.01 0.99 0.99]);
        if startsWith(electrodes_analysis{e},'ecog')
            area_name= string( electrode.HCPMMP1_area_fullname( strcmp(electrode.electrode,electrodes_analysis{e}) ) );
            sgtitle(strcat(SUBJECT," ",strrep(electrodes_analysis{e},'_',' ')," ",area_name))
        else
            sgtitle(strcat(SUBJECT," ",strrep(electrodes_analysis{e},'_',' ')))
        end

        for m=1:numel(measures)
            disp(['-- analysis: ',measures{m}])
            subplot(2,3,m)
           
            
            xpos=nan(1,numel(windows));
            measure_windows=cell(1,numel(windows));
            for w=1:numel(windows)
                % mask the window
                start_window= timing_windows.(windows{w}).start;
                end_window= timing_windows.(windows{w}).end;
                [~,idx_start]=min(abs(time-start_window));
                [~,idx_end]=min(abs(time-end_window));
                window_time=time(idx_start:idx_end);
        
                % mask for each trial
                switch measures{m}
                    case 'probability'
                    if w==1
                        matrix_mask=nan( numel(trial_ids),length(time) );
                        % continuous probabiliy
                        for j=1:numel(trial_ids)
                            mask=false(size(time));
                            subtab=tab( strcmp(tab.chan_id,electrodes_analysis(e)) & tab.trial_id==trial_ids(j),:);
                            for pp=1:height(subtab)
                                [~,idx_start]=min(abs(time-subtab.starts(pp)));
                                [~,idx_end]=min(abs(time-subtab.ends(pp)));
                                mask(idx_start:idx_end)=1;
                            end
                            matrix_mask(j,:) = mask;
                        end
                        std = sqrt(var(matrix_mask, 0, 1,'omitnan'));
                        std_error = std/sqrt(numel(trial_ids));
                        prob_continuous = mean(matrix_mask);
                        uppbound = prob_continuous + std_error;
                        lowbound = prob_continuous - std_error;
                        plot(time,prob_continuous)
                        %legend('Burst probability')
                        %legend('AutoUpdate','off')
                        ylim([0,0.9])
                        hold on
                        patch([time, fliplr(time)], [uppbound, fliplr(lowbound)],'b','FaceAlpha', 0.1, 'EdgeColor', 'none')
                    end

                    measure_trials=nan(1,numel(trial_ids));
                    for j=1:numel(trial_ids)
                        subtab=tab( strcmp(tab.chan_id,electrodes_analysis(e)) & tab.trial_id==trial_ids(j)  & strcmp(tab.(tag),windows(w)),:);
                        mask=false(size(window_time));
                        try
                            for pp=1:height(subtab)
                                [~,idx_start]=min(abs(window_time-subtab.starts(pp)));
                                [~,idx_end]=min(abs(window_time-subtab.ends(pp)));
                                mask(idx_start:idx_end)=1;
                            end 
                        catch
                        end
                        measure_trials(j) = sum(mask) / length(mask);
                    end


                    case {'duration','volt_amp','time_rdsym','time_ptsym','frequency'}
                    if w==1
                        matrix_it_mean=nan( numel(trial_ids),length(time) );
                        % continuous duration
                        for j=1:numel(trial_ids)
                            it_mean=zeros(size(time));
                            n_elements=zeros(size(time));
                            subtab=tab( strcmp(tab.chan_id,electrodes_analysis(e)) & tab.trial_id==trial_ids(j),:);
                            for pp=1:height(subtab)
                                [~,idx_start]=min(abs(time-subtab.starts(pp)));
                                [~,idx_end]=min(abs(time-subtab.ends(pp)));
                                for idx=idx_start:idx_end
                                    % iterative mean
                                    % https://dev.to/cipharius/calculating-mean-iteratively-3aah
                                    if it_mean(idx)==0
                                        n_elements(idx)=2;
                                        it_mean(idx)=subtab.(measures{m})(pp); 
                                    else
                                        n_elements(idx)=n_elements(idx)+1;
                                        it_mean(idx) = it_mean(idx) * (n_elements(idx)-1)/n_elements(idx) + subtab.(measures{m})(pp)/n_elements(idx);
                                    end
                                end
                            end
                            it_mean(it_mean==0)=nan;
                            matrix_it_mean(j,:) = it_mean;
                        end
                        std = sqrt(var(matrix_it_mean, 0, 1,'omitnan'));
                        std_error = std/sqrt(numel(trial_ids));
                        mean_continuous = mean(matrix_it_mean,'omitnan');
                        uppbound = mean_continuous + std_error;
                        lowbound = mean_continuous - std_error;
                        validIndices = ~isnan(uppbound) & ~isnan(lowbound);
                        uppbound=uppbound(validIndices);lowbound=lowbound(validIndices);
                        plot(time,mean_continuous)
                        hold on
                        if sum(isnan(uppbound))~=0
                            ciao=1;
                        end
                        patch([time(validIndices), fliplr(time(validIndices))], [uppbound, fliplr(lowbound)],'b','FaceAlpha', 0.1, 'EdgeColor', 'none')
                    end  
                    measure_trials=cell(1,numel(trial_ids));
                    for j=1:numel(trial_ids)
                        subtab=tab( strcmp(tab.chan_id,electrodes_analysis(e)) & tab.trial_id==trial_ids(j)  & strcmp(tab.(tag),windows(w)),:);
                        m_trial=nan(1,height(subtab));
                        try
                            for pp=1:height(subtab)
                                m_trial(pp)=subtab.(measures{m})(pp);
                            end 
                        catch
                        end
                        if isempty(m_trial);measure_trials{j}=nan;
                        else;measure_trials{j}=m_trial;
                        end
                    end
                    measure_trials=cell2mat(measure_trials);


                    case 'power'
                    measure_trials=nan(1,numel(trial_ids));
                    for j=1:numel(trial_ids)
                        cfg=[];
                        cfg.trials=[ 1:numel(Echannel.trial)]==trial_ids(j);
                        Etrial=ft_selectdata(cfg,Echannel);
    
                        sig_trial=Etrial.trial{1}(Etrial.time{1}>start_window & Etrial.time{1}<end_window);
    
                        window_pxx=256;noverlap_pxx=128;nfft_pxx=512;fs_pxx=Etrial.fsample;
                        [pxx,freq_pxx]=pwelch(sig_trial,window_pxx,noverlap_pxx,nfft_pxx,fs_pxx);
                        idx_freq=freq_pxx>12 & freq_pxx<35;
                        measure_trials(j)=mean(pxx(idx_freq),'omitnan');
                    end

                end
                measure_windows{w}=measure_trials;
                xpos(w)=end_window-(end_window-start_window)/2;
                if m==numel(measures)
                    outliers=quantile(measure_trials,0.95);% 3.5*(quantile(measure_trials_high,0.75) - quantile(measure_trials_high,0.25));
                    boxplot(measure_trials(abs(measure_trials)<outliers) ,'Positions',xpos(w));hold on
                else
                    boxplot(measure_trials,'Positions',xpos(w))
                end
                
            end

            ylimits=ylim;
            xlim([time(1) time(end)])
            fill([timing_windows.baseline.start timing_windows.baseline.end timing_windows.baseline.end timing_windows.baseline.start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], [1 0.8 0.4], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            fill([timing_windows.stimulus.start timing_windows.stimulus.end timing_windows.stimulus.end timing_windows.stimulus.start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], [0.3 0.4 0.9], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            fill([timing_windows.prespeech.start timing_windows.prespeech.end timing_windows.prespeech.end timing_windows.prespeech.start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)],[0.5 0.8 0.4], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            fill([timing_windows.speech.start timing_windows.speech.end timing_windows.speech.end timing_windows.speech.start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], [1 0.3 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            fill([timing_windows.rebound.start timing_windows.rebound.end timing_windows.rebound.end timing_windows.rebound.start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], [1 0.6 0.4], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            set(gca, 'XTick', xpos, 'XTickLabel', windows);

            marker_pos = ylimits(1) + 0.55*(ylimits(2)-ylimits(1)); 
            step_line = 0.03*(ylimits(2)-ylimits(1));
            step_marker = 0.005*(ylimits(2)-ylimits(1));

            dataConcat = [];
            windowClass = [];
            for ii = 1:length(measure_windows)
                dataConcat = [dataConcat, measure_windows{ii}];
                windowClass = [windowClass, ii * ones(1, length(measure_windows{ii}))];
            end
            [p_value, tbl, stats] = anovan(dataConcat, {windowClass}, 'model', 'linear', 'varnames', {'Time window'},'display', 'off');
            F_value_all = tbl{2, 6};
            p_value_all=p_value;

%                 dataConcat = [];
%                 windowClass = [];
%                 trialClass = [];
%                 for ii = 3:length(measure_windows)
%                     dataConcat = [dataConcat, measure_windows{ii}];
%                     windowClass = [windowClass, ii*ones(1, length(measure_windows{ii}))];
%                 end
%                 [p_value, tbl] = anovan(dataConcat, {windowClass}, 'model', 'linear', 'varnames', {'Time windows'},'display', 'off');
%                 F_value_speech = tbl{2, 6};
%                 p_value_speech = p_value;
            
            %if p_value_all<0.05 && p_value_speech<0.05; back_color=[1,0.8,0.6]; % orange 
            if p_value_all<0.05; back_color=[1,0.8,0.8]; % red
            %elseif p_value_speech<0.05;back_color=[1,1,0.6]; % yellow
            else; back_color='white';
            end
            
%                 text(-2.9,ylimits(1) + 0.82*(ylimits(2)-ylimits(1)), sprintf('Anova1: \nAll: F = %.2f _ p = %.4f \nSpeech: F = %.2f _ p = %.4f',...
%                     F_value_all,p_value_all,F_value_speech,p_value_speech),...
%                     'FontSize',9, 'FontWeight','bold','BackgroundColor', back_color, 'EdgeColor', 'black')
            text(-2.9,ylimits(1) + 0.85*(ylimits(2)-ylimits(1)), sprintf('Anova1: \nF = %.2f _ p = %.4f ',...
                F_value_all,p_value_all),...
                'FontSize',9, 'FontWeight','bold','BackgroundColor', back_color, 'EdgeColor', 'black')

            % Tukey-Kramer post hoc test, has mean, confidence interval and significance
            % columns: 1,2 --> groups in comparison, 
            %           3 --> mean difference, 
            %           4 --> CI for mean difference,
            %           5 --> t value (estimate difference / std)
            %           6 --> p value 
            if p_value < 0.05
                disp('Significant ANOVA. Post-hoc test:');
                tk_results = multcompare(stats, 'CType', 'tukey-kramer','Dimension', [1],'display', 'off');
                
                % plot a marker if significant comparison is found
                for tk=1:height(tk_results)
                    if tk_results(tk,6) < 0.05
                        if tk_results(tk,6) < 0.001;n_markers=3;
                        elseif tk_results(tk,6) < 0.01;n_markers=2;
                        else; n_markers=1;
                        end
                        x1=xpos(tk_results(tk,1));
                        x2=xpos(tk_results(tk,2));

                        plot([x1, x2], [marker_pos, marker_pos], 'k-', 'LineWidth',0.2)
                        markers=repmat('*',1,n_markers);
                        text(mean([x1,x2]), marker_pos+step_marker ,markers,'HorizontalAlignment', 'center', 'FontSize',12)
                        
                        marker_pos = marker_pos + step_line;
                    end
                end  
            end
            title(strrep(measures{m},'_',' '))
        end
        saveas(fig1,strcat('...',SUBJECT," ",strrep(electrodes_analysis{e},'_',' ')," ",measures{m},'.png'));
        saveas(fig1,strcat('...',SUBJECT," ",strrep(electrodes_analysis{e},'_',' ')," ",measures{m},'.fig'));
        close all
    end   
end
