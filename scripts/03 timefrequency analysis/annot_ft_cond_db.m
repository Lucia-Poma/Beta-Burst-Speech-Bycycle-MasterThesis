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
TF_RATE = 20; %Hz Sampling rate of time frequency plot
TF_FOI = round(10.^(0.30:0.05:2.4),2,'signif');
TF_BASELINE_WIDTH = 0.5;
beta=[11,22,36];
saltati=0;

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;


for i=1:numel(SUBJECTS)  
    if i==6;continue;end
    %open a subject dir
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))

    % paths
    NAME_RAW_SIGNAL=strcat(SUBJECT,'_ft_raw_session');
    PATH_SIGNAL=strcat(..., 'Preprocessed data\FieldTrip\');
    PATH_ANNOT=strcat(..., 'Preprocessed data\Sync\annot');

    % take sessions with DBS LFP data
    session=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_session'));
    id_session=session.id(strcmpi(session.type, 'LEAD'));

    % upload useful annot tables
    coding=bml_annot_read(strcat(...,SUBJECT,'_coding')); coding=coding(coding.session_id==id_session,:);
    cue=bml_annot_read(strcat(...SUBJECT,'_cue_precise')); cue=cue(cue.session_id==id_session, :);
    prod_triplet=bml_annot_read(strcat(...,SUBJECT,'_produced_triplet'));prod_triplet=prod_triplet(prod_triplet.session_id==id_session,:);
    electrode=bml_annot_read(strcat(...,SUBJECT,'_electrode'));
    electrode=electrode(:, {'id', 'starts','ends','duration','electrode','connector','port','HCPMMP1_label_1','HCPMMP1_weight_1'});
    
    tokeepCod=(~isnan(coding.syl1_onset));
    tokeepCue=(~isnan(cue.stim1_starts));
    switch string(SUBJECTS(i))
        case {'3003', '3010', '3015', '3020', '3025', '3027'};tokeepCue=tokeepCod;
        case {'3014','3012'};tokeepCod=tokeepCue;
        case {'3024'};tokeepCod=tokeepCod & tokeepCue;tokeepCue=tokeepCod;
        case {'3021'};idx=find(~tokeepCue)+1;idx=idx(1:2);tokeepCod(idx)=0;
    end
    
    % load data
    load(fullfile(...,SUBJECT+"_clean_syl.mat"));a=-(A+LEFT);b=B+RIGHT;center=coding.syl1_onset(tokeepCod);
    
    % define epoch, baseline and rebound for ft
    tf_epoch_t0=E.time{:};tf_epoch_t0=[tf_epoch_t0(1) tf_epoch_t0(end)];
    tf_epoch_t0=round(tf_epoch_t0 .* TF_RATE) ./ TF_RATE;

    x=center - cue.stim1_starts(tokeepCue); %TRIAL SPECIFIC
    % x=median(coding.syl1_onset - cue.stim1_starts,'omitnan'); %MEDIATING
    tf_baseline_t0=[-(x+0.9) -(x+0.1)];
    tf_baseline_t0=round(tf_baseline_t0 .* TF_RATE) ./ TF_RATE;

    % x=coding.syl3_offset(tokeepCod) - center; %TRIAL SPECIFIC
    x=mean(coding.syl3_offset(tokeepCod) - center,'omitnan'); %MEDIATING
    tf_rebound_t0=[x+0.1 x+0.9];
    tf_rebound_t0=round(tf_rebound_t0 .* TF_RATE) ./ TF_RATE;

    % define stim and syl x
    x_stim1=median(cue.stim1_starts(tokeepCue)-center,'omitnan');
    x_stim2=median(cue.stim2_starts(tokeepCue)-center,'omitnan');
    x_stim3=median(cue.stim3_starts(tokeepCue)-center,'omitnan');
    x_stim3_end=median(cue.stim3_ends(tokeepCue)-center,'omitnan');
    x_stim1_end=median(cue.stim1_ends(tokeepCue)-center,'omitnan');
    x_stim2_end=median(cue.stim2_ends(tokeepCue)-center,'omitnan');
    x_syl1=median(coding.syl1_onset(tokeepCod)-center,'omitnan');
    x_syl2=median(coding.syl2_onset(tokeepCod)-center,'omitnan');
    x_syl3=median(coding.syl3_onset(tokeepCod)-center,'omitnan');
    x_syl1_end=median(coding.syl1_offset(tokeepCod)-center,'omitnan');
    x_syl2_end=median(coding.syl2_offset(tokeepCod)-center,'omitnan');
    x_syl3_end=median(coding.syl3_offset(tokeepCod)-center,'omitnan');

    %% select ecog and referencing
    cfg=[];
    cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                % threshold_subj or threshold_weight
    electrode=bml_getEcogArea(cfg,electrode);
    
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
    cfg.channel     ={'dbs*'};
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
    electrode=electrode(ismember(electrode.electrode,chan_selected),:); 
    
    %% ft
    disp('*************************************************************')
    fprintf('Time Frequency Analysis for subject %s \n',SUBJECT)
    nTrials=numel(E.trial);
    criteria={'Volume','RT','Accuracy'};
    cond_definition={'positive','negative','difference'};
    powtab=table();
    powtab_volp=table();
    powtab_voln=table();
    powtab_vold=table();
    powtab_rtp=table();
    powtab_rtn=table();
    powtab_rtd=table();
    powtab_accp=table();
    powtab_accn=table();
    powtab_accd=table();
 
    % remasking NaNs with zeros
    cfg=[];
    cfg.value           =0;
    cfg.remask_nan      =true;
    cfg.complete_trial  =true;
    E=bml_mask(cfg,E);

    nElec=numel(chan_selected);
    
    baselinetype='db'; % 'difference', 'percentage', 'db', 'zscore'
    for e=1:nElec
        disp('---------------------------------------------------------')
        fprintf('%s %s',SUBJECT,chan_selected{e})
        
        % select single channel 
        cfg=[];
        cfg.channel     ={chan_selected{e}}; 
        Esingle=ft_selectdata(cfg,E);

        data=cell2mat(Esingle.trial');
        if sum(data(:)==0)/numel(data)>=0.7;disp('Electrode rejected, going to the next');continue;end

        % condition 
        for j=1:numel(criteria)
            switch criteria{j}
            case 'Volume'    
                data=prod_triplet.stim_volume(tokeepCod);
                positive=find(data); % high volume
                negative=find(~data); % low
            case 'RT'
                data=coding.syl1_onset(tokeepCod)-cue.stim1_ends(tokeepCue);
                positive=find(data>=median(data,'omitnan')); % high RT
                negative=find(data<median(data,'omitnan')); % low
            case 'Accuracy'
                try 
                    data=prod_triplet.phonetic_ontarget(tokeepCod);
                    positive=find(data); % high accuracy
                    negative=find(~data); % low
                catch 
                    saltati=saltati+1;
                    continue
                end
            end
            disp(strcat(SUBJECT,": ",chan_selected{e}," ",criteria{j}))
            for k=1:numel(cond_definition)
                if k<3
                    warning off
                    cfg=[];
                    if k==1;cfg.trials=positive';tit=strcat('N trials: ',num2str(length(positive)));
                    elseif k==2;cfg.trials=negative';tit=strcat('N trials: ',num2str(length(negative)));end
                    Econd=ft_selectdata(cfg,Esingle); 
                    
                    % tf analysis
                    
                    cfg=[];
                    cfg.foi     =TF_FOI;
                    cfg.dt      =1/TF_RATE;
                    cfg.toilim  =tf_epoch_t0;
                    width=7;
                    cfg.width   =width;
                    Econd_mor=bml_freqanalysis_power_wavelet(cfg,Econd);
                    warning on
                    
                    % normalize data
                    cfg=[];
                    cfg.baselinetype=baselinetype;
                    cfg.baselinedef=tf_baseline_t0;
                    Econd_norm=baseline_normalization(cfg,Econd_mor);
                    
                    % xtime to plot (only original epoch)
                    [~,aidx]=min(abs(Econd_norm.time-a));
                    [~,bidx]=min(abs(Econd_norm.time-b));
                    Econd_norm.time=Econd_norm.time(aidx:bidx);
                    Econd_norm.powspctrm=Econd_norm.powspctrm(:,aidx:bidx);
                    
                    if k==1;Epos=Econd_norm;
                    elseif k==2;Eneg=Econd_norm;end
                    Etoplot=Econd_norm;
                else
                    Etoplot=Epos;Etoplot.powspctrm=Epos.powspctrm-Eneg.powspctrm;
                end
                
                % mean power in beta during task and rebound
                pow_beta=Etoplot.powspctrm( Etoplot.freq>=beta(1)&Etoplot.freq<=beta(3) ,:);
                pow_highbeta=Etoplot.powspctrm( Etoplot.freq>=beta(2)&Etoplot.freq<=beta(3) ,:);
                pow_lowbeta=Etoplot.powspctrm( Etoplot.freq>=beta(1)&Etoplot.freq<=beta(2) ,:);
                powtab.id=e;
                powtab.label=chan_selected(e);
                    %baseline 
                [~,idx1]=min(abs(Etoplot.time-a));
                [~,idx2]=min(abs(Etoplot.time-b));
                powtab.base_beta=mean(pow_beta(:, idx1:idx2),'all','omitnan');
                powtab.base_highbeta=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
                powtab.base_lowbeta=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');
                    % stim
                [~,idx1]=min(abs(Etoplot.time-x_stim1));
                [~,idx2]=min(abs(Etoplot.time-x_stim3_end));
                powtab.stim_beta=mean(pow_beta(:, idx1:idx2),'all','omitnan');
                powtab.stim_highbeta=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
                powtab.stim_lowbeta=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');
                    % prespeech
                [~,idx1]=min(abs(Etoplot.time-(x_stim3_end)));
                [~,idx2]=min(abs(Etoplot.time-(x_syl1)));
                powtab.prespeech_beta=mean(pow_beta(:, idx1:idx2),'all','omitnan');
                powtab.prespeech_highbeta=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
                powtab.prespeech_lowbeta=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');
                    % speech
                [~,idx1]=min(abs(Etoplot.time-x_syl1));
                [~,idx2]=min(abs(Etoplot.time-x_syl3_end));
                powtab.speech_beta=mean(pow_beta(:, idx1:idx2),'all','omitnan');
                powtab.speech_highbeta=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
                powtab.speech_lowbeta=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');
                    % rebound
                [~,idx1]=min(abs(Etoplot.time-(x_syl3_end)));
                [~,idx2]=min(abs(Etoplot.time-(x_syl3_end+1)));
                powtab.rebound_beta=mean(pow_beta(:, idx1:idx2),'all','omitnan');
                powtab.rebound_highbeta=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
                powtab.rebound_lowbeta=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan'); 

                if j==1 && k==1
                    powtab_volp(e,:)=powtab;
                elseif j==1 && k==2
                    powtab_voln(e,:)=powtab;     
                elseif j==1 && k==3
                    powtab_vold(e,:)=powtab;
                elseif j==2 && k==1
                    powtab_rtp(e,:)=powtab;
                elseif j==2 && k==2
                    powtab_rtn(e,:)=powtab;
                elseif j==2 && k==3 
                    powtab_rtd(e,:)=powtab;
                elseif j==3 && k==1   
                    powtab_accp(e,:)=powtab;
                elseif j==3 && k==2
                    powtab_accn(e,:)=powtab;
                elseif j==3 && k==3
                    powtab_accd(e,:)=powtab;
                end
            end
        end
    end
    powtab_volp=powtab_volp(~powtab_volp.id==0, :);powtab_voln=powtab_voln(~powtab_voln.id==0, :);powtab_vold=powtab_vold(~powtab_vold.id==0, :);
    powtab_rtp=powtab_rtp(~powtab_rtp.id==0, :);powtab_rtn=powtab_rtn(~powtab_rtn.id==0, :);powtab_rtd=powtab_rtd(~powtab_rtd.id==0, :);
    try powtab_accp=powtab_accp(~powtab_accp.id==0, :);powtab_accn=powtab_accn(~powtab_accn.id==0, :);powtab_accd=powtab_accd(~powtab_accd.id==0, :); catch;end

    powtab_volp.id=(1:height(powtab_volp))';powtab_voln.id=(1:height(powtab_voln))';powtab_vold.id=(1:height(powtab_vold))';
    powtab_rtp.id=(1:height(powtab_rtp))';powtab_rtn.id=(1:height(powtab_rtn))';powtab_rtd.id=(1:height(powtab_rtd))';
    try powtab_accp.id=(1:height(powtab_accp))';powtab_accn.id=(1:height(powtab_accn))';powtab_accd.id=(1:height(powtab_accd))';catch;end

    writetable(powtab_volp, strcat('annot/conditions CTAR/',SUBJECT,'_speech_normdb_volume_pos.txt'), 'Delimiter', '\t');
    writetable(powtab_voln, strcat('annot/conditions CTAR/',SUBJECT,'_speech_normdb_volume_neg.txt'), 'Delimiter', '\t');
    writetable(powtab_vold, strcat('annot/conditions CTAR/',SUBJECT,'_speech_normdb_volume_diff.txt'), 'Delimiter', '\t');
    writetable(powtab_rtp, strcat('annot/conditions CTAR/',SUBJECT,'_speech_normdb_reacttime_pos.txt'), 'Delimiter', '\t');
    writetable(powtab_rtn, strcat('annot/conditions CTAR/',SUBJECT,'_speech_normdb_reacttime_neg.txt'), 'Delimiter', '\t');
    writetable(powtab_rtd, strcat('annot/conditions CTAR/',SUBJECT,'_speech_normdb_reacttime_diff.txt'), 'Delimiter', '\t');
    try writetable(powtab_accp, strcat('annot/conditions CTAR/',SUBJECT,'_speech_normdb_accuracy_pos.txt'), 'Delimiter', '\t');catch;end
    try writetable(powtab_accn, strcat('annot/conditions CTAR/',SUBJECT,'_speech_normdb_accuracy_neg.txt'), 'Delimiter', '\t');catch;end
    try writetable(powtab_accd, strcat('annot/conditions CTAR/',SUBJECT,'_speech_normdb_accuracy_diff.txt'), 'Delimiter', '\t');catch;end
end

disp('finishh :)')
%  pcolor(Esingle_norm.time, Esingle_norm.freq,pow)
%  set(gca,'yscale','log')      




