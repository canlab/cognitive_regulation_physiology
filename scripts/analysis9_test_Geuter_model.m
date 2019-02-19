%% Testing Geuter SCR model 

clear;
close all;

basedir = '/Users/clinpsywoo/Dropbox/github/MWW_inprep/scripts';
cd(basedir);

load('../data/CRB_dataset_SCR_lpf5Hz_DS25Hz_011516.mat');

figdir = '../figures';
datdir = '../data';

%% Step 1. Creating Y: pain ratings
%%
for i = 1:numel(D.Event_Level.data)
    y_int{i} = D.Event_Level.data{i}(:,12);
    y_unp{i} = -2*(D.Event_Level.data{i}(:,13)-50);
    xx{i} = [D.Event_Level.data{i}(:,11) D.Event_Level.data{i}(:,16) scale(D.Event_Level.data{i}(:,11),1).*scale(D.Event_Level.data{i}(:,16),1)];
    reg{i} = D.Event_Level.data{i}(:,16);
    temp{i} = D.Event_Level.data{i}(:,11);
end

u_temp = unique(temp{1});
u_reg = unique(reg{1});

for subj = 1:numel(temp)
    for i = 1:numel(u_temp)
        for j = 1:numel(u_reg)
            dat_int{i,j}(subj,1) = nanmean(y_int{subj}(temp{subj}==u_temp(i) & reg{subj}==u_reg(j)));
            dat_unp{i,j}(subj,1) = nanmean(y_unp{subj}(temp{subj}==u_temp(i) & reg{subj}==u_reg(j)));
        end
    end
end

%% Step 2. Creating X: SCR data (averaged over three trials)
%%
for j = 1:numel(D.Continuous.data)      % loop through subjects
    
    scr = D.Continuous.data{j};         % SCR data for each person
    
    ons = D.Event_Level.data{j}(:,5);   % Heat onset 
    
    for ii = 1:6
        for jj = 1:3
            signal{j}.cond{ii, jj} = [];  % pre-allocation 
        end
    end
        
    for i = 1:numel(ons)
        temp_lev = D.Event_Level.data{j}(i,11)-43;  % temperature level 1-6
        reg_lev = D.Event_Level.data{j}(i,16)+2;    % regulation level 1-3
        a = round((ons(i)-3)*25);                   % onset in 25 Hz
        signal{j}.cond{temp_lev, reg_lev}(end+1,:) = scr(a:(a+(23*25))); 
                                                    % get data for 23 seconds
    end
end

clear signal_m;
k = 0;

for i = 1:6
    for j = 1:3
        k = k + 1;
        for subj = 1:numel(signal)
            signal_m{i,j}(subj,:) = mean(signal{subj}.cond{i,j} - repmat(mean(signal{subj}.cond{i,j}(:,1:75),2), 1, size(signal{subj}.cond{i,j},2)));
                                                     % subtracting the baseline, and then average
        end
    end
end

%% Step 3. Test Geuter et al. (2014)'s 17-second SCR model
%%
scr = load(fullfile(datdir, 'SCR_prediction_dat_112816.mat'));
load(fullfile(datdir, 'Geuter_SCR_model.mat'));

% Applying the model on regulation trials using leave-one-participant-out cross validation
for i = 1:3
    for j = 1:6
        for subj = 1:size(scr.signal_m{j,i},1)
            test_scr_int{j,i}(subj,1) = scr.signal_m{j,i}(subj,76:575)*Geuter_SCR_model(1:500);
            test_scr_unp{j,i}(subj,1) = scr.signal_m{j,i}(subj,76:575)*Geuter_SCR_model(1:500);
        end
    end
end

temp_int = reshape(scr.pcr_stats.int.yfit, 41, 6);
temp_unp = reshape(scr.pcr_stats.unp.yfit, 41, 6);

for j = 1:6
    test_scr_int{j,2} = temp_int(:,j);
    test_scr_unp{j,2} = temp_unp(:,j);
end

% scatter plot
y_int_passive = cat(2,scr.dat_int{:,2});
yfit_int_passive = cat(2,test_scr_int{:,2});

y_unp_passive = cat(2,scr.dat_unp{:,1});
yfit_unp_passive = cat(2,test_scr_int{:,1});

y_int = [cat(2,scr.dat_int{:,1}) cat(2,scr.dat_int{:,3})];
yfit_int = [cat(2,test_scr_int{:,1}) cat(2,test_scr_int{:,3})];

y_unp = [cat(2,scr.dat_unp{:,1}) cat(2,scr.dat_unp{:,3})];
yfit_unp = [cat(2,test_scr_int{:,1}) cat(2,test_scr_int{:,3})];

%% Step 4. correlations
%%
clear test_por;
for i = 1:size(y_int,1)
    x = y_int(i,:);
    y = yfit_int(i,:);
    test_por_int(i) = corr(x',y');
    
    x = y_unp(i,:);
    y = yfit_unp(i,:);
    test_por_unp(i) = corr(x',y');
    
    x = y_int_passive(i,:);
    y = yfit_int_passive(i,:);
    test_por_int_passive(i) = corr(x',y');
    
    x = y_unp_passive(i,:);
    y = yfit_unp_passive(i,:);
    test_por_unp_passive(i) = corr(x',y');
end

fprintf('\nTest results: mean prediction_outcome_r for intensity (passive) = %1.3f', mean(test_por_int_passive));
fprintf('\nTest results: mean prediction_outcome_r for unpleasantness (passive) = %1.3f', mean(test_por_unp_passive));


fprintf('\nTest results: mean prediction_outcome_r for intensity (regulation) = %1.3f', mean(test_por_int));
fprintf('\nTest results: mean prediction_outcome_r for unpleasantness (regulation)= %1.3f', mean(test_por_unp));
%% Step 5. Test the effects of temperature and regulation on SCR pattern response
% 5-1. multilevel GLM analysis for SCR intensity 

temp = repmat(1:6, 1, 3)';
reg = [ones(6,1)*-1; zeros(6,1); ones(6,1)*1];
for i = 1:size(test_scr_int{1},1)
    xx{i} = [temp reg scale(temp,1).*scale(reg,1)];
    
    kk = 0;
    for k = 1:3
        for j = 1:6
            kk = kk + 1;
            yy_int{i}(kk,1) = test_scr_int{j,k}(i);
            yy_unp{i}(kk,1) = test_scr_unp{j,k}(i);
        end
    end
end
glm_scr_int = glmfit_multilevel(yy_int, xx, [], 'names', {'intcp', 'temp', 'reg', 'intera'}, ...
    'weighted', 'verbose', 'boot', 'nresample', 10000);

glm_scr_unp = glmfit_multilevel(yy_unp, xx, [], 'names', {'intcp', 'temp', 'reg', 'intera'}, ...
    'weighted', 'verbose', 'boot', 'nresample', 10000);
%% 
% 5-2. Test up vs. passive, passive vs. down

temp = repmat(1:6, 1, 3)';
reg = [ones(6,1)*-1; zeros(6,1); ones(6,1)*1];

%% passive vs. down

clear yy_int yy_unp;

for i = 1:size(test_scr_int{1},1)
    reg_pass_down = reg(1:12);
    reg_pass_down(reg_pass_down==0) = 1;
    xx{i} = [temp(1:12) reg_pass_down scale(temp(1:12),1).*scale(reg_pass_down,1)];
    
    kk = 0;
    for k = 1:2
        for j = 1:6
            kk = kk + 1;
            yy_int{i}(kk,1) = test_scr_int{j,k}(i);
            yy_unp{i}(kk,1) = test_scr_unp{j,k}(i);
        end
    end
end

glm_scr_int = glmfit_multilevel(yy_int, xx, [], 'names', {'intcp', 'temp', 'reg', 'intera'}, ...
    'weighted', 'verbose', 'boot', 'nresample', 10000);

glm_scr_unp = glmfit_multilevel(yy_unp, xx, [], 'names', {'intcp', 'temp', 'reg', 'intera'}, ...
    'weighted', 'verbose', 'boot', 'nresample', 10000);

%%
%% up vs. passive

clear yy_int yy_unp;

for i = 1:size(test_scr_int{1},1)
    
    reg_up_pass = reg(7:18);
    reg_up_pass(reg_up_pass==0) = -1;
    xx{i} = [temp(7:18) reg_up_pass scale(temp(7:18),1).*scale(reg_up_pass,1)];
    
    kk = 0;
    for k = 2:3
        for j = 1:6
            kk = kk + 1;
            yy_int{i}(kk,1) = test_scr_int{j,k}(i);
            yy_unp{i}(kk,1) = test_scr_unp{j,k}(i);
        end
    end
end

glm_scr_int = glmfit_multilevel(yy_int, xx, [], 'names', {'intcp', 'temp', 'reg', 'intera'}, ...
    'weighted', 'verbose', 'boot', 'nresample', 10000);

glm_scr_unp = glmfit_multilevel(yy_unp, xx, [], 'names', {'intcp', 'temp', 'reg', 'intera'}, ...
    'weighted', 'verbose', 'boot', 'nresample', 10000);