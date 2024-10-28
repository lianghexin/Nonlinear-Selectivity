%% global params
ref_text = {'target','','feedback'};

% PRL CPD

FullRegressors = {'Loc','PreLoc','RL','PRL','LocInter','Col','PreCol','RC','PRC','ColInter',....
    'Rwd','PreRwd','POS','ChosenMag','UnchosenMag','LMag','HVL','SwitchHVL'};

for ref = 1

    if ref==1
        epoch=[1000 1001];
    else
        epoch=[0 1];
    end

    filename = ['PRL_' ref_text{ref} '_' num2str(epoch(1)) '_1000_1000_' num2str(epoch(2)) '_std0_PopSummary.mat'];

    if ~size(dir(filename),1)
        AllCoef=PopFractionPRL(ref,epoch,[0 1 2],[0 1],FullRegressors,'std',0, ...
            'smoothbin',1000,'samplebin',1000);
    else
        load(filename)
        fprintf('Fullcoef exists')
    end

    %
    CPD_val = cell(3,2,length(FullRegressors));

    for vv=1:length(FullRegressors)
        vv
        varOI=1:length(FullRegressors)~=vv;
        reduced_var = FullRegressors{vv};
        regressors = FullRegressors(varOI);

        filename=['C:\Users\liang\Documents\GitHub\Nonlinear-Selectivity\CPD\PRL_' ref_text{ref} '_' reduced_var '_1000_std0_PopSummary.mat'];

        if ~size(dir(filename),1)
            Coef=PopFractionPRL(ref,epoch,[0 1 2],[0 1],regressors,'std',0,'save',0, ...
                'smoothbin',1000,'samplebin',1000);
            save(filename,'Coef')
        else
            load(filename)
            fprintf('reduced coef exists')
        end

        for reg = 1:3
            ind = ~ismember(1:size(AllCoef{reg,1}.Beta,3),AllCoef{reg,1}.Errorcell);
            for volatility = 1:2
                SSE_full = AllCoef{reg,volatility}.SSE(:,ind);
                SSE_red = Coef{reg,volatility}.SSE(:,ind);
                CPD_val{reg,volatility,vv} = (SSE_red-SSE_full)./SSE_red*100;
            end
        end

    end

    save(['PRL_' ref_text{ref} '_1000_std0_CPD'],'CPD_val','FullRegressors');
end

%% MP

model=1;
FullRegressors = {'Loc','PreLoc','RL','PRL','LocInter','Rwd','PreRwd'};

ref_text = {'target','','feedback'};
% ref = 1;
% FullCoef = PopFractionMP(ref,-500:1500,[0 1 2 3],['model' num2str(model)],FullRegressors,'std',1);
ref = 3;
FullCoef = PopFractionMP(ref,[0 1],[0 1 2 3],FullRegressors,'std',0, ...
    'smoothbin',1000,'samplebin',1000);

%%
ref=3;
load('MP_feedback_0_1000_1000_1_std0_PopSummary.mat'); % PopFractionMP(ref,-500:1500,[0 1 2 3],['model' num2str(model)],FullRegressors,'save',1);
FullCoef=AllCoef;
FullRegressors = {'Loc','PreLoc','RL','PRL','LocInter','Rwd','PreRwd'};
ref_text = {'target','','feedback'};

%
CPD_mean = cell(1,4);
CPD_sem = cell(1,4);
CPD_val=cell(4,length(FullRegressors));

for vv=1:length(FullRegressors)
    vv
    varOI=1:length(FullRegressors)~=vv;
    reduced_var = FullRegressors{vv};
    regressors = FullRegressors(varOI);

    filename = ['MP_' ref_text{ref} '_model1_' reduced_var '_1000_std0_PopSummary.mat'];
    if ~size(dir(filename),1)
        Coef = PopFractionMP(ref,[0 1],[0 1 2 3],regressors,'save',0, ...
            'smoothbin',1000,'samplebin',1000);
    else
        Coef = load(filename);
        Coef = Coef.AllCoef;
    end

    for reg = 1:4
        SSE_full = FullCoef{reg}.SSE;
        SSE_red = Coef{reg}.SSE;
        CPD_val{reg,vv} = (SSE_red-SSE_full)./SSE_red*100;
        %         CPD_mean{reg}(:,vv) = nanmean(CPD_val,2);
        %         CPD_sem{reg}(:,vv) = nanstd(CPD_val,0,2)/sqrt(size(CPD_val,2));
    end

end

save(['MP_' ref_text{ref} '_1000_std0_CPD'],'CPD_val','FullRegressors');


%% looking at choice signals in reduced models, MP

varsOI={'Loc','PreLoc'};
ts=2;

data = load('C:\Users\liang\Documents\GitHub\Nonlinear-Selectivity\CPD\MP_feedback_model1_LocInter_-500_500_250_1500_std0_PopSummary');
AllCoef = data.AllCoef;
Var = AllCoef{1,1}.Var;
Allsig = cell(4,1);
Fraction = zeros(4,3);

for region=1:4
    sig = zeros(length(varsOI),length(ts),size(AllCoef{region}.Beta,3));
    for iid=1:length(varsOI)
        var = strcmp(Var,varsOI{iid});
        sig(iid,:,:) = squeeze(AllCoef{region}.sig(var,ts,:));
    end
    sig = squeeze(sig);
    sig(3,:) = sig(1,:).*sig(2,:);
    Allsig{region}=sig;
    Fraction(region,:)=mean(sig,2);
end

Fraction(:,4)=Fraction(:,1).*Fraction(:,2);

number=[322 185 154 205];
zscore_add=zeros(4,1);
p_add=zeros(4,1);
for ii=1:4
    [zscore_add(ii),p_add(ii)]=prop_test([Fraction(ii,3) Fraction(ii,4)],number(ii),'z','two');
end


%% PRL, remove all interaction terms, ColInter, LocInter, switchHVL

ref_text = {'target','','feedback'};

for ref = [1 3]

    regressors = {'Loc','PreLoc','RL','PRL','Col','PreCol','RC','PRC',....
        'Rwd','PreRwd','POS','ChosenMag','UnchosenMag','LMag','HVL'};
    
    filename=['C:\Users\liang\Documents\GitHub\Nonlinear-Selectivity\CPD\PRL_' ref_text{ref} '_interaction_-500_500_250_1500_std0_PopSummary.mat'];
    
    Coef=PopFractionPRL(ref,[-500 1500],[0 1 2],[0 1],regressors,'std',0,'save',0);
    save(filename,'Coef')

end

%% PRL target
% Col 

data = load('C:\Users\liang\Documents\GitHub\Nonlinear-Selectivity\CPD\PRL_target_interaction_-500_500_250_1500_std0_PopSummary');

varsOI={'Col','PreCol'};
ts=7;
AllCoef = data.Coef;
Var = AllCoef{1,1}.Var;
Allsig = cell(3,1);
Fraction = zeros(3,3);

for region=1:3
    sig = zeros(length(varsOI),length(ts),size(AllCoef{region,2}.Beta,3));
    for iid=1:length(varsOI)
        var = strcmp(Var,varsOI{iid});
        sig(iid,:,:) = squeeze(AllCoef{region,2}.sig(var,ts,:));
    end
    sig = squeeze(sig);
    sig(3,:) = sig(1,:).*sig(2,:);
    Allsig{region}=sig;
    Fraction(region,:)=mean(sig,2);
end

Fraction(:,4)=Fraction(:,1).*Fraction(:,2);

% Loc
varsOI={'Loc','PreLoc'};
Allsig = cell(3,1);
Fraction = zeros(3,3);
for region=1:3
    sig = zeros(length(varsOI),length(ts),size(AllCoef{region,2}.Beta,3));
    for iid=1:length(varsOI)
        var = strcmp(Var,varsOI{iid});
        sig(iid,:,:) = squeeze(AllCoef{region,2}.sig(var,ts,:));
    end
    sig = squeeze(sig);
    sig(3,:) = sig(1,:).*sig(2,:);
    Allsig{region}=sig;
    Fraction(region,:)=mean(sig,2);
end

Fraction(:,4)=Fraction(:,1).*Fraction(:,2);


%% PRL feedback
% Col 

data = load('C:\Users\liang\Documents\GitHub\Nonlinear-Selectivity\CPD\PRL_feedback_interaction_-500_500_250_1500_std0_PopSummary');

varsOI={'Col','PreCol'};
ts=3;
AllCoef = data.Coef;
Var = AllCoef{1,1}.Var;
Allsig = cell(3,1);
Fraction = zeros(3,3);

for region=1:3
    sig = zeros(length(varsOI),length(ts),size(AllCoef{region,2}.Beta,3));
    for iid=1:length(varsOI)
        var = strcmp(Var,varsOI{iid});
        sig(iid,:,:) = squeeze(AllCoef{region,2}.sig(var,ts,:));
    end
    sig = squeeze(sig);
    sig(3,:) = sig(1,:).*sig(2,:);
    Allsig{region}=sig;
    Fraction(region,:)=mean(sig,2);
end

Fraction(:,4)=Fraction(:,1).*Fraction(:,2);

% Loc
varsOI={'Loc','PreLoc'};
Allsig = cell(3,1);
Fraction = zeros(3,3);
for region=1:3
    sig = zeros(length(varsOI),length(ts),size(AllCoef{region,2}.Beta,3));
    for iid=1:length(varsOI)
        var = strcmp(Var,varsOI{iid});
        sig(iid,:,:) = squeeze(AllCoef{region,2}.sig(var,ts,:));
    end
    sig = squeeze(sig);
    sig(3,:) = sig(1,:).*sig(2,:);
    Allsig{region}=sig;
    Fraction(region,:)=mean(sig,2);
end

Fraction(:,4)=Fraction(:,1).*Fraction(:,2);

%%

number=[400 135 135];
zscore_add=zeros(3,1);
p_add=zeros(3,1);

% value=[0.0100	0.0120;	
% 0.0074	0.0099;	
% 0.0074	0.0054];

value=[0.0450	0.0371;	
0.0000	0.0025;	
0.0074	0.0014];

for ii=1:3
    [zscore_add(ii),p_add(ii)]=prop_test([value(ii,1) value(ii,2)],number(ii),'z','two');
end

%% look at temporal evolution
FullRegressors = {'Loc','PreLoc','RL','PRL','LocInter','Col','PreCol','RC','PRC','ColInter',....
    'Rwd','PreRwd','POS','ChosenMag','UnchosenMag','LMag','HVL','SwitchHVL'};
Coef=PopFractionPRL(1,[500 1500],[0],[1],FullRegressors,'std',0,'save',1, ...
                'smoothbin',200,'samplebin',20);

w_col=squeeze(Coef{1,1}.Beta(7,:,:));
w_switch=squeeze(Coef{1,1}.Beta(11,:,:));

[~,ts] = max(Coef{1,1}.Fraction(:,11));

[S,I] = sort(w_switch(ts,:)); 
map = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1);ones(128, 1), linspace(1, 0, 128)', linspace(1,0, 128)'];
h = heatmap(w_col(:,I)');
h.Colormap = map;
h.GridVisible = 'off';
h.ColorLimits = [-0.2 0.2];

coef_hm=zeros(51,51);
for t1=1:51
    for t2=1:51
        cf = corrcoef(w_switch(t1,:),w_switch(t2,:));
        coef_hm(t1,t2) = cf(1,2);
    end
end
close all
map = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1);ones(128, 1), linspace(1, 0, 128)', linspace(1,0, 128)'];
h = heatmap(coef_hm);
h.Colormap = map;
h.GridVisible = 'off';
h.ColorLimits = [-1 1];