function [Data,Var]=generatePRLData(ref,window,varargin)

p=inputParser;
addParameter(p,'saveFlag',1);
addParameter(p,'smoothbin',500);
addParameter(p,'samplebin',50);
addParameter(p,'method','gaussian');
addParameter(p,'norm',0);
parse(p,varargin{:});

%% Data 2015

% B{session#}
% col 1: Absolute trial count (counting aborted trials)
% col 2: Successful trial count (only increment when successful)
% col 3: Correct  = 1; Aborted = 0
% col 4: Target onset time
% col 5: Magnitude onset time
% col 6: Feedback onset time
%%%%%%%: 'ref' value range 1-3, 1 means target onset, whose column number is 4

% Behavior data from 2015 binnedData_old
% 4.	Chosen Side: (L=0; R=1)
% 5.	Chosen Color: (G=0; R=1)
% 6.	Reward: (no reward = 0; reward = 1)
% 7.	Green Probability
% 8.    Red Probability
% 9.	Green Magnitude  %%%%    Magnitude combinations are [1 1; 1 2; 1 4; 1 8; 2 1; 2 4; 4 1; 4 2; 4 4; 8 1]
% 10.	Red Magnitude
% 11.	Volatility (0 = Low; 1 = High)
% 12.	Monkey (1 =Oscar; 2 = Uba)
% 13.	Good Tr (1 = good trial; 0 = previous fixation break, targets were revealed)

load C:\Users\liang\Documents\Thesis\NeuralAnalysis\OldData\Donahue_2015\Donahue_2015\donahue2015.mat
load C:\Users\liang\Documents\Thesis\NeuralAnalysis\OldData\Donahue_2015\Donahue_2015\binnedData_old\PT_neurData.mat
Data2015=[];
cellcount=1;
for ses=1:length(spkMtx)
    
    beh=B{1,ses}(B{1,ses}(:,3)==1,:);  % only include completed trials here
    
    for cell=1:length(spkMtx{1,ses})
        cellTrialRange=trialRange{1,ses}{1,cell};
        cellbeh=beh(cellTrialRange(1):cellTrialRange(2),:);
        sesData=allData(allData(:,2)==ses & allData(:,1)==cellcount,1:13);
        goodTrialID=sesData(:,13)==1;
        spkTime=spkMtx{1,ses}{1,cell};
        
        if (cellTrialRange(2)-cellTrialRange(1)+1) ~= size(sesData,1)
            fprintf(['cell' num2str(cell) ':trial number inconsistent'])
        end
        
        Loc = sesData(:,4)*2-1; Col = sesData(:,5)*2-1; Rwd = sesData(:,6)*2-1;
        PreLoc =[nan;Loc(1:end-1)]; PreCol =[nan;Col(1:end-1)]; PreRwd =[nan;Rwd(1:end-1)];
        Loc = Loc(goodTrialID,:);
        Col = Col(goodTrialID,:);
        Rwd = Rwd(goodTrialID,:);
        PreLoc = PreLoc(goodTrialID,:);
        PreCol = PreCol(goodTrialID,:);
        PreRwd = PreRwd(goodTrialID,:);
        PRL = PreLoc.*PreRwd;
        RL = Loc.*Rwd;
        POS = Loc.*Col;
        PreCol_POS = PreCol.*POS;
        HVL = PreCol_POS.*PreRwd;
        PRC = PreCol.*PreRwd;
        RC = Col.*Rwd;
        ColInter = Col.*PreCol;
        LocInter = Loc.*PreLoc;
        SwitchHVL = Loc.*HVL;
        
        RedMag = sesData(goodTrialID,10); 
        GrnMag = sesData(goodTrialID,9); 
        LMag = nan(sum(goodTrialID),1); 
        LMag(POS==-1,:) = RedMag(POS==-1); LMag(POS==1,:) = GrnMag(POS==1);
        RMag=RedMag+GrnMag-LMag;
        ChosenMag = nan(sum(goodTrialID),1);
        ChosenMag(Col==-1,:) = GrnMag(Col==-1,:);
        ChosenMag(Col==1,:) = RedMag(Col==1,:);
        UnchosenMag=RedMag+GrnMag-ChosenMag;
                
        for trial=1:length(sesData)
            spk(trial,:)= ismember((window(1)-p.Results.smoothbin):(window(end)+p.Results.smoothbin),spkTime-cellbeh(trial,ref+3));
        end
        spk = double(spk(goodTrialID,:));
        
        if strcmp(p.Results.method,'spkcounts')
            smthFR = movsum(spk,p.Results.smoothbin,2);
        else
            smthFR = smoothdata(spk,2,p.Results.method,p.Results.smoothbin)*1000;
        end
        
        smthFR = smthFR(:,p.Results.smoothbin+1:end-p.Results.smoothbin);
        smthFR = downsample(smthFR',p.Results.samplebin)';
        
        if p.Results.norm
            allFRmean=mean(smthFR,'all');
            allFRstd=std(smthFR,0,'all');
            smthFR = (smthFR-allFRmean)/allFRstd;
        end  
        
        Data2015 = [Data2015; zeros(sum(goodTrialID),1) cellcount*ones(sum(goodTrialID),1) ones(sum(goodTrialID),1)...
            Loc PreLoc RL PRL LocInter Col PreCol RC PRC ColInter Rwd PreRwd POS RMag LMag RedMag GrnMag ChosenMag UnchosenMag HVL SwitchHVL PreCol_POS smthFR];
        cellcount=cellcount+1;
        
    end
end

clearvars -except Data2015 p ref window

%% Data 2018

% B{session#}
% col 1: Absolute trial count (counting aborted trials)
% col 2: Successful trial count (only increment when successful)
% col 3: Correct  = 1; Aborted = 0
% col 4: Target onset time
% col 5: Magnitude onset time
% col 6: Feedback onset time
%%%%%%%: 'ref' value range 1-3, 1 means target onset, whose column number is 4

% Behavior data from 2018 behMtx_session for S.mtx:
% 1. area (0: dlpfc, 1: ofc, 2: acc).
% 2. Monkey (1: Uba, 2: Xavier)
% 3. cell number (count resets for each region)
% 4. Session Number (count resets for each region)
% 5. Trial Number (starts at 1, resets for each neuron)
% 6: Trial in Block (starts at 1, resets every block change)
% 7. Trial Type (0: Low Color, 1: Low Shape, 2: High Vol)
% 8. Volatility (0: Low, 1: High)
% 9. rewarded (0: NR, 1: Rew).
% 10. Chosen Side (0: Left, 1: Right).
% 11. Chosen Color (1: Green, 2: Red, 3: orange/square, 4: cyan/diamond)
%     [Note: This is flipped for Xavier]. Both monkeys are in terms of high value target (Lval = 3; Hval = 4).
% 12. Green (or Low) probability
% 13. Red (or High) probability
% 14. Green (or Low) Magnitude
% 15. Red (or High) Magnitude
% 16. Color of Feedback Ring (0: Grey, 1: Cyan, 2: Green, 3: Red)
% 17. chamber radius
% 18. chamber angle
% 19. Channel
% 20. Unit
% 21. Recording Quality (1: Poor, 2: Fair, 3: Good, 4: Excellent)
% 22. Recording Depth

load C:\Users\liang\Documents\Thesis\NeuralAnalysis\OldData\Massi2018_behavior\dataset_old\behMtx_neuron.mat

if ref==1
    rawspikes=S.spikes.targ;
elseif ref==3
    rawspikes=S.spikes.fb;
end

spk=cellfun(@(x) ismember((window(1)-p.Results.smoothbin):(window(end)+p.Results.smoothbin),x),rawspikes,'UniformOutput',0);
spk=double(cell2mat(spk'));

smthFR = [];
allcellid = findgroups(S.mtx(:,1),S.mtx(:,3));
for cell=1:length(S.neuron)
    
    cellSpk = spk(allcellid==cell,:);
    if strcmp(p.Results.method,'spkcounts')
        cell_smthFR = movsum(cellSpk,p.Results.smoothbin,2);
    else
        cell_smthFR = smoothdata(cellSpk,2,p.Results.method,p.Results.smoothbin)*1000;
    end
    
    cell_smthFR = cell_smthFR(:,p.Results.smoothbin+1:end-p.Results.smoothbin);
    cell_smthFR = downsample(cell_smthFR',p.Results.samplebin)';
    
    if p.Results.norm
        allFRmean=mean(cell_smthFR,'all');
        allFRstd=std(cell_smthFR,0,'all');
        cell_smthFR = (cell_smthFR-allFRmean)/allFRstd;
    end
    
    smthFR = [smthFR;cell_smthFR];
    
end

blockInd=[0; diff(S.mtx(:,6))]~=1;
Loc=S.mtx(:,10)*2-1; Rwd=S.mtx(:,9)*2-1;
highVolInd = S.mtx(:,8)==1;
Col(highVolInd,1)=S.mtx(highVolInd,11)*2-3;
Col(~highVolInd,1)=S.mtx(~highVolInd,11)*2-7;

PreLoc =[nan;Loc(1:end-1)]; PreCol =[nan;Col(1:end-1)]; PreRwd =[nan;Rwd(1:end-1)];
PreLoc(blockInd)=nan; PreCol(blockInd)=nan; PreRwd(blockInd)=nan;
PRL = PreLoc.*PreRwd;
RL = Loc.*Rwd;
POS = Loc.*Col;
PreCol_POS = PreCol.*POS;
HVL = PreCol.*PreRwd.*POS;
PRC = PreCol.*PreRwd;
RC = Col.*Rwd;
ColInter = Col.* PreCol;
LocInter = Loc.* PreLoc;
SwitchHVL = Loc.*HVL;

RedMag = S.mtx(:,15); GrnMag = S.mtx(:,14); 
LMag = nan(length(RedMag),1); 
LMag(POS==-1,:) = RedMag(POS==-1); LMag(POS==1,:) = GrnMag(POS==1);
RMag=RedMag+GrnMag-LMag;
ChosenMag = nan(length(RedMag),1);
ChosenMag(Col==-1,:) = GrnMag(Col==-1,:);
ChosenMag(Col==1,:) = RedMag(Col==1,:);
UnchosenMag=RedMag+GrnMag-ChosenMag;
        
Data2018 = [S.mtx(:,[1 3 8]) Loc PreLoc RL PRL LocInter Col PreCol RC PRC ColInter ...
    Rwd PreRwd POS RMag LMag RedMag GrnMag ChosenMag UnchosenMag HVL SwitchHVL PreCol_POS smthFR];

%% combine two datasets

Data2018(Data2018(:,1)==0,2) = Data2018(Data2018(:,1)==0,2) + length(unique(Data2015(:,2)));
Data = [Data2015; Data2018];

Var={'area','cellid','Volatility','Loc', 'PreLoc', 'RL', 'PRL', 'LocInter', 'Col', 'PreCol', 'RC', 'PRC', 'ColInter', ...
    'Rwd', 'PreRwd', 'POS', 'RMag', 'LMag', 'RedMag', 'GrnMag', 'ChosenMag', 'UnchosenMag', 'HVL', 'SwitchHVL', 'PreCol_POS'};

if p.Results.saveFlag
    ref_text={'target','magnitude','feedback'};
    save(['PRL_' ref_text{ref} '_' num2str(window(1)) '_' num2str(p.Results.smoothbin) '_' num2str(p.Results.samplebin) '_'...
         num2str(window(end)) '_' p.Results.method '_norm' num2str(p.Results.norm) '.mat'],...
        'Data','Var');
end

