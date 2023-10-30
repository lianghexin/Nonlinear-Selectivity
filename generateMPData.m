function [Data,Var]=generateMPData(ref,window,varargin)

% Re-organize MP data, ref=1 is target onset, ref=3 is feedback onset
p=inputParser;
addParameter(p,'saveFlag',1);
addParameter(p,'smoothbin',500);
addParameter(p,'samplebin',50);
addParameter(p,'method','spkcounts');
addParameter(p,'norm',0);
parse(p,varargin{:});

% MPtable, col1 area, col2 animal, col3 filename, col4 oldareacode
% 0 DLPFC, 1 SEF, 2 ACC, 3 LIP

%%
load('C:\Users\liang\Documents\GitHub\ThesisProject\NeuralAnalysis\matlab\MP_tables\alldata.mat','newmtx')
MPtable = readtable('C:\Users\liang\Documents\GitHub\ThesisProject\NeuralAnalysis\matlab\MP_tables\MPtable.xlsx');
MPtable = table2cell(MPtable);
ref_text = {'target','','feedback'};
Data = [];
for region=0:3
    areaData = [];
    region_ind = cell2mat(MPtable(:,1))==region;
    animals = cell2mat(MPtable(region_ind,2));
    
    for animal=1:length(animals)
        
        fileind = region_ind & cell2mat(MPtable(:,2))==animals(animal);
        file = MPtable{fileind,3};
        Stable = struct2array((load(['C:\Users\liang\Documents\GitHub\ThesisProject\NeuralAnalysis\matlab\MP_tables\' file], [file(1:end-9) 'Stable'])));
        SPtable = struct2array((load(['C:\Users\liang\Documents\GitHub\ThesisProject\NeuralAnalysis\matlab\MP_tables\' file], [file(1:end-9) 'SPtable'])));
        
        cid = unique(newmtx(newmtx(:,1)==region & newmtx(:,26)==animals(animal),20)); % only include the cells in Hyojung's preprocessed dataset
        trial_id = ismember(Stable(:,1),cid) & Stable(:,3)==2 & Stable(:,18)~=-1; % only include choice task, and exclude choice(t-1) is -1
        
        Stable = Stable(trial_id,:);
        SPtable = SPtable(trial_id,:);
        
        spkMtx = [];
        Allcells = unique(Stable(:,1));
        
        for cell=1:length(Allcells)
            
            trialInd = Stable(:,1)==Allcells(cell);
            cData = SPtable(trialInd,:);
            spk=[];
            
            for trial=1:length(cData)
                spkTime = cData{trial,3};
                spk(trial,:)= ismember((window(1)-p.Results.smoothbin):(window(end)+p.Results.smoothbin),spkTime-Stable(trial,ref+6));
            end
            
            spk = double(spk);
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
            
            spkMtx = [spkMtx; smthFR];
        end
        Loc = Stable(:,21)*2-1;
        PreLoc = Stable(:,18)*2-1;
        Rwd = Stable(:,23)*2-1;
        PreRwd = Stable(:,20)*2-1;
        RL = Rwd.*Loc;
        PRL = PreRwd.*PreLoc;
        LocInter = Loc.*PreLoc;
        areaData = [areaData; ones(length(spkMtx),1)*[region animals(animal)]  Stable(:,1:3) ...
            Loc PreLoc Rwd PreRwd PRL RL LocInter spkMtx];
    end
    cellid = findgroups(areaData(:,2),areaData(:,3));
    areaData(:,3) = cellid;
    Data = [Data; areaData];
end

Var = {'area','animal','cellid','TN','task','Loc','PreLoc','Rwd','PreRwd','PRL','RL','LocInter'};
save(['MP_' ref_text{ref} '_' num2str(window(1)) '_' num2str(p.Results.smoothbin) '_' num2str(p.Results.samplebin) '_'...
    num2str(window(end)) '_' p.Results.method '_norm' num2str(p.Results.norm) '.mat'],'Data','Var');
