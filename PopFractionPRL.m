function AllCoef=PopFractionPRL(ref,window,region,volatility,regressors,varargin)

p=inputParser;
addParameter(p,'std',0);
addParameter(p,'smoothbin',500);
addParameter(p,'samplebin',250);
addParameter(p,'save',1);
parse(p,varargin{:});

%% load FR data
ref_text = {'target','','feedback'};

filename = ['PRL_both_' ref_text{ref} '_' num2str(window(1)) '_' num2str(p.Results.smoothbin) ...
    '_' num2str(p.Results.samplebin) '_' num2str(window(end)) '_spkcounts_norm' num2str(p.Results.std) '.mat'];

if ~size(dir(filename),1)
    [Data,Var]=generatePRLData(ref,window,'smoothbin',p.Results.smoothbin, ...
        'samplebin',p.Results.samplebin,'method','spkcounts','norm',p.Results.std);
else
    load(filename)
end

%% Regression

for iid=1:length(regressors)
    varOI(iid) = find(strcmp(Var,regressors{iid}));
end

for ii = 1:length(region)
    for jj = 1:length(volatility)
        
        iData = Data(Data(:,1)==region(ii) & Data(:,3)==volatility(jj),:);
        Cells = unique(iData(:,2));
        binNum = size(Data,2)-length(Var);%ceil(length(window)/p.Results.samplebin);

        Coef.Beta = zeros(length(varOI)+1,binNum,length(Cells));
        Coef.pVal = zeros(length(varOI)+1,binNum,length(Cells));
        Coef.SSE = zeros(binNum,length(Cells));
        
        errorcell=[];
        for c=1:length(Cells)
            ind=iData(:,2)==Cells(c);
            F = iData(ind,varOI);
            FR = iData(ind,length(Var)+1:end);
            lastwarn('', '');
            
            for bin=1:binNum
                FR_t = FR(:,bin);
                mdl = fitlm(F,FR_t);
                Coef.Beta(:,bin,c) = mdl.Coefficients.Estimate;  % mdl put intercept as the first
                Coef.pVal(:,bin,c) = mdl.Coefficients.pValue;
                Coef.SSE(bin,c) = mdl.SSE;
            end
                           
            [~, warnId] = lastwarn();
            if(~isempty(warnId))
                errorcell=[errorcell;c];
            end
            
        end

        Coef.sig = Coef.pVal<0.05;

        for var =1:length(varOI)+1
            var_sig=squeeze(Coef.sig(var,:,~ismember(1:length(Cells),errorcell)));
            if size(var_sig,2)>1
                sigFrac(:,var)=mean(var_sig,2);
            else
                sigFrac(:,var)=mean(var_sig);
            end
        end

        Coef.Fraction = sigFrac;
        
        Coef.Var = [{'intercept'} regressors];
        Coef.Errorcell = errorcell;
        AllCoef{ii,jj}=Coef;
        
    end
end

if p.Results.save
    save(['PRL_' ref_text{ref} '_' num2str(window(1)) '_' num2str(p.Results.smoothbin) '_' ...
        num2str(p.Results.samplebin) '_' num2str(window(end)) '_std' num2str(p.Results.std) '_PopSummary'],'AllCoef');
end
