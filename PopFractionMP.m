function AllCoef=PopFractionMP(ref,window,region,regressors,varargin)

p=inputParser;
addParameter(p,'std',0);
addParameter(p,'smoothbin',500);
addParameter(p,'samplebin',250);
addParameter(p,'save',1);
parse(p,varargin{:});

%% load & preprocess data
ref_text = {'target','','feedback'};
filename = ['MP_' ref_text{ref} '_' num2str(window(1)) '_' num2str(p.Results.smoothbin) ...
    '_' num2str(p.Results.samplebin) '_' num2str(window(end)) '_spkcounts_norm' num2str(p.Results.std) '.mat'];

if ~size(dir(filename),1)
    [Data,Var]=generateMPData(ref,window,'samplebin',p.Results.samplebin,'smoothbin',p.Results.smoothbin, ...
        'method','spkcounts','norm',0);
else
    load(filename)
end

%% Regression

for iid=1:length(regressors)
    varOI(iid) = find(strcmp(Var,regressors{iid}));
end

for ii = 1:length(region)
            
            iData = Data(Data(:,1)==region(ii),:);
            Cells = unique(iData(:,3));
            binNum = ceil(length(window)/p.Results.samplebin);
            Coef.Beta = zeros(length(varOI)+1,binNum,length(Cells));
            Coef.pVal = zeros(length(varOI)+1,binNum,length(Cells));
            Coef.SSE = zeros(binNum,length(Cells));
            
            for c=1:length(Cells)
                
                ind=iData(:,3)==Cells(c);
                F = iData(ind,varOI);
                FR = iData(ind,length(Var)+1:end);
                
                for bin=1:binNum
                    FR_t = FR(:,bin);
                    mdl = fitlm(F,FR_t);                    
                    Coef.Beta(:,bin,c) = mdl.Coefficients.Estimate;  % mdl put intercept as the first
                    Coef.pVal(:,bin,c) = mdl.Coefficients.pValue;
                    Coef.SSE(bin,c) = mdl.SSE;
                end
                
            end
            
            Coef.sig = Coef.pVal<0.05;
            
            for var =1:length(varOI)+1
                var_sig=squeeze(Coef.sig(var,:,:));
                if size(var_sig,2)>1
                    sigFrac(:,var)=mean(var_sig,2);
                else
                    sigFrac(:,var)=mean(var_sig);
                end
            end
            
            Coef.Fraction = sigFrac;
            Coef.Var = [{'intercept'} regressors];
            AllCoef{ii}=Coef;
        
end

save(['MP_' ref_text{ref} '_' num2str(window(1)) '_' num2str(p.Results.smoothbin) '_' ...
    num2str(p.Results.samplebin) '_' num2str(window(end)) '_std' num2str(p.Results.std) '_PopSummary'],'AllCoef');

