function ValRunNKL(f)
% Runs the model after fitting to produce estimates for saturation function f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% f - the saturation function to be used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load([FilenameFor parameterEstimates]); 
           
            
load('DataNorthKivu.mat'); % Load data
NW=length(IData); % Determine how long to run the model for
clear IData; % Clear the data

% Determine the start time of vaccination
startDateofSim = datenum('04-30-2018'); % The start of the first epiweek is April 30, 2018
endDate = datenum('08-12-2018'); % The start date of vaccination is Aug. 12 2018
TVac=endDate-startDateofSim;  % Determine the start time for vaccination to be used in the model

% Melding process
SS=10^3; % The number of sampled wanted to be produced
LL=par(:,end); % Set the log-likelihood values from the parameter fitting
w=exp((LL))./sum(exp((LL))); % Define the weights based on the likelihood (i.e. likelihood is exp(LL))
wc=cumsum(w); % Produce the cumulative sum for we can randomly select the paramter sets       
ind=zeros(SS,1); % The index vector for the samples
r=rand(SS,1); % Generate random numbers for selecting paramters

for ii=1:SS
    gg=find(r(ii)<=wc); % Find the which weights are greater than the random number
    ind(ii)=(gg(1)); % Choose the weight that is first and record the index
end

% Set up the paramters to be used
TTE=par(:,1:48); % Matrix of the attacks included in the model simulation ( there are 48)
par=par(:,49:end); % Truncate the paramter vector
WET=par(:,1:5); % The weights of the different types of events 
% 1- Ville morte
% 2- Helathcare worker attack
% 3- HCW protest
% 4- Attack on ETC
% 5 - other
par=par(:,6:end); % Truncate parameters further
n=par(ind,5); % Value for the hill functiion stauration function
R=par(ind,1); % R effective
Kappa=par(ind,2); % Saturatino value for the saturation function
gamma=par(ind,3); % The rate of removal for infectious class in absence of conflict
epsv=par(ind,6); % The effectiveness of vaccination
mk=par(ind,8); % The rate time to isolation returns to baseline after an attack
k=par(ind,4); % the hyper paramter for the likelihood function
mv=par(ind,7); % % The rate effectiveness of vaccination returns to baseline after an attack
DR=par(ind,9); % Days before attack public health response impacted

% Set up to run only the unique points slected
indu=unique(ind); % Determine the unique index from the sampling
SSu=length(indu); % Number of unique sampled
LL=par(ind,end); % log-likelihood of the sample
mle=find(LL==max(LL)); % the maximum likelihood estimate
mle=mle(1); % The mle
TES=NW*7; % The number of days to run the model for as NW is the number of weeks

% Allocate space for the results
CIdvEq=zeros(SSu,length([0:7:(TES)])); % Cumulative infections

EIdvEq=zeros(SSu,length([0:7:(TES)])); % Cumulative latent infections

CIdvEqR=zeros(SSu,length([0:7:(TES)])); % Cumulative removal/recovery
ReffW=zeros(SSu,length([0:7:(TES)-7])); % Reff for the week
epsW=zeros(SSu,length([0:7:(TES)-7])); % Effectiveness of vaccination for the week
gammaW=zeros(SSu,length([0:7:(TES)-7])); % Time to isoaltino for the week

Avgeps=zeros(SSu,1); % Average values of effectiveness of vaccination
Avggam=zeros(SSu,1); % Avg value for the time to isolation
WW=[0:7:TES]; % Index for the weeks 

% Run through unique samples
 for cl=1:SSu
        ne=2; % Set number latent classes 
        ni=1; % Set number of infectious classes
        a=1/9.4; % Set avg. rate to symptoms for latent infection
        [T,Y,Re,CC,CV]=EpiSimDCVGREff(par(indu(cl),3),par(indu(cl),3).*par(indu(cl),1),par(indu(cl),2),a,f,ne,ni,par(indu(cl),5),TVac,par(indu(cl),6),par(indu(cl),7),TES+7,par(indu(cl),8),TTE(indu(cl),:),WET(indu(cl),:),par(indu(cl),9));
         CI=pchip(T,Y(:,ne+ni+1),[0:7:(TES)]); % Determine cumulative incidence at each week
         EI=pchip(T,Y(:,ne+ni+3),[0:7:(TES)]); % Determine cumulative latent infection at each week
         RepI=pchip(T,Y(:,ne+ni+2),[0:7:(TES)]); % Determine cumulative recovery at each week
         CIdvEq(cl,:)=[CI(1) CI(2:end)-CI(1:end-1)]; % Determine the weekly incidence    (Included the week prior to April 30 which will have cumulative value of zero)  
         EIdvEq(cl,:)=[EI(1) EI(2:end)-EI(1:end-1)]; % Determine the weekly incidence for latent infections    (Included the week prior to April 30 which will have cumulative value of zero) 
         CIdvEqR(cl,:)=[RepI(1) RepI(2:end)-RepI(1:end-1)]; % Determine the weekly incidence for the recovery (Included the week prior to April 30 which will have cumulative value of zero)

         % Calculate weekly averages
         for nn=1:(length(WW)-1)
             ReffW(cl,nn)=sum(pchip(T,Re,[WW(nn):0.01:WW(nn+1)-0.01]).*0.01)./7; % Average Reff for the week

             epsW(cl,nn)=sum(par(indu(cl),6).*pchip(T,CV,[WW(nn):0.01:WW(nn+1)-0.01]).*0.01)./7; % Average effectiveness of vaccination for the week
             gammaW(cl,nn)=sum(pchip(T,(1./par(indu(cl),3)+4.7.*CC),[WW(nn):0.01:WW(nn+1)-0.01]).*0.01)./7; % Average time to isolation (days)
         end  
         Avgeps(cl)=par(indu(cl),6).*sum(pchip(T,CV,[TVac:0.01:(NW*7-0.01)]).*0.01)./(NW*7-TVac); % Average effectiveness of vaccination for the outbreak from the time of vaccination
         Avggam(cl)=(sum(pchip(T,(1./par(indu(cl),3)+4.7.*CC),[TVac:0.01:(NW*7-0.01)]).*0.01))./(NW*7-TVac); % Average time to isolation for the outbreak from the time of vaccination
 end
         
% Transform unique values back to all SS samples         
RID=(zeros(SS,length([0:7:(TES)]))); % Weekly recovery incidence (Has one more index as we have the week before April 30)
CID=(zeros(SS,length([0:7:(TES)]))); % Weekl incidence (Has one more index as we have the week before April 30)
EID=(zeros(SS,length([0:7:(TES)]))); % Weekly incidence ofr latent infections (Has one more index as we have the week before April 30)
RW=(zeros(SS,length([0:7:(TES)-7]))); % Weekly Reff
EW=(zeros(SS,length([0:7:(TES)-7]))); % Weekly effectiveness of vaccination
GW=(zeros(SS,length([0:7:(TES)-7]))); % Weekly time to isolation
AvgepsA=zeros(SS,1); % Avergae effectiveness of vaccination for the outbreak
AvggamA=zeros(SS,1); % Avergae time to isolation for the outbreak

for ii=1:length(ind) % Loop through sampled values
   fg=find(ind(ii)==indu); % Find index (ii) equal to the unique index
   
   CID(ii,:)=CIdvEq(fg,:); % Weekly incidence
   EID(ii,:)=EIdvEq(fg,:); % Weekly incidence latent infection
   RID(ii,:)=CIdvEqR(fg,:); % Weekly incidence recovery/removal
   RW(ii,:)=ReffW(fg,:); % Weekly Reff
   EW(ii,:)=epsW(fg,:); % Weekly effectiveness of vaccination
   GW(ii,:)=gammaW(fg,:); % Weekly time to isoaltion (days)
   AvgepsA(ii)=Avgeps(fg); % Avg. effectivenss of vaccination
   AvggamA(ii)=Avggam(fg); % Avg time to isolation
end

% Compute statsitics
IQRI=prctile(CID,[25 75]); % IQR for weekly incidence
IQRRW=prctile(RW,[25 75]); % IQR for weekly Reff
IQREW=prctile(EW,[25 75]); % IQR for effectiveness of vaccination
IQRGW=prctile(GW,[25 75]); % IQR for time to isolation

UBI=max(CID(find(LL>prctile(LL,5)),:)); % Upper bound for 95% CI for weekly incidence
LBI=min(CID(find(LL>prctile(LL,5)),:)); % Lower bound for 95% CI for weekly incidence
UBE=max(EID(find(LL>prctile(LL,5)),:)); % Upper bound for 95% CI for weekly incidence latnet  infections
LBE=min(EID(find(LL>prctile(LL,5)),:)); % Lower bound for 95% CI for weekly incidence latent infections
UBRW=max(RW(find(LL>prctile(LL,5)),:)); % Upper bound for 95% CI for weekly Reff
LBRW=min(RW(find(LL>prctile(LL,5)),:)); % Lower bound for 95% CI for weekly Reff
UBEW=max(EW(find(LL>prctile(LL,5)),:)); % Upper bound for 95% CI for weekly effectiveness of vaccination
LBEW=min(EW(find(LL>prctile(LL,5)),:)); % Lower bound for 95% CI for weekly effectiveness of vaccination
UBGW=max(GW(find(LL>prctile(LL,5)),:)); % Upper bound for 95% CI for weekly time to isolation
LBGW=min(GW(find(LL>prctile(LL,5)),:)); % Lower bound for 95% CI for weekly time to isolation



meanRI=mean(RID); % mean weekly incidence recovery/removal
medRI=median(RID); % median weekly incidence recovery/removal
mleRI=RID(mle,:); % mle weekly incidence recovery/removal


meanI=mean(CID); % mean weekly incidence
medI=median(CID); % meidan weekly incidence 
mleI=CID(mle,:); % mle weekly incidence

meanE=mean(EID); % mean weekly incidence latent infections
medE=median(EID); % median weekly incidence latent infections
mleE=EID(mle,:); % mle weekly incidence latent infections


meanRW=mean(RW); % mean weekly Reff
medRW=median(RW); % median weekly Reff
mleRW=RW(mle,:); % mle weekly Reff


meanGW=mean(GW); % mean weekly time to isolation
medGW=median(GW); % median weekly time to isolation
mleGW=GW(mle,:); % mle weekly time to isolation


meanEW=mean(EW); % mean weekly effectiveness of vaccination
medEW=median(EW); % median weekly time to isolation
mleEW=EW(mle,:); % mle weekly time to isolation


minEW=min(EW(:,16:end)')'; % Produces the minimum effectivness of vaccination for the outbreak for each sample

TTEP=TTE(ind,:); % The attacks included from the smaple
WETP=WET(ind,:); % The weights of the different events
save(['ModelFit.mat'],'TTEP','mv','Kappa','UBGW','LBGW','UBEW','LBEW','UBRW','LBRW','UBI','LBI','meanI','medI','mleI','meanEW','medEW','mleEW','meanGW','medGW','mleGW','meanRW','medRW','mleRW','epsv','gamma','AvgepsA','AvggamA','LL','mk','IQRI','IQRRW','IQREW','IQRGW','k','CID','RID','n','R','WETP','meanRI','medRI','mleRI','DR','EID','UBE','LBE','meanE','medE','mleE','minEW');
    
    
    
end

