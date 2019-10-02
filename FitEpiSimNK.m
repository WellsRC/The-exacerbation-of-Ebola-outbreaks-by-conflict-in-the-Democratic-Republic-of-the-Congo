function F = FitEpiSimNK(x,f,a,ne,ni,IData,TT,TTE,WET)
% Fits the model using daily conflict data and dates of attacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x - vector for the parameters
%         R=x(1); %Basic reproductive value
%         K=x(2); %saturation constant 
%         g=x(3); % duration of infection
%         pn=R*g; % infection rate
%         k=x(4); % the hyper paramter for the likelihood function
%         n=x(5); % The hill coefficnet
%         epsv=x(6); % the effectiveness of vaccination
% 
%         mV=x(7); % the extendend effects of an attack on the effectiveness of vaccination
%         mTI=x(8); % the extendend effects of an attack on the the time to isolation
%         DB=x(9); % THe time it takes for the hostility to reach 100% before the attack
% f - the stauration function to be used
% a - rate of transitioning from latent to infectious
% ne - the number of exposed classes in the model
% ni - the number of infectious classes in the model
% IData - Weekly incidence data 
% TT - the Weeks that the are being used to fit the model (In case want to
% use partial data or evaluate a certain point)
% TTE - the attacks to be included in the model
% WET - the weights based on the type of attack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F - the negative log-likelihood

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define paramters
R=x(1); %Basic reproductive value
K=x(2); %saturation constant 
g=x(3); % duration of infection
pn=R*g; % infection rate
k=x(4); % the hyper paramter for the likelihood function
n=x(5); % The hill coefficnet
epsv=x(6); % the effectiveness of vaccination

mV=x(7); % the extendend effects of an attack on the effectiveness of vaccination
mTI=x(8); % the extendend effects of an attack on the the time to isolation
DB=x(9); % THe time it takes for the hostility to reach 100% before the attack
    
startDateofSim = datenum('04-30-2018'); % Start date of the simulation
endDate = datenum('08-12-2018'); % THe date that vaccinayion started
TVac=endDate-startDateofSim; % Comput the time that vaccinatino started
TE=1+7.*length(IData); % Time to run the simulation (IData is in weeks so need to multiply by 7. Also ran for one extra day

% Run the model to obtain output
[T,Y] = EpiSimNKV(g,pn,K,a,f,ne,ni,n,TVac,epsv,mV,TE,mTI,TTE,WET,DB);

CI=pchip(T,Y(:,ne+ni+1),[0:7:TE]); % Determine cumulative incidence at each week
CI(2:end)=CI(2:end)-CI(1:end-1); % Calculate the weekly incidence (NOTE: CI(1)=0 and is the week before April 30)
        
I=[IData];

EI=[1;2;2;1;3;4;4;4;4;4;2;2;4;4;5;3;4;4;4;4;4;4;4;2;3;2;5;5;4;4;3;5;2;4;2;4;4;4;4;4;2;1;1;1;4;4;2;2]; % The type of attack
TC=TT+1; % adding one since the cumulative infection is zero at the first point (Week prior to April 30) and Incidence for April 30 is at index 2 for CI
pnb=k./(k+CI(TC)'); % Compute value for negative binomial
PP=nbinpdf(I(TT),k,pnb); % Evalaute pdf for negative binomial at the time points of interest
PP(isnan(PP))=0; % Ensure there are no Nan value and if so set to zero
F=-sum(log(PP))-length(TTE).*log(sum(TTE.*WET(EI))./length(TTE)); % The negative is is we want to minimizing if we ever want to use fmin con (Data log likelihood sum(log(PP)) and maximize impact of estimated attacks length(TTE).*log(sum(TTE.*WET(EI))./length(TTE))
end

