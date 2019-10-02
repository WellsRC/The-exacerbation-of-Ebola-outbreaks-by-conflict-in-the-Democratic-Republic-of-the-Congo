function BM(IData)
%  Run the Bayesian melding for the conflict model for a primary search of
%  the paramter space (Run BMFocus later)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IData- Data to be used in the fitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle');
for ns=1:100
    nn=1;
    f=3; % Specify which saturatino functino to be used
    % Create file to save to
    while (exist(['BM_M' num2str(f) '_' num2str(nn) '.mat'],'file')==2)
         nn=nn+1;
    end
    save(['BM_M' num2str(f) '_' num2str(nn) '.mat']);
    
    % Melding sampling
    SS=5.*10^4; % Number of samples to run

    TT=[1:length(IData)]; % The weeks we want to fit the model to
    L2=zeros(SS,1); % Allocate space for log-likelihood
    lhs=lhsdesign(SS,8); % Use latin-hypercube sampling
    ne=2; % Define the number of latent phases
    ni=1; % Define the numebr of infectious phases
    a=1./9.4; % Avg rate of latent period 
    R=1.3+(2-1.3).*lhs(:,1); % Sample Reff
    K=10.^(log10(2*10^(-4))+(log10(8*10^(-4))-log10(2*10^(-4))).*lhs(:,2));  % Sample saturatino constant
    k=20+30.*lhs(:,3); % Sample hyper parameter
    g=1/(4.69+2.63)+(1/(4.69-2.63)-1/(4.69+2.63)).*lhs(:,4); % Sample rate of removal/recovery in absence of conflic
    n=1+(2.5-1).*lhs(:,5);  % Hill coefficient
    mv=0.1+(1.4-0.1).*lhs(:,6);  % Rate effectiveness of vaccinatino returns to baseline
    epsv=0.15+0.5.*lhs(:,7);  % Effectienss of vaccination
    mk=0.1+(3-0.1).*lhs(:,8);  % Rate time to isoaltion returns to baseline
    DR=7+(randi(22,SS,1)-1);  % Duration before attack public health response impeded
    % Vector for x
%         R=x(1); %daily rate of infection per infectious case
%         K=x(2); %saturation constant 
%         g=x(3); % duration of infection
%         k=x(4); % the hyper paramter for the likelihood function
%         n=x(5); % The hill coefficnet
%         epsv=x(6); % the effectiveness of vaccination
%         mV=x(7); % the extendend effects of an attack on the effectiveness of vaccination
%         mTI=x(8); % the extendend effects of an attack on the the time to isolation
%         DR = x(9); % THe duration it take for the hostility to grow to
%         100%
    x=[R K g k n epsv mv mk DR];
    
    TTE=round(rand(SS,48)); % Randomly include attacks and round to obtian 0,1  

    WET=0.6+0.4.*rand(SS,5); % The weights for the 5 types of attacks

    % Compute the log-likelihood for the samples
    for jj=1:5
        ts=1+10^4*(jj-1);
        te=10^4*jj;
        parfor ii=ts:te   
            L2(ii)=-FitEpiSimNK(x(ii,:),f,a,ne,ni,IData,TT,TTE(ii,:),WET(ii,:));   % Calculates the likelihood for the paramterset (the negative is there becase it returns the negative log likelihood)
        end
        save(['BM_M' num2str(f) '_' num2str(nn) '.mat']); % Save the output
    end
end
end