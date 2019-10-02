function BMFocus(IData)
%  Run the Bayesian melding for the conflict model focused around the
%  better smapled vlaues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IData- Data to be used in the fitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           
rng('shuffle');
for ns=1:100
    f=3; % The saturatino function to be used
    load(['NK_M' num2str(f) '.mat']); % load paramter sets already fit to
    LL=par(:,end); % Set the log-likelihood
    w=exp((LL))./sum(exp((LL))); % Calculate the weights of the likelihood (i.e. exp((LL)))
    wc=cumsum(w);  % Calculate the cumulative wieght for the selection of paramters

    r=rand(1); % Generate a random number
    gg=find(r<=wc); % Find set above randomly slected number
    ind=(gg(1)); % Take the first instance
    par=par(ind,49:end); % truncate the paramters (i.e. remove the attacks included) and to the one selected
    WWP=par(:,1:5); % The weight
    par=par(:,6:end); % truncate the remaining set of par 
    nn=1; 
    % Create file to save to
    while (exist(['BMFocus_M' num2str(f) '_' num2str(nn) '.mat'],'file')==2)
         nn=nn+1;
    end
    save(['BMFocus_M' num2str(f) '_' num2str(nn) '.mat']);
    % Melding sampling
    SS=5.*10^4; % Number of samples to run

    TT=[1:length(IData)]; % The weeks we want to fit the model to
    L2=zeros(SS,1); % Allocate space for log-likelihood
    lhs=lhsdesign(SS,8); % Use latin-hypercube sampling
    ne=2; % Define the number of latent phases
    ni=1; % Define the numebr of infectious phases
    a=1./9.4; % Avg rate of latent period 
    % allocate space for the parameters 
    R=zeros(size(lhs(:,1)));
    K=zeros(size(lhs(:,1)));
    k=zeros(size(lhs(:,1)));
    g=zeros(size(lhs(:,1)));
    n=zeros(size(lhs(:,1)));
    mv=zeros(size(lhs(:,1)));
    epsv=zeros(size(lhs(:,1)));
    mk=zeros(size(lhs(:,1)));
    DR=zeros(size(lhs(:,1)));
    for mm=1:5
        R([1:10^4]+(mm-1)*10^4)=par(1).*(1+0.3.*(0.5-lhs([1:10^4]+(mm-1)*10^4,1))); % Sampled +/- %15 around R_0
        K([1:10^4]+(mm-1)*10^4)=10.^(log10(par(2)).*(1+0.3.*(0.5-lhs([1:10^4]+(mm-1)*10^4,2)))); % Sampled +/- %15 around saturation constant
        k([1:10^4]+(mm-1)*10^4)=par(4).*(1+0.3.*(0.5-lhs([1:10^4]+(mm-1)*10^4,4))); % Sampled +/- %15 around hyperparamter
        g([1:10^4]+(mm-1)*10^4)=par(3).*(1+0.3.*(0.5-lhs([1:10^4]+(mm-1)*10^4,3))); % Sampled +/- %15 around rate to recovery/removal for infection
        n([1:10^4]+(mm-1)*10^4)=par(5).*(1+0.3.*(0.5-lhs([1:10^4]+(mm-1)*10^4,5))); % Sampled +/- %15 around hill coefficeint
        epsv([1:10^4]+(mm-1)*10^4)= par(6).*(1+0.3.*(0.5-lhs([1:10^4]+(mm-1)*10^4,6))); % Sampled +/- %15 around effectiveness of vaccination
        mv([1:10^4]+(mm-1)*10^4)= par(7).*(1+0.3.*(0.5-lhs([1:10^4]+(mm-1)*10^4,7))); % Sampled +/- %15 around rate effectiveness of vaccination returns to baseline
        mk([1:10^4]+(mm-1)*10^4)= par(8).*(1+0.3.*(0.5-lhs([1:10^4]+(mm-1)*10^4,8))); % Sampled +/- %15 around rate time to isolation retirns to baseline
        DR([1:10^4]+(mm-1)*10^4)= par(9)+(5 -randi(9,10^4,1)); % Sampled +/- 4 days around time before attask
        for jj=1:5
            WET([1:10^4]+(mm-1)*10^4,jj)=WWP(jj).*(1+0.3.*(0.5-rand(10^4,1))); % Sampled +/- %15 around the weights for the different attacks
        end
    end  
    % Ensure weights between 0 and 1
    WET(WET>1)=1; 
    WET(WET<0)=0;
    % Ensure R_0>=1
    R(R<1)=1;
    % Ensure hyperparamter not fall below calibrated lower bound
    k(k<20)=20;
    % Ensure time to isolation in range
    g(g<1/(4.69+2.63))=1/(4.69+2.63);
    g(g>1/(4.69-2.63))=1/(4.69-2.63);
    % ensure effectiveness of vaccinatino between 0 ans 1
    epsv(epsv<0)=0;
    epsv(epsv>1)=1;
    % Ensure rate return to baseline is non-negative
    mv(mv<0)=0;
    % Ensure rate return to baseline is below calibrated upperbound
    mv(mv>3)=3;
    % Ensure rate return to baseline is non-negative
    mk(mk<0)=0;
    % Ensure rate return to baseline is below calibrated upperbound
    mk(mk>3)=3;
    % Ensure day before not below lower bound
    DR(DR<7)=7;
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
    
    TTE=round(rand(SS,48)); % Randomly select attacks to include 
   
    
    % Compute the log-likelihood for the samples
    for jj=1:5
        ts=1+10^4*(jj-1);
        te=10^4*jj;
        parfor ii=ts:te   
            L2(ii)=-FitEpiSimNK(x(ii,:),f,a,ne,ni,IData,TT,TTE(ii,:),WET(ii,:));   % Calculates the likelihood for the paramterset 
        end

    save(['BMFocus_M' num2str(f) '_' num2str(nn) '.mat'],'R','K','g','k','n','epsv','mv','mk','DR','TTE','WET','L2');
    end
end
end


