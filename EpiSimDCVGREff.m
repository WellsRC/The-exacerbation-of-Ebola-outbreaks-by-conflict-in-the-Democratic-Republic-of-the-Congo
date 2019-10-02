function [T,Y,RE,CC,CV] = EpiSimDCVGREff(g,pn,K,a,f,ne,ni,n,TVac,epsv,mV,TE,mTI,TEE,WET,DR) 
% Used to provide esitmates from the simulation
%===============================
% Input
%===============================
%g - the rate of isolation in the absence of conflict
%pn the infection rate
%K - the saturation constant 
%a - the rate a latent infected case becomes infectious
%f - The saturnation fucntion to be used in the simulation
%ne - the number of latent satges for the disease
%ni - the number of infectious stages for the disease
%n - the hill couefficent for the one stauration function
%TVac - the time that vaccination is started
%epsv - the efficacy of the vaccine
%mV - rate time to vaccination returns to baseline after attack
%TE - The time to run the simulation to
%mTI - rate time to isolation returns to baseline after attack
%mk - the long-term effect of attacks 
%TEE- the attacks to be used in the simuation
% daily conflict
%WET - The weights of the different types of attacks
%DR- days for ramp-up before attack
%===============================
% Output
%===============================
% T - the vector of time
% Y - A matrix size ne+ni+3 x length(T) of the dynamics of the outbreak
 % 1:ne - latent stage of infection
 % ne+1:ni+ne - infectious stage
 % ne+ni+1 - cumulative infections
 % ne+ni+2 - cumulative recovery
 % ne+ni+3 - cumulative exposed
% CC - % Effect of conflit on the time to sioalation
% CV - % Effect of conflit on vaccination

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Run 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=0.01; % Specify the time step
% Initialize the simulation
Y=zeros(ne+ni+3,TE./dt+1); % Output from system of equations
RE=zeros(1,TE./dt+1); % R_eff for the course of the outbreak
CC=zeros(1,TE./dt+1); % The extent conflict increases the time to isolation
CV=zeros(1,TE./dt+1); % The extent effectiveness of vaccination is reduced by conflict

Y(1,1)=1; % Start with a single latent infected individual
T=zeros(1,TE./dt+1); % Time vector


nn=1; % Starting index
while((nn-1)*dt<TE)   % run while less than time desired (nn-1) as we are starting at t=0 
    
    % Calculate the effects of attacks  on the time to isolation
    CeTI=TCDC((nn-1)*dt,mTI,TEE,0,WET,DR);
    % Calculate the effects of attacks  on effectiveness of vaccination
    CeEV=TCDC((nn-1)*dt,mV,TEE,1,WET,DR);


    %The impact on the time to isolation
    CF=(1-CeTI);
    % The impact on the effectiveness of vaccination
    CF2=(1-CeEV);

    
    if((nn-1)*dt<TVac)% Time before vaccination
        RE(nn)=(pn)./(1/(1/g+4.7*(1-CF)))*SatF(f,K,Y(end,nn),(nn-1)*dt,n); % Calculation of Reff
        CV(nn)=0; % Effect of conflict on vaccination has no effect since vaccinatino has yet to occur
    else                       % time after vaccination
        RE(nn)=(pn)./(1/(1/g+4.7*(1-CF)))*SatF(f,K,Y(end,nn),(nn-1)*dt,n)*(1-epsv*CF2); % Calculation of Reff
        CV(nn)=CF2; % Effect of conflict on vaccination
    end
    CC(nn)=(1-CF); % Effect of conflit on the time to sioalation
    
    % Upated the status of the system
    Y(:,nn+1)=Y(:,nn)+dt.*EpiSimDVGAlt((nn-1)*dt,Y(:,nn),1/(1/g+4.7*(1-CF)),pn,K,a,f,ne,ni,n,TVac,epsv*CF2);                
    T(nn+1)=nn*dt; % update the time vector
    nn=nn+1; % increase the index               
end 

% Calcuate for the last time point

% Calculate the effects of attacks  on the time to isolation
CeTI=TCDC((nn-1)*dt,mTI,TEE,0,WET,DR);
% Calculate the effects of attacks  on effectiveness of vaccination
CeEV=TCDC((nn-1)*dt,mV,TEE,1,WET,DR);

%The impact on the time to isolation
CF=(1-CeTI);
% The impact on the effectiveness of vaccination
CF2=(1-CeEV);

RE(nn)=(pn)./(1/(1/g+4.7*(1-CF)))*SatF(f,K,Y(end,nn),(nn-1)*dt,n)*(1-epsv*CF2); % Calculation of Reff
CV(nn)=(1-CF);       % Effect of conflit on the time to sioalation    
CC(nn)=CF2;  %Effect of conflict on vaccination

Y=Y';
end
