function [T,Y] = EpiSimNKV(g,pn,K,a,f,ne,ni,n,TVac,epsv,mV,TE,mTI,TEE,WET,DB)     
% Uses the forward Euler method to simulate the an outbreak with conflict
% for the fitting process
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

% Set the step size
dt=0.1;
% Initialize the outbreak
Y=zeros(ne+ni+1,TE./dt+1); % Output from system of equations
Y(1,1)=1; % Start with a single latent infected individual
T=zeros(1,TE./dt+1); % Time vector
%Set index for counting to one
nn=1;
while((nn-1)*dt<TE)    
    % Calculate the effects of attacks  on the time to isolation
    CeTI=TCDC((nn-1)*dt,mTI,TEE,0,WET,DB);
    % Calculate the effects of attacks  on effectiveness of vaccination
    CeEV=TCDC((nn-1)*dt,mV,TEE,1,WET,DB);
    % % Calculate the effects of daily conflict
    %Cg = muC((nn-1)*dt);

    %The impact on the time to isolation
    CF=(1-CeTI);
    % The impact on the effectiveness of vaccination
    CF2=(1-CeEV);
 
    
    
    %Update the Y vector
    Y(:,nn+1)=Y(:,nn)+dt.*EpiSimDVG((nn-1)*dt,Y(:,nn),1/(1/g+4.7*(1-CF)),pn,K,a,f,ne,ni,n,TVac,epsv*CF2);                
    % Update the time vector
    T(nn+1)=nn*dt;
    % increase the time step
    nn=nn+1;                
end

Y=Y';
end