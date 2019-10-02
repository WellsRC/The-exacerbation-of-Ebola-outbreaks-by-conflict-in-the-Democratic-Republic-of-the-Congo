function dx = EpiSimDVGAlt(t,x,g,p,K,a,f,ne,ni,n,TVac,epsv)
% Simulates epidemic and retruns the state of the
% various latent class, infectious classes, and cumulative infections with
% respect to time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t - time
% x - Exposed, infection, cumulative infection
% g - rate of recovery/death
% p - daily rate of infection per infectious case
% K - saturation constant
% a - rate of transitioning from latent to infectious
% f - Saturation function to use
% ne - number of latent phases
% ni - number of infectious phases
% n - value for saturatino function
% TVac - time vaccinatino is started
% epsv - value for the effectiveness of vaccination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dx - rate of change for the different classes
 % 1:ne - latent stage of infection
 % ne+1:ni+ne - infectious stage
 % ne+ni+1 - cumulative infections
 % ne+ni+2 - cumulative recovery
 % ne+ni+3 - cumulative exposed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

epsvt=0; %  Set effectiveness of vaccination to zero as it may not have happened
if(t>=TVac) % if time is greater than or equal to the time of vaccination use epsv
    epsvt=epsv; % adjust the effectiveness of vaccinatino such that vaccination is implemeted
end
dx=zeros(ne+ni+3,1); % define the size of the dx vector to be returned

% Erlang rates
ag=a.*ne; % adjust the rate of transition for the incubation period for the elang distribution
gg=g.*ni; % adjust the rate of recovery/removal for thei nfectious periods for the erlang distribution


% Define index values for call values
E=1;% The starting index for the latent classs
I=ne+1; % The starting index for the infectious class
IC=ne+ni+1; % The index for cumulative infections
RC=ne+ni+2; % The index for cumulative recovery
EC=ne+ni+3; % The index for cumulative latent infections

% Latent period
dx(E)=p*(sum(x(ne+1:ni+ne)))*SatF(f,K,x(IC),t,n).*(1-epsvt)-ag.*x(E); % initial latent period contains the force of infection summing over all infectious classes (i.e. sum(x(ne+1:ni+ne)))
for ii=2:ne
    dx(ii)=ag.*x(ii-1)-ag.*x(ii); % Transition through the other latent classes
end

% Infectious periods
dx(I)=ag.*x(ne)-gg.*x(I); % Last latent class enters first infectious class
for ii=ne+2:ni+ne
    dx(ii)=gg.*x(ii-1)-gg.*x(ii); % Transition through the other infectious  classes
end

%Cumulative infections
dx(IC)=ag.*x(ne); % Count infections onece they are symptomatic (i.e. leave the last latent class and enter first infectious class)
% Cumulative recovery
dx(RC)=gg.*x(ni+ne); % Count recovery onece they are removed from the infectious class (i.e. leave the last infectious class)
% Cumulative latent 
dx(EC)=p*(sum(x(ne+1:ni+ne)))*SatF(f,K,x(IC),t,n).*(1-epsvt); % Count new latent infections (i.e. enter the first latent class)
end

