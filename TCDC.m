function Ce = TCDC(t,W,TE,MechI,WTE,DB)
% Returns the value of the conflict function
%===============================
% Input
%===============================
% t - The time during the outbreak 
% W - Rate of decay in the outbreak
% TE - A binary vector specifying which events will be used in the course
% of the simulation
% Mech I is the mechanism conflict is impacting 
%       MechI=0 is for time to isolation
%       MechI=1 is effectivenes sof vaccination
% WTE is the weight based on the type of event happenening
% DB is the number of days before an attack that we would see an effect on
% disease control efforts
%=================================
% Output
%=================================
% Ce - The value of the conflict function due to the attacks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hardcoded the attacks to be used
% TCtemp
% Column 1 (TCtemp(:,1)) - Day of the attack (i.e t=0 April 30)
% Column 2 (TCtemp(:,2)) - Duratino of hightened effect (We did not use this in the model)
% Column 3 (TCtemp(:,3)) - Location of the attack
        % 1- Mabalako
        % 2- Mandima
        % 3 - Beni
        % 4 - Butembo and Katwa
        % 5 - Musienene
        % 6- Other
% Column 4 (TCtemp(:,4)) - Type of the attack    
        % 1- Ville morte
        % 2- Attack Healthcare worker
        % 3 - Health care worker protest
        % 4 - Attack ETC
        % 5 - Other
TCtemp=[147,5,3,1;155,1,4,2;173,1,4,2;173,2,3,1;190,1,4,3;241,4,3,4;241,4,3,4;241,4,3,4;241,4,6,4;241,4,3,4;256,1,6,2;269,1,6,2;281,1,6,4;287,1,4,4;291,1,6,5;294,1,6,3;300,1,6,4;300,1,4,4;303,1,4,4;312,1,3,4;313,1,4,4;318,1,6,4;326,1,4,4;328,1,4,2;329,1,6,3;340,1,3,2;350,1,4,5;353,1,4,5;354,1,4,4;355,1,4,4;359,1,4,3;365,1,6,5;368,1,4,2;372,1,4,4;372,1,6,2;377,1,4,4;378,1,4,4;379,1,6,4;390,1,3,4;391,1,4,4;401,1,3,2;407,3,6,1;409,3,6,1;409,2,3,1;412,1,6,4;414,1,3,4;420,1,3,2;421,1,6,2];

% Find the events to be included which are indicated by values of 1
f=find(TE==1);
TC=TCtemp(f,:); % Restructure the Matric of TCtemp to be of the ones oncluded

% Weights for the location of the attacks
if(MechI==0) % The impact on the time to isolation
    WADC=[1	0.4390755702	0.54941052	0.6142427696	0.8560627894	0.6727183894]; %Index are for attacks in: Mabalako Mandima Beni Butembo&Katwa Musienene Other
else % the impact on the effectiveness of vaccination
    WADC=[0.6411571435	0.6884600079	0.8314556328	1	0.7593736245	0.4630178471]; %Index are for attacks in: Mabalako Mandima Beni Butembo&Katwa Musienene Other
end

%Set the growth of effect before the attack to increase linear
NR=1;

% NOT Accounting for how long the ramifcations of the attack took place (i.e. DO NOT consider that a ville morte that lasted 5 days) 
NP = TriDist(t,TC(:,1),DB,W,ones(size(TC(:,1))),NR); % Set Duratino of the attacks to be one as could not identify consistent informatino for all

Ce=1-prod(1-WTE(TC(:,4))'.*WADC(TC(:,3))'.*NP); % Compund the effects of the attacks WTE(TC(:,4))' is the weight based on type of attack WADC(TC(:,3))' is weigth based on the location and NP is the impact of the attack relative to time t     

end

