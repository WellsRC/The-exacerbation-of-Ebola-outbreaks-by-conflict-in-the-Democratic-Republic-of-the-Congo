function Y = TriDist(t,tE,DB,W,DE,NR)
% The effects of an attack(s) relative to time t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t- time
% tE - time of events
% DB - day before the attack that it has impact on response
% W - The decay rate of the event
% DE - Duratino of the event such that it influences decay rate (We assumed DE =1 as duratino was not well defined for all events)
% NR - polynomial value for growth (We assumed NR=1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y - imapct of each attack relative to time t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y=zeros(size(tE));  % Size of the Y vector. set to zero

if(DB>0) % Days before must be greater than zero to compute
    Y(t<tE)=((t-tE(t<tE))./DB+1); % The linear growth rate for any event not yet to happen (i.e. t<tE) increases from zero to 1
    Y(Y<0)=0; % Ensure that any negative is zeor
    Y=Y.^NR; % Introduce any non-linearaity
end

% Decay affects
Y(t>=tE)=exp(-W.*(t-tE(t>=tE))./DE(t>=tE)); % Rate of effect decays exponentially after for attacks that have happened (i.e. t>=tE)
end
