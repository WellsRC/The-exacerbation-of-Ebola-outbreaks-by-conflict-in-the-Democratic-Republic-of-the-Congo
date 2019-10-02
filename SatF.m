function Y = SatF(f,K,CI,t,n)
% Returns the value of the saturatrion function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% f- Saturatino function to use
%   1 - Exponential decay w.r.t. cumulative incidence
%   2 - Exponential decay w.r.t. time
%   3 - Hill function saturating w.r.t. cumulative incidence
%   4 - Hill function saturating w.r.t. time
%   5 - Discoutn function saturating w.r.t. cumulative incidence
%   6 - Discoutn function saturating w.r.t. time
% K - saturatino constant
% CI - umulative incidence
% t - time in outbreak
% n - Hill coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Y - Value of saturatino function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch f
    case 1 %   1 - Exponential decay w.r.t. cumulative incidence
        Y=exp(-K*(max(CI,0))); % Max is used incase of numerical error early on in a simulation using ode45
    case 2 %   2 - Exponential decay w.r.t. time
        Y=exp(-K*t);
    case 3 %   3 - Hill function saturating w.r.t. cumulative incidence
        Y=1./(1+(K*(max(CI,0)))^n); % Max is used incase of numerical error early on in a simulation using
    case 4 %   4 - Hill function saturating w.r.t. time
        Y=1./(1+(K*t)^n);
    case 5 %   5 - Discoutn function saturating w.r.t. cumulative incidence
        Y=1./(1+K)^max(CI,0); % Max is used incase of numerical error early on in a simulation using
    case 6 %   6 - Discoutn function saturating w.r.t. time
        Y=1./(1+K)^t;
end

end

