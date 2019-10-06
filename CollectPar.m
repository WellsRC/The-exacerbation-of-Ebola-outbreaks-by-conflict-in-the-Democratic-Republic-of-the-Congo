% CollectPar collects the parameters from the sampling and puts them into
% one file
            

f=3; % Saturation function one is using
% f- Saturatino function to use
%   1 - Exponential decay w.r.t. cumulative incidence
%   2 - Exponential decay w.r.t. time
%   3 - Hill function saturating w.r.t. cumulative incidence
%   4 - Hill function saturating w.r.t. time
%   5 - Discoutn function saturating w.r.t. cumulative incidence
%   6 - Discoutn function saturating w.r.t. time

par=[]; % Specify an empty parameter set to add to (otherwise use past collected paramaters)

% Code to load past file if desired

%load(['ModelPar_' num2str(f) '.mat'],'par'); % Loads past collected parameter file if waneted
%par=par(par(:,end)~=0,:); % Ensuring that there are no paramter sets added
%where the log-lilelihood was not calculated
%par=unique(par,'rows'); % Ensures that the parameter sets are unique
%par=flip(sortrows(par,length(par(1,:)))); % Need to flip the matrix to ensure the best are at the top for truncation
%par=par(1:min([2.5.*10^6 length(par(:,1))]),:); % Truncate the paramater sets if they exceed 2.5.*10^6

nfs=1000; % Number of files to run
NFR=floor(nfs/50); % Paritions the files into loop size of 50. Found that ran into memory issues otherwise

% Run loop to collect parameters
for mm=1:NFR
    for nn=1+50*(mm-1):50*mm   
        load(['BM_M' num2str(f) '_' num2str(nn) '.mat']); % load(['BMFocus_M' num2str(f) '_' num2str(nn) '.mat']);                               
        gg=find(L2~=-Inf); % Find the ones that are not negative infinity
        hh=find(L2(gg)~=0); % Find ones that are non-zero (in case the melding process was stopped prematurley)
        gg=gg(hh); % Update the index for desired paramter sets
        L=L2(gg); % specify the log-likelihood                                           
        par=[par;TTE(gg,:) WET(gg,:) R(gg) K(gg) g(gg) k(gg) n(gg) epsv(gg) mv(gg) mk(gg) DR(gg) L]; % Save the parameters to the file with the likelihood at the end
    end
    par=unique(par,'rows'); % Ensures that the parameter sets are unique
    par=flip(sortrows(par,length(par(1,:))));  % Need to flip the matrix to ensure the best are at the top for truncation
    par=par(1:min([2.5.*10^6 length(par(:,1))]),:); % Truncate the paramater sets if they exceed 2.5.*10^6
end

% Run remianing files to add
for nn=(NFR*50)+1:nfs   
    load(['BM_M' num2str(f) '_' num2str(nn) '.mat']); % load(['BMFocus_M' num2str(f) '_' num2str(nn) '.mat']);                               
    gg=find(L2~=-Inf); % Find the ones that are not negative infinity
    hh=find(L2(gg)~=0); % Find ones that are non-zero (in case the melding process was stopped prematurley)
    gg=gg(hh); % Update the index for desired paramter sets
    L=L2(gg); % specify the log-likelihood                                           
    par=[par;TTE(gg,:) WET(gg,:) R(gg) K(gg) g(gg) k(gg) n(gg) epsv(gg) mv(gg) mk(gg) DR(gg) L]; % Save the parameters to the file with the likelihood at the end
end
par=unique(par,'rows'); % Ensures that the parameter sets are unique
par=flip(sortrows(par,length(par(1,:))));  % Need to flip the matrix to ensure the best are at the top for truncation
par=par(1:min([2.5.*10^6 length(par(:,1))]),:); % Truncate the paramater sets if they exceed 2.5.*10^6

save(['NK_M' num2str(f) '.mat'],'par'); % Save the paramter file