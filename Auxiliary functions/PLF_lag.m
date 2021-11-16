function [plf,lag] = PLF_lag(S1,S2)
% Summary of this function goes here
% -------------------------------------------------------------------------
%   [plf,lag] = PLF_lag(S1,S2)
%--------------------------------------------------------------------------
% FUNCTION
%   Extract functional Connectivity (phase synchrony) b/w 2 signals 
%   using Phase Locking Factor Synchrony model (pure phase synchrony) 
%   sometimes it is called Phase Locking Value 
% -------------------------------------------------------------------------
% INPUTS:
%   S1, S2 - two 1D signals [1 x time points]
%
% OUTPUT:
%   plf - phase locking factor/value UNdirected index
%   lag - lag between the two signals
% -------------------------------------------------------------------------
% REFERENCES:
%  Lachaux et al., (1999); 
%  Petkov, G., Goodfellow, M., Richardson, M. P., & Terry, J. R. (2014). 
%    A critical role for network structure in seizure onset: 
%      a computational modeling approach. Front Neurol, 5, 261. 
%      http://doi.org/10.3389/fneur.2014.00261
% -------------------------------------------------------------------------
% VERSION
% GPetkov, 02.06.2013
% (K&P lab)
% Modified - M.Lopes @ 2017
plf = [];
    if nargin<2
        disp( 'usage: plf = fg_PLF_13( S1, S2 );' );             
        return;     
    end
    Cm = angle(hilbert(S1'))-angle(hilbert(S2'));    
    arg = mean(exp(1i*Cm));    
    plf = abs(arg);
    lag = angle(arg);     
return
end
