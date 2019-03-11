%SPCOMNAV - J.A. Del Peral (24/03/16)
%==========================================================================
% Provides symbol and subcarrier allocation of the LTE pilot signals.
%   - Pilots supported: PSS, SSS and CRS.
%
% Inputs:
%       nIDc            Cell ID
%       nRB             Number of resource blocks       
%       N               Number of subcarriers
%       N_symb_frame    Number of symbols per radio frame
%       nCP             CP configuration
%       ind_SS          Indexes for mapping of synchronization signals
%
% Outputs:
%       d_S             pilot sequence for every symbol and subcarrier       
%       pilot_signal    type of pilot
%
%==========================================================================
function [d_S,pilot_signal] = sdr_LTE_2Dpilots(nIDc,nRB,N,N_symb_frame,nCP,ind_SS)

nID1 = floor(nIDc/3);
nID2 = nIDc - 3*nID1;

d_S = zeros(N,N_symb_frame); pilot_signal = d_S;
for nsymb = 1:N_symb_frame
    nSlot = mod(floor((nsymb-1)/7),20);
    SF = mod(nsymb,N_symb_frame) > 71;
    % Synchronisation signals
    if mod(nSlot,10) == 0
        if mod(nsymb-1,7) == 5  % SSS
            d_S(ind_SS,nsymb) = genSSS(nID1,nID2,SF);
            pilot_signal(ind_SS,nsymb) = 2;
        elseif mod(nsymb-1,7) == 6  % PSS
            d_S(ind_SS,nsymb) = genPSS(nID2);
            pilot_signal(ind_SS,nsymb) = 1;
        end
    end
    % Cell-specific reference signals
    [pilots,pos] = genCRS(nSlot,mod(nsymb-1,7),nCP,nIDc,nRB);
    if not(isempty(pos))
        ind_CRS = pos + (N-nRB*12)/2; 
        d_S(ind_CRS,nsymb) = pilots;
        pilot_signal(ind_CRS,nsymb) = 3;
    end
end