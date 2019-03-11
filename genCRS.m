%SPCOMNAV - J.A. Del Peral (14/03/11)
%==========================================================================
% GENERATION OF CELL-SPECIFIC REFERENCE SIGNALS
%
% [d,pos] = genCRS(nSlot,nSymbol,nCP,nIDc,nRB)
%
% nSlot   [0,1]       Slot number within a radio frame
% nSymbol [0,..6]     OFDM symbol number within the slot
% nCP     [0,1]       Normal CP = 1, extemded CP = 0
% nIDc    [0,...,503] Cell-ID
% nRB     [1,...,110] Maximum number of resource blocks
%   
%==========================================================================
function [d,pos] = genCRS(nSlot,nSymbol,nCP,nIDc,nRB)

% %Default
% nSlot = 1;
% nSymbol = 1;
% 
% nIDc = 0;
% nCP  = 1;
% nRB = 31;

Nsymb = (nCP==1)*7 + (nCP==0)*6; % OFDM symbols per slot

%==========================================================================
%% Generate sequence
%==========================================================================

% Pseudo-random sequence according to Section 7.2 in TS36.211
N  = 31;    % Pseudo-random sequences are defined by a length-31 Gold sequence
Nc = 1600;  % The first 1600 samples are discarded, fast-forward

nRBmax = 110; % Maximum # of RB

c_init = 2^10 * (7*(nSlot + 1) + nSymbol + 1)*(2*nIDc + 1) + 2*nIDc + nCP;

% Initialization___________________________________________________________

x1 = [1 zeros(1,N-1)];

% Algorithm from DEC2BIN function (without char result)
[f,e]=log2(max(c_init));
x2=rem(floor(c_init*pow2(0:-1:(1-max(N,e)))),2);
%c_init==sum(x2.*2.^(0:30)) % CHECKING

% The two m-sequences x1(n) and x2(n) are respectively generated by
% feedback polynomials D31 +D3 +1 and D31 +D3 + D2 + D + 1

for n = 1:(4*nRBmax+Nc)
    x1(n+N) = mod(x1(n + 3) + x1(n),2);
    x2(n+N) = mod(x2(n + 3) + x2(n + 2) + x2(n + 1) + x2(n),2);
end

c= zeros(1,4*nRBmax);
for n = 1:4*nRBmax
    c(n) = mod(x1(n + Nc) + x2(n + Nc),2);
end

% Reference-signal sequence generation_____________________________________

d = zeros(1,2*nRB);
k = 0:2*nRB-1;
m = k + nRBmax - nRB;
d(k+1) = 1/sqrt(2)*(1-2*c(2*m+1)) + 1i*1/sqrt(2)*(1-2*c(2*m+2));
% length(r);
% plot(cxcorr(r,r))

%==========================================================================
%% Mapping
%==========================================================================

p = 0; % Antenna port 0, 1, 2, 3
if p < 2, l = [0 Nsymb-3];
else l = 1; end

m = 0:2*nRB-1;
if sum(nSymbol == l)
    switch p
        case 0, v = 0*(nSymbol==0) + 3*(nSymbol~=0);
        case 1, v = 3*(nSymbol==0) + 0*(nSymbol~=0);
        case 2, v = 3*mod(nSlot,2);
        case 3, v = 3 + 3*mod(nSlot,2);
    end
    v_shift = mod(nIDc,6);
    pos = 6*m + mod(v + v_shift,6) + 1; % +1 for MATLAB indexing
    pos(nRB+1:end) = pos(nRB+1:end)+1; % Avoid DC subcarrier
else
    pos = [];
end
