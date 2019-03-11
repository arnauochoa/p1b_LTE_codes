%SPCOMNAV - J.A. Del Peral (12/03/11)
%==========================================================================
% GENERATION OF PRIMARY SYNCHRONIZATION SIGNALS
%
% [d] = genPSS(nID2)
%
% nID2 [0,1,2]      Physical-layer identity within the cell-ID group       
%
%==========================================================================
function [d] = genPSS(nID2)

% Root indices for the primary synchronization signal

switch nID2
    case 0, u = 25;
    case 1, u = 29;
    case 2, u = 34;
end

FORM = 2;    % [1 = Standard, 2 = LTE book]

% Frequency-domain Zadoff-Chu squences
switch FORM
    case 1  % STANDARD
        d = zeros(1,62);
        n = 0:30;
        d(n+1) = exp(-j*pi*u*n.*(n+1)/63);
        n = 31:61;
        d(n+1) = exp(-j*pi*u*(n+1).*(n+2)/63);
    case 2  % LTE BOOK
        d = zeros(1,63);
        n = 0:62;
        d(n+1) = exp(-j*pi*u*n.*(n+1)/63);
        d = [d(1:31) d(33:end)];
end

% *with L = 62, the result is the same for both procedures 