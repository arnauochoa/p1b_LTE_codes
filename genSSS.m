%SPCOMNAV - J.A. Del Peral (11/03/11)
%==========================================================================
% GENERATION OF SECONDARY SYNCHRONIZATION SIGNALS
%
% d =  genSSS(nID1,nID2,nSF,SHOW)
%
% nID1 [0,...,167]  Physical-layer cell-identity group 
% nID2 [0,1,2]      Physical-layer identity within the cell-ID group   
% nSF  [0,1]        Subframe 0 = 0, subframe 10 = 1;
% SHOW (optional)   To show [nID1,m0,m1]
%
%==========================================================================
function d = genSSS(nID1,nID2,nSF,varargin)

% Generation of indices m0 and m1 derived from nID

q0 = floor(nID1/30);
q1 = floor((nID1+q0*(q0+1)/2)/30);
m  = nID1 + q1*(q1+1)/2;

m0 = mod(m,31);
m1 = mod(m0+floor(m/31)+1,31);

% Generation of the sequences s0 and s1

x  = [0 0 0 0 1 zeros(1,26)];    % Initial conditions
for i  = 1:26
    x(i+5) = mod(x(i+2)+x(i),2);
end
s  = 1 - 2*x;
s_ref = [1 1 1 1 -1 1 1 -1 1 -1 -1 1 1 -1 -1 -1 -1 ...
        -1 1 1 1 -1 -1 1 -1 -1 -1 1 -1 1 -1];

n  = 0:30;
s0 = s(mod(n + m0,31)+1);
s1 = s(mod(n + m1,31)+1);

% Generation of the scrambling sequences c0 and c1

x  = [0 0 0 0 1 zeros(1,26)];    % Initial conditions
for i  = 1:26
    x(i+5) = mod(x(i+3)+x(i),2);
end
c = 1 - 2*x;
c_ref = [1 1 1 1 -1 1 -1 1 -1 -1 -1 1 -1 -1 1 1 1 ...
        -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 1 -1];

c0 = c(mod(n + nID2,31)+1);
c1 = c(mod(n + nID2 + 3,31)+1);

% Generation of the scrambling sequences z0 and z1

x  = [0 0 0 0 1 zeros(1,26)];    % Initial conditions
for i  = 1:26
    x(i+5) = mod(x(i+4) + x(i+2) + x(i+1) + x(i),2);
end
z = 1 - 2*x;
z_ref = [1 1 1 1 -1 -1 -1 1 1 -1 -1 1 -1 -1 -1 -1 -1 ...
         1 -1 1 1 1 -1 1 1 -1 1 -1 1 -1 -1];

z0 = z(mod(n + mod(m0,8),31)+1);
z1 = z(mod(n + mod(m1,8),31)+1);

% Combination of two length-31 sequences to define the SSS
d = zeros(1,62);
if not(nSF)    % SUBFRAME 0
    d(2*n+1)   = s0 .* c0;
    d(2*(n+1)) = s1 .* c1 .* z0;
else           % SUBFRAME 5
    d(2*n+1)   = s1 .* c0;
    d(2*(n+1)) = s0 .* c1 .* z1;
end

if length(varargin)
    [nID1,m0,m1]
    % To check sequences -> Book: LTE for 4G Mobile Broadband
    [sum(s==s_ref) sum(c==c_ref) sum(z==z_ref)] 
end