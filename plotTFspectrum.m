%==========================================================================
% SPCOMNAV research group - Universitat Autònoma de Barcelona (UAB)
%   
%   Plot Time-Frequency LTE Spectrum
%   
%   function plotTFspectrum(s,BOF,BWsys)
%
%   Inputs:
%       s                  Vector of LTE signal samples
%       BOF                Sample index of the beginning of frame (BOF)
%       BWsys       [MHz]  System bandwidth configuration
%                          (e.g. 1.4, 3, 5, 10, 15 and 20 MHz)
%
%   Example to plot the time-frequency spectrum of an LTE signal vector
%   from the first sample and using a 1.4 MHz of system bandwidth:
%   	plotTFspectrum(s,0,1.4)
%
%   Author: J.A. del Peral-Rosado                 [05/03/2019, version 1.0]
%==========================================================================
function plotTFspectrum(s,BOF,BWsys)

% LTE system bandwidth
switch BWsys    % [MHz]  System bandwidth configuration
    case 1.4,   nRB = 6;   % # of resource blocks
    case 3,     nRB = 15;  % # of resource blocks
    case 5,     nRB = 25;  % # of resource blocks
    case 10,    nRB = 50;  % # of resource blocks
    case 15,    nRB = 75;  % # of resource blocks
    case 20,    nRB = 100; % # of resource blocks
    otherwise, error(['=> LTE system bandwidth is here defomed in MHz and'...
                      'needs to be one of the following values: '...
                      '1.4, 3, 5, 10, 15, 20']);
end

% Core LTE parameters
Fsc    = 15e3;      % [Hz]   Sub-carrier spacing (15 kHz)
NcRB   = 12;        %       # of sub-carriers per RB  
Ts     = 1/Fsc;     % [s]   OFDM time duration (66.667 us)
N      = 2^ceil(log2(NcRB * nRB));  % # of subcarriers / samples per symbol                    
Fs     = N/Ts;      % [Hz]  Sampling frequency (MHz): 1.92, 3.84, 7.68, 15.36, 23.04, 30.72
T_cp    = 4.7e-6;   % [s]   Period of the cyclic prefix (CP)
T_cpe   = 5.2e-6;   % [s]   Period of the extended CP
L_cp    = round(Fs*T_cp); % # of samples of the CP
L_cpe   = round(Fs*T_cpe);% # of samples of the extended CP

% Define samples per radio frame
N_symb_frame   = 140;       % Symbols per radio frame
N_samp_slot    = 7*N + 6*L_cp + L_cpe;  % Samples per slot
N_slot_frame   = 20;        % Slots per radio frame
N_samp_frame   = N_slot_frame*N_samp_slot;  % Samples per radio frame

% Check length of signal vector
if numel(s) <= N_samp_frame
    error(['=> There is not sufficient number of samples to plot '...
           'the time-frequency spectrum of one LTE radio frame, i.e. ' ...
           num2str(N_samp_frame) ' samples for a system bandwidth of ' ...
           num2str(BWsys) ' MHz.']);
end
 
% Generate FFT samples of the LTE radio frame
NSymb = 1:N_symb_frame;
nSlot = mod(floor((NSymb-1)/7),20);
nInitSamp = L_cpe+N_samp_slot*nSlot+(N+L_cp)*mod(NSymb-1,7);
indSamples = BOF + ones(N,1)*nInitSamp +(1:N).'*ones(1,N_symb_frame);
    
% Apply FFT to indexed signal samples
s_frame = s(indSamples);
S_frame = fftshift(fft(s_frame),1)/sqrt(N);
      
% Compute power spectral density (PSD)
S_frame_psd = 10*log10((abs(S_frame*sqrt(N)/Fs).^2)*Fsc);

% Normalize PSD
S_frame_psd_norm = S_frame_psd-max(max(S_frame_psd));

% Plot time-frequency LTE spectrum
textSize = 12; % default text size
figure;
hP = pcolor((0:N_symb_frame-1)-.5,(-N/2:N/2-1)+.5,S_frame_psd_norm);
if BWsys > 5, set(hP, 'EdgeColor', 'none'); end
xlabel('OFDM symbols [Time]','FontSize',textSize);
ylabel('Subcarriers [Frequency]','FontSize',textSize);
hCol = colorbar; ylabel(hCol,'Power spectral density (dB/Hz)','FontSize',textSize);
colormap(flipud(hot)),
caxis([-60 40])
title('LTE radio frame','FontSize',textSize)
set(gca,'FontSize',textSize);