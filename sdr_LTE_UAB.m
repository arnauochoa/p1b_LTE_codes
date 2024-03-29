%[UAB TA18-19] SPCOMNAV - J.A. del Peral-Rosado (18/03/16, rev. 05/03/2019)
%==========================================================================
% - Read IQ samples from USRP file
% - Detect cyclic and acquisition
%==========================================================================
clearvars; close all; clc;

% Settings
N_skip        = 300000;     %       # of samples to skip
N_RadioFrames = 100;        %       # of radio frames to read

% Core LTE parameters
nRB     = 6;                %       # of resource blocks
Fsc     = 15e3;             % [Hz]  Sub-carrier spacing (15 kHz)
NcRB    = 12;               %       # of sub-carriers per RB  
Ts      = 1/Fsc;            % [s]   OFDM time duration (66.667 us)
nCP     = 1;                % [0,1]	Normal CP = 1, extemded CP = 0
N       = 2^ceil(log2(NcRB*nRB));%  # of subcarriers / samples per symbol                    
Fs_0    = 2e6;              % [Hz]  Sampling frequency of the dongle
Fs      = N/Ts;             % [Hz]  Sampling frequency (MHz): 1.92, 3.84, 7.68, 15.36, 23.04, 30.72
T_cp    = 4.7e-6;           % [s]   Period of the cyclic prefix (CP)
T_cpe   = 5.2e-6;           % [s]   Period of the extended CP
L_cp    = round(Fs*T_cp);   %       # of samples of the CP
L_cpe   = round(Fs*T_cpe);  %       # of samples of the extended CP

% Define samples per radio frame
N_symb_frame   = 140;       %       Symbols per radio frame
N_symb_slot    = 6 + nCP;   %       Symbols per slot
N_samp_slot    = 7*N + 6*L_cp + L_cpe; %    Samples per slot
N_slot_frame   = 20;        %       Slots per radio frame
N_samp_frame   = N_slot_frame*N_samp_slot;% Samples per radio frame
N_cellSectors  = 3;         %       # of sectors per cell 
N_cellGroups   = 168;       %       # of groups per cell

% Define indexes of the synchronisation subcarriers
ind_SS = [(N-62)/2:N/2-1 N/2+1:(N+62)/2]+1;  

% Define filename to open
FileName = 'USRP_OneBS_Fs_2MHz_Gain_31dB_Fc_806MHz_301014.dat';

% Open file
fID = fopen(FileName);
if fID < 0, error('Invalid file name.'); end
fprintf('------- File name: %s\n',FileName); 

% Skip 'N_skip' samples to avoid recording problems at the initialization
fseek(fID,2*round(N_skip*Fs_0/Fs/2)*2,'bof');

% Read file
Tread   = N_RadioFrames*10e-3;
Nss_0   = Tread*Fs_0;   
Nss     = Tread*Fs;
S       = fread(fID,2*Nss_0,'short');

% Deinterleave I/Q samples
x = S(1:2:end)+1i*S(2:2:end);

% Resample Fs_0 to Fs
InterpolationFactor = Fs/gcd(Fs_0,Fs);
DecimationFactor = Fs_0/gcd(Fs_0,Fs); 
r = resample(x,InterpolationFactor,DecimationFactor);

%% UNCOMMENT ONLY for exercise 2.g) on the validation of the CFO estimation
% addCFO = 0.02; % Additional freq. offset
% r = r.*exp(j*(2*pi* addCFO *(0:length(r)-1)/N)).'; 

%% Exercise 1
%   1.a) Plot real part of the first N_samp_slot samples of the received 
%        signal r
    figure;
    plot(1:N_samp_slot, real(r(1:N_samp_slot)));
    xlabel("Time [sample]");
    ylabel("Received signal [?V]");

%   1.b) Use fft and fftshift to transform received signal in the frequency 
%        domain  
    R_symb  =   fftshift(fft(r(1:N)));
    R_slot  =   fftshift(fft(r(1:N_samp_slot)));

%   1.c) Compute the power spectrum of the frequency received signal over 
%        one symbol and over one slot in dB
    S_symb  =   10 * log10((abs(R_symb).^2)/N);
    S_slot  =   10 * log10((abs(R_slot).^2)/N_samp_slot);

%   1.d) Compare both power spectrum in the same plot
    B       =   Fs/2;
    f_symb  =   linspace(-B, B, length(S_symb))/1e6;
    f_slot  =   linspace(-B, B, length(S_slot))/1e6;

    figure;
    plot(f_symb, S_symb, 'b'); hold on;
    plot(f_slot, S_slot, 'r');
    legend("S_{symb}", "S_{slot}");
    xlabel("Frequency [MHz]");
    ylabel("Power spectrum [dB?/Hz]");

%   1.e) Identify the LTE pilot signal in the time-frequency spectrum
    plotTFspectrum(r, 0, 1.4);

%% Exercise 2. Implementation of van de Beek algorithm
%   2.a) Define index vector of samples over a slot
    m       =   1:N_samp_slot;

%   2.b) Compute rho
    num     =   mean(r(m) .* conj(r(m+N)));
    denom   =   sqrt(mean(abs(r(m)).^2) * mean(abs(r(m+N)).^2));
    rho     =   abs(num/denom);

%   2.c) Compute gamma and (capital) phi
% rho     =   zeros(1, N_samp_slot);
    gamma   =   zeros(1, N_samp_slot);
    phi     =   zeros(1, N_samp_slot);
    for m = 1:N_samp_slot
        k           =   m:m+L_cp-1;
        gamma(m)    =   sum(r(k) .* conj(r(k+N)));
        phi(m)      =   1/2 * sum(abs(r(k)).^2 + abs(r(k+N)).^2);
    end

    figure;
    plot(abs(gamma));
    title('$\gamma(m)$', 'Interpreter', 'latex');
    figure;
    plot(phi);
    title('$\Phi(m)$', 'Interpreter', 'latex');

%   2.d) Compute and plot the log-likelihood function
    lambda = abs(gamma) - (rho .* phi);
    figure;
    plot(lambda);
    title('Log-likelihood function $\Lambda(\theta, \hat{\varepsilon}_{ML}(\theta))$ ', 'Interpreter', 'latex');

%   2.e) Compute estimation of the coarse timing offset (phi)
%   Name this variable: 'startCP'
    [~, startCP]    =   max(lambda);

%   2.f) Compute estimation of residual frequency offset
%   Name this variable: 'freq_res'
    freq_res        =   -1/(2*pi) * angle(gamma(startCP));

% for i = 1:6
%     saveas(figure(i), sprintf('figure %d', i));
% end
    close all;

    fprintf('Coarse timing offset = %d\n',startCP);
    fprintf('Normalized residual frequency = %.4f\n',freq_res);

%% Exercise 3. Compensate timing and residual carrier frequency 
%   Name the signal with compensted errors as 's';
    m   =   startCP:(length(r) - 1);
    s   =   r(m) .* exp(-j * (2 * pi * freq_res * m/N)).';

%% Exercise 4. Cell sector detector by correlating PSS
%   4.a) Define a vector with the first sample of the OFDM symbols over a
%   radio frame.
% Indices of the first samples of every symbol in a slot
% i.e. 11 : 128+9 : 11+6*(128+9)
    first_ind_slot         =   L_cpe+1 : N+L_cp : (L_cpe+1)+(N_symb_slot-1)*(N+L_cp);

    first_ind   =   nan(1, N_symb_slot*N_slot_frame/2);
    for nSlot = 0:(N_slot_frame/2) - 1
        % Indices indicating symbols on a slot, i.e. 1 2 3 4 5 6 7, 8 9 10 ...
        slot_samp   =   1 + (N_symb_slot*nSlot):N_symb_slot + (N_symb_slot*nSlot);
        % Indices of the first samples of every symbol in half a frame
        % i.e. 11 148 285 422 559 696 833, 971 1108 ...
        first_ind(slot_samp) = first_ind_slot + (N_samp_slot*nSlot);
    end

%   4.b) Correlate in the frequency domain the received samples with the
%   PSS sequence generated with the genPSS(nCell) function
    % Initialize vector with N_cellSectors and N_symb_slot*N_slot_frame/2 zeros
    maxPeak     =   zeros(N_cellSectors, N_symb_slot*N_slot_frame/2);

    for nCell = 1:N_cellSectors
        % Initialize PSS sequence with N zeros
        PSS_seq         =   zeros(N, 1);
        % Use ind_SS vector to index PSS sequence vector with the PSS chips
        %   obtained by using �genPSS.m� function        
        PSS_seq(ind_SS) =   genPSS(nCell-1);
        
        for nInit = 1:N_symb_slot*N_slot_frame/2
            % Apply fft and fftshift to OFDM received symbols of N samples 
            sym_t   =   s(first_ind*nInit:first_ind*nInit+N-1);
            sym_f   =   fftshift(fft(sym_t));
            % Multiply subcarrier by subcarrier the frequency received 
            %   signal with the conjugate of the PSS sequence, by using �.*�    
            corr_f  =   sym_f .* conj(PSS_seq);
            % Apply ifft and fftshift to correlated signal in order to
            %   obtain the cross-correlation function
            corr_t  =   fftshift(ifft(corr_f));
            % Find maximum and absolute value of the CCF and save it in 
            %   �maxPeak(nCell,nInit)�
            maxPeak(nCell,nInit) = max(abs(corr_t));
        end
    end
    
    figure;
    plot(maxPeak.');
    legend("ID = 0", "ID = 1", "ID = 2");
    title("Cross-correlation functions of signal and PSS sequences");
    xlabel("Samples"); ylabel("Correlation");
    
%   4.c) Find cell ID sector with maximum of 'maxPeak
%   Name the cell ID sector variable as nID2
%   4.d) Find first sample of the OFDM symbol corresponding to the PSS
%   Name this index 'startPSS'
    [M, I]          =   max(maxPeak, [], 1);
    
    %- Sample where PSS starts
    [~, startPSS]   =   max(M);
    %- Cell ID
    nID2            =   I(startPSS) - 1;

    fprintf('Cell Sector ID = ID2 = %d\n',nID2);
    fprintf('Sample where PSS starts = %d\n',startPSS);

%   4.e) Estimate integer CFO
    sc_offset       =   -5:5;
    PSS_seq         =   zeros(N, 1);    
    PSS_seq(ind_SS) =   genPSS(nCell-1);
    pss_sym_t       =   s(startPSS:startPSS + N - 1);
    m               =   0:length(pss_sym_t)-1;
    maxPeak         =   zeros(length(sc_offset), 1);
    for nSubCarr =  1:length(sc_offset)
        shift_sym_t     =   pss_sym_t .* exp(-j * (2 * pi * sc_offset(nSubCarr) * m/N)).';
        shift_sym_f     =   fftshift(fft(shift_sym_t, N));
        corr_f          =   shift_sym_f .* conj(PSS_seq);
        corr_t          =   fftshift(ifft(corr_f));
        maxPeak(nSubCarr) = max(abs(corr_t));
        
        figure;
        plot(abs(shift_sym_f)/max(abs(shift_sym_f)));
        hold on
        plot(abs(conj(PSS_seq)))
        title(sprintf("subcarr = %d", sc_offset(nSubCarr)));
    end

%     for i = 2:length(sc_offset)+1
%         saveas(figure(i), sprintf("figure %d.png", i));
%         close(figure(i));
%     end
    
    figure;
    plot(sc_offset, maxPeak.');
    title("Cross-correlation functions of OFDM symbol and PSS sequence");
    xlabel("Frequency offset (subcarrier)"); ylabel("Correlation");
    
    [~, offset_samp]    =   max(maxPeak);
    integerCFO          =   sc_offset(offset_samp) * Fsc;
    fprintf('Subcarrier offset = %d subcarriers\n', sc_offset(offset_samp));
    
%   4.f) Obtain CFO
%   Name CFO variable 'freq_offset'
    freq_offset         =   integerCFO + freq_res;

    fprintf('Normalized frequency offset = %.4f\n',freq_offset);
    
%   4.g) Compensate frequency offset
    pss_sym_t           =   pss_sym_t .* exp(-j * (2 * pi * sc_offset(offset_samp) * m/N)).';
  
    

