%% Parameters
N = 64;             % Number of subcarriers
CP = 16;            % Cyclic prefix length
D = 4;              % Number of OFDM symbols
EbN0 = 20;          % Signal-to-noise ratio in dB

%% Generate random data
data = randi([0 1], 1, N*D);

%% OFDM modulation
ofdm_signal = ofdm_mod(data, N, CP, D);

%% Add noise to the signal
SNR = 10^(EbN0/10);     % Convert dB to linear scale
sigma = sqrt(1/(2*SNR)); % Compute noise standard deviation
noise = sigma*(randn(1,length(ofdm_signal))+1i*randn(1,length(ofdm_signal)));
ofdm_signal_noisy = ofdm_signal + noise;

%% OFDM demodulation
received_data = ofdm_demod(ofdm_signal_noisy, N, CP, D);

%% Calculate bit error rate (BER)
ber = sum(xor(data,received_data))/(N*D);

%% Display results
fprintf('Eb/N0 = %d dB\n', EbN0);
fprintf('Bit error rate = %e\n', ber);

%% Plot the signals
figure;
subplot(2,1,1);
plot(real(ofdm_signal));
hold on;
plot(real(ofdm_signal_noisy));
hold off;
title('Transmitted and Received Signals (Real Part)');
legend('Transmitted Signal', 'Received Signal (with Noise)');
xlabel('Time');
ylabel('Amplitude');

subplot(2,1,2);
plot(imag(ofdm_signal));
hold on;
plot(imag(ofdm_signal_noisy));
hold off;
title('Transmitted and Received Signals (Imaginary Part)');
legend('Transmitted Signal', 'Received Signal (with Noise)');
xlabel('Time');
ylabel('Amplitude');

%% OFDM modulation function
function ofdm_signal = ofdm_mod(data, N, CP, D)
% OFDM modulation with cyclic prefix

% Reshape the data into D blocks of length N
data_matrix = reshape(data, N, D).';

% Perform the FFT on each block
subcarriers = fft(data_matrix, [], 2);

% Add the cyclic prefix to each block
subcarriers_CP = [subcarriers(:,end-CP+1:end) subcarriers];

% Convert the subcarriers to a time-domain signal
ofdm_signal = ifft(subcarriers_CP, [], 2);

% Reshape the time-domain signal into a vector
ofdm_signal = ofdm_signal(:).';
end

%% OFDM demodulation function
function received_data = ofdm_demod(ofdm_signal, N, CP, D)
% OFDM demodulation with cyclic prefix removal

% Reshape the OFDM signal into D blocks of length N+CP
ofdm_signal_matrix = reshape(ofdm_signal, N+CP, D).';

% Remove the cyclic prefix from each block
subcarriers_received = ofdm_signal_matrix(:,CP+1:end);

% Perform the IFFT on each block
received_data_matrix = ifft(subcarriers_received, [], 2);

% Convert the received data to a vector
received_data = received_data_matrix(:).';

% Threshold the received data to get the bits
received_data = real(received_data) > 0.5;
end
