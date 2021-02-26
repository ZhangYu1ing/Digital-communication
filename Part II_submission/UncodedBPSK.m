close all
clear all

%%%%%%%%%%%%%%%% Uncoded BPSK %%%%%%%%%%%%%%%%
EbN0 = [-1:8];
for ii = 1:length(EbN0)
    %% Generation of the message
    N = 1000000;
    bits = randi([0 1], 1, N);

    %% BPSK Mapping
    symbols = 1-2*bits;

    %% AWGN Noise and Channel
    Es = mean(abs(symbols).^2);     % Energy per symbol
    sigma2 = Es/(2*(10^(EbN0(ii)/10)));     % DSP of the noise per branch
    noise = sqrt(sigma2)*(randn(1,size(symbols,2))+1j*randn(1,size(symbols,2)));
    signal = symbols + noise;               % Noisy symbols
    
    %% Symbol Decoder
    bits_dec = real(signal) < 0;      % Detector : the bit is 1 if the condition is true
    
    %% BER
    errors = numel(find(bits_dec ~= bits));     % Count the bit errors
    BER(ii) = errors/N;
    BER_theory(ii) = qfunc(sqrt(2*10^(EbN0(ii)/10)));

end

%% Results
figure('name', 'BER of an uncoded BPSK')
semilogy(EbN0, BER, 'r-+')
hold on
semilogy(EbN0, BER_theory, 'b-*')
ylim([1e-4 1])
xlabel('Eb/N0 [dB]');
ylabel('Bit Error Rate');
legend('Simulated Uncoded BPSK', 'Theoretical Uncoded BPSK');
grid on;