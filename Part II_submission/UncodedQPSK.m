close all
clear all

%%%%%%%%%%%%%%%% Uncoded QPSK %%%%%%%%%%%%%%%%
EbN0 = [-1:8];
for ii = 1:length(EbN0)
    %% Generation of the message
    N = 1000000;
    bits = randi([0 1], 1, N);

    %% QPSK Mapping
    ak = 1-2*bits(1:2:N);  % Real part
    bk = 1-2*bits(2:2:N);  % Imaginary part
    symbols = ak+1j*bk;

    %% AWGN Noise and Channel
    Es = mean(abs(symbols).^2);     % Energy per symbol
    sigma2 = Es/(4*(10^(EbN0(ii)/10)));     % DSP of the noise per branch
    noise = sqrt(sigma2)*(randn(1,size(symbols,2))+1j*randn(1,size(symbols,2)));
    signal = symbols + noise;               % Noisy symbols
    
    %% Symbol Decoder (using Gray Mapping)
    ak_dec = real(signal) < 0;      % Detector : the bit is 1 if the condition is true
    bk_dec = imag(signal) < 0;
    bits_dec = zeros(1,N);     % Detected bits
    bits_dec(1:2:N) = ak_dec;
    bits_dec(2:2:N) = bk_dec;
    
    %% BER
    errors = numel(find(bits_dec ~= bits));     % Count the bit errors
    BER(ii) = errors/N;
    BER_theory(ii) = qfunc(sqrt(2*10^(EbN0(ii)/10)));

end

%% Results
figure('name', 'BER of an uncoded QPSK')
semilogy(EbN0, BER, 'r-+')
hold on
semilogy(EbN0, BER_theory, 'b-*')
ylim([1e-4 1])
xlabel('Eb/N0 [dB]');
ylabel('Bit Error Rate');
legend('Simulated Uncoded QPSK', 'Theoretical Uncoded QPSK');
grid on;