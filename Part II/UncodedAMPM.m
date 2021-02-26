close all
clear all

%%%%%%%%%%%%%%%% Uncoded AMPM %%%%%%%%%%%%%%%%
EbN0 = [-1:12];
for ii = 1:length(EbN0)
    %% Generation of the message
    N = 1000000*3;
    bits = randi([0 1], 1, N);

    %% AMPM Mapping
    ak = 2*(1-2*bits(3:3:N));
    bk = 2*(2*(mod(bits(2:3:N)+bits(3:3:N),2))-1);
    symbols = (-1+2*bits(1:3:N)-1j*(-1+2*bits(1:3:N))) + ak+1j*bk;
    
    %% AWGN Noise and Channel
    Es = mean(abs(symbols).^2);     % Energy per symbol
    %Es = 1;
    sigma2 = Es/(6*(10^(EbN0(ii)/10)));     % DSP of the noise per branch
    noise = sqrt(sigma2)*(randn(1,size(symbols,2))+1j*randn(1,size(symbols,2)));
    signal = symbols + noise;               % Noisy symbols

    %% Symbol Decoder by Minimum Eucledian distance detector
   constellation = [(1 - 1i) (-3 + 3i) (1 + 3i) (-3 - 1i) (3 - 3i) (-1 + 1i) (3 + 1i) (-1 - 3i)]; 
   metric = abs(repmat(signal.',1,1) - repmat(constellation, length(signal), 1)).^2; % compute the distance to each possible symbol
   [dis index] = min(metric, [], 2); % find the closest for each received symbol
   index = index - 1;
   bits_de = de2bi(index,3,'left-msb');
   bits_de = bits_de';
   bits_dec = bits_de(:)';
   
    %% BER
    errors = numel(find(bits_dec ~= bits));     % Count the bit errors
    BER(ii) = errors/N;
    BER_theory(ii) = (4/3)*qfunc(sqrt((6/5)*10^(EbN0(ii)/10)));
end

%% Results
figure('name', 'BER of an uncoded AMPM')
semilogy(EbN0, BER, 'r-+')
hold on
semilogy(EbN0, BER_theory, 'b-*')
ylim([1e-4 1])
xlim([-1 12])
xlabel('Eb/N0 [dB]');
ylabel('Bit Error Rate');
legend('Simulated Uncoded AMPM','Theretical Uncoded AMPM');
grid on;