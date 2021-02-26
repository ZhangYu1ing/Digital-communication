% ======================================================================= %
% SSY125 Project
% ======================================================================= %
clc
clear all
close all

% ======================================================================= %
% Simulation Options
% ======================================================================= %
N = 1e5;  % simulate N bits each transmission (one block)
maxNumErrs = 100; % get at least 100 bit errors (more is better)
maxNum = 1e6; % OR stop if maxNum bits have been simulated
EbN0 = -1:8; % power efficiency range


% ======================================================================= %
% Other Options
% ======================================================================= %
% ...

% ======================================================================= %
% Simulation Chain
% ======================================================================= %
BER = zeros(1, length(EbN0)); % pre-allocate a vector for BER results
BER_theory = zeros(1,length(EbN0));
const = [1+1i 1-1i -1+1i -1-1i]/sqrt(2);

% figure(1)
% plot(const,'ro')  % QPSK
% hold on
for i = 1:length(EbN0) % use parfor ('help parfor') to parallelize  
  totErr = 0;  % Number of errors observed
  totErr_un = 0;
  num = 0; % Number of bits processed
  num_un = 0 ;
  while((totErr < maxNumErrs) && (num < maxNum))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  %% [SRC] generate N information bits 
  bits = randi([0 1],1,N); 

  %% [ENC] convolutional encoder
  rate = 1;
  switch rate
      case 1
          input = 1;
          output = 2;
      case 2
          input = 2;
          output = 3;
  end
  
  number = 1;  % switch different type of encoder
  switch number
      case 1
          poly = [5 7];    %G = [1 0 1; 1 1 1];
      case 2
          poly = [23 22];  %G = [1 0 1 1 1; 1 0 1 1 0];
      case 3
          poly = [23 33];  %G = [1 0 0 1 1; 1 1 0 1 1];   
      
  end
  trellis = generateTrellis(input,output,poly);
  codeMessage = convolutional(bits,trellis); % codeMessage: encoded message
  
  %% [MOD] symbol mapper [Grey Mapping in QPSK]
  ak = 1-2*codeMessage(1:2:2*N);  % Real part
  bk = 1-2*codeMessage(2:2:2*N);  % Imaginary part
  symbols = ak+1j*bk;
  

  %% [CHA] add Gaussian noise
  %coded 
  %Es = mean(abs(symbols).^2);
  Es = 1;
  %sigma2 = 1*10^(-EbN0(i)/10);
  sigma2 = Es/((10^(EbN0(i)/10)));     % DSP of the noise per branch
  noise = sqrt(sigma2)*(randn(1,size(symbols,2))+1j*randn(1,size(symbols,2)));
  signal = symbols + noise;               % Noisy symbols
  
  
  % [HR] Hard Receiver
  ak_dec = real(signal) < 0;      % The decoded bit is 1 if the condition is true
  bk_dec = imag(signal) < 0;
  bits_dec = zeros(1,2*N);     % Decoded bits
  bits_dec(1:2:2*N) = ak_dec;
  bits_dec(2:2:2*N) = bk_dec;
  
  hard_decoded = viterbiDecoding(bits_dec,trellis);
  BitErrs = numel(find(hard_decoded ~= bits));
  totErr = totErr + BitErrs;
  
  num = num + N;
  %% Uncoded 
  ak_un = 1-2*bits(1:2:N);  % Real part
  bk_un = 1-2*bits(2:2:N);  % Imaginary part
  symbols_un = ak_un+1j*bk_un;
  
    % [CHA] add Gaussian noise
  % uncoded
  % Es = mean(abs(symbols).^2);     % Energy per symbol
  sigma2_un = Es/(2*(10^(EbN0(i)/10)));     % DSP of the noise per branch
  noise_un = sqrt(sigma2_un)*(randn(1,size(symbols_un,2))+1j*randn(1,size(symbols_un,2)));
  signal_un = symbols_un + noise_un;
  
  ak_dec_un = real(signal_un) < 0;      % The decoded bit is 1 if the condition is true
  bk_dec_un = imag(signal_un) < 0;
  bits_dec_un = zeros(1,N);     % Decoded bits
  bits_dec_un(1:2:N) = ak_dec_un;
  bits_dec_un(2:2:N) = bk_dec_un;
  
  BitErrs_un = numel(find(bits_dec_un ~= bits)); % count the bit errors and evaluate the bit error rate
  totErr_un = totErr_un + BitErrs_un;
  
  %num_un = num_un + N;
  
  %% [SR] Soft Receiver
  % ...
  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== % 
  
  disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
  num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
  num2str(totErr/num, '%10.1e') '. +++']);
  end 
  disp(['Eb/N0 = ' num2str(i-2)]);
  BER_un(i) = totErr_un/num; 
  BER_HDD(i) = totErr/num;
  BER_theory(i) = qfunc(sqrt(2*10^(EbN0(i)/10)));
  
end
% title('Constellation of QPSK')
% legend('QPSK theory','symnols with noise','symbols without noise')

figure('name', 'BER of an uncoded and coded QPSK')
semilogy(EbN0, BER_un,'b*')
hold on
semilogy(EbN0, BER_HDD)
hold on
semilogy(EbN0, BER_theory)
ylim([1e-4 1])
xlabel('Eb/N0 [dB]');
ylabel('Binary Error Rate');
legend('Uncoded QPSK','Hard receiver coded QPSK', 'Theoretical Uncoded QPSK');
grid on;
% ======================================================================= %
% End
% ======================================================================= %