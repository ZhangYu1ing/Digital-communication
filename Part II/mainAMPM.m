% ======================================================================= %
% SSY125 Project
% ======================================================================= %
clc
clear

% ======================================================================= %
% Simulation Options
% ======================================================================= %
N = 3e5;  % simulate N bits each transmission (one block)
maxNumErrs = 100; % get at least 100 bit errors (more is better)
maxNum = 3e6; % OR stop if maxNum bits have been simulated
EbN0 = -1:8; % power efficiency range

% ======================================================================= %
% Other Options
% ======================================================================= %
% ...

% ======================================================================= %
% Simulation Chain
% ======================================================================= %
M = 8; % Modulation order
BER_soft = zeros(1, length(EbN0)); % pre-allocate a vector for BER results (soft)

for i = 1:length(EbN0) % use parfor ('help parfor') to parallelize  
  totErr = 0;  % Number of errors observed (soft)
  num = 0; % Number of bits processed

  while((totErr < maxNumErrs) && (num < maxNum))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  % [SRC] generate N information bits 
  bits = randi([0 1], 1, N);

  % [ENC] convolutional encoder
  input = 2;
  output = 3;
  Rc = input/output; % Code rate
  capacity = 10*log10((2^(log2(M)*Rc)-1)/(log2(M)*Rc));
  [codedMessage, trellis] = E4(bits);
  
  % [MOD] symbol mapper -> AMPM
  ak = 2*(1-2*codedMessage(3:3:length(codedMessage)));
  bk = 2*(2*(mod(codedMessage(2:3:length(codedMessage))+codedMessage(3:3:length(codedMessage)),2))-1);
  symbols = (-1+2*codedMessage(1:3:length(codedMessage))-1j*(-1+2*codedMessage(1:3:length(codedMessage)))) + ak+1j*bk;

  % [CHA] add Gaussian noise
  %Es = mean(abs(symbols).^2);     % Energy per symbol
  Es = 1;
  %sigma2 = Es/(6*Rc*(10^(EbN0(i)/10)));     % DSP of the noise per branch
  sigma2 = Es/(6*Rc*(10^(EbN0(i)/10))); 
  noise = sqrt(sigma2)*(randn(1,size(symbols,2))+1j*randn(1,size(symbols,2)));
  signal = symbols + noise;      % Noisy symbols

  % [SR] Soft Receiver
  codedSequence_s = zeros(1,2*length(signal));
  codedSequence_s(1:2:end) = real(signal);
  codedSequence_s(2:2:end) = imag(signal);
  softDec = viterbiDecodingSoft(codedSequence_s,trellis);     % Viterbi decoder
  
  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== %
  BitErrs = numel(find(softDec ~= bits)); % count the bit errors and evaluate the bit error rate
  totErr = totErr + BitErrs;
  num = num + N;

  disp(['Eb/N0 = ' num2str(EbN0(i))]);
  disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr/num, '%10.1e') '. +++']);

 end 
  BER_soft(i) = totErr/num;
end

% ======================================================================= %
% Results
% ======================================================================= %
figure('name', 'BER of coded system')
semilogy(EbN0, BER_soft, 'r-+')
hold on
xline(capacity, 'r-.');
ylim([1e-4 1])
xlabel('Eb/N0 [dB]');
ylabel('Bit Error Rate');
legend('Simulated Coded AMPM (soft)', 'Capacity');
grid on;

% ======================================================================= %
% End
% ======================================================================= %
