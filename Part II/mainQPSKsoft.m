% ======================================================================= %
% SSY125 Project
% ======================================================================= %
clc
clear

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
M = 4; % Modulation order
BER_soft = zeros(1, length(EbN0)); % pre-allocate a vector for BER results (soft)
UpperBound = zeros(1, length(EbN0));

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
  rate = 1;
  switch rate
      case 1
          input = 1;
          output = 2;
      case 2
          input = 2;
          output = 3;
  end
  
  number = 3;  % switch different type of encoder
  switch number
      case 1
          poly = [5 7];    %G = [1 0 1; 1 1 1];
      case 2
          poly = [26 27];  %G = [1 0 1 1 1; 1 0 1 1 0];
      case 3
          poly = [23 33];  %G = [1 0 0 1 1; 1 1 0 1 1];   
      
  end
  Rc = input/output; % Code rate
  capacity = 10*log10((2^(log2(M)*Rc)-1)/(log2(M)*Rc));
  trellis = generateTrellis(input,output,poly);
  codedMessage = convolutional(bits,trellis);

  % [MOD] symbol mapper -> QPSK (Gray mapping)
  ak = 1-2*codedMessage(1:2:end);  % Real part
  bk = 1-2*codedMessage(2:2:end);  % Imaginary part
  symbols = ak+1j*bk;

  % [CHA] add Gaussian noise
  Es = mean(abs(symbols).^2);     % Energy per symbol
  sigma2 = Es/(4*Rc*(10^(EbN0(i)/10)));     % DSP of the noise per branch
  noise = sqrt(sigma2)*(randn(1,size(symbols,2))+1j*randn(1,size(symbols,2)));
  signal = symbols + noise;      % Noisy symbols

  % [SR] Soft Receiver
  codedSequence_s = zeros(1,length(codedMessage));    % Sequence of the previous detected bits (coded)
  codedSequence_s(1:2:end) = real(signal);
  codedSequence_s(2:2:end) = imag(signal);
  softDec = viterbiDecodingSoft(codedSequence_s,trellis);
  
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
  
  % [UB] Upper Bound
%   Nb = 25;
%   spect = distspec(trellis,Nb);
%   for ll=1:Nb
%     UpperBound(i) = UpperBound(i) + spect.weight(ll)*qfunc(sqrt(2*(spect.dfree-1+ll)*Rc*10^(EbN0(i)/10)));
%   end
end

% ======================================================================= %
% Results
% ======================================================================= %
figure('name', 'BER of the coded system')
semilogy(EbN0, BER_soft, 'r-+')
% hold on
% semilogy(EbN0, UpperBound, 'r--')
hold on
xline(capacity, 'r-.');
ylim([1e-4 1])
xlabel('Eb/N0 [dB]');
ylabel('Bit Error Rate');
legend('Simulated Coded QPSK (soft)', 'Capacity');
grid on;

% ======================================================================= %
% End
% ======================================================================= %
