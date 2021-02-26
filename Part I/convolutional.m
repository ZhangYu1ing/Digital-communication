%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates the encoded message given the trellis of a
% convolutional code.
%
% Parameters : message : the message to be encoded.
%              trellis : the trellis of the convolutional code
% 
% Return : codedMessage : the encoded message.
%
% Example : message = [0 1 0 1 1 1 0 0 1]
%            trellis = generateTrellis(1,2,[5 7])
%            -> Rate 1/2 convolutional code with G = (1+D^2,1+D+D^2)
%            codedMessage = convolutional(message, trellis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function codedMessage = convolutional(message, trellis)
    n = log2(trellis.numOutputSymbols);     % Length of the codewords
    currentState = 0;      % Current state initialized to all-0
    codedMessage = [];
    for ii=1:length(message)
        inputBit = message(ii);
        codedMessage = [codedMessage oct2dec(optiDe2Bi(trellis.outputs(currentState+1,inputBit+1),n))]; % Coded sequence corresponding to inputBit
        currentState = trellis.nextStates(currentState+1,inputBit+1);     % Update of the current state
    end
end