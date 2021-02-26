%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates the encoded message using the E4 encoder and
% computes its trellis.
%
% Parameters : message : the message to be encoded.
% 
% Return : codedMessage : the encoded message.
%          trellis : the trellis of E4
%
% Example : message = [0 1 0 1 1 1 0 0 1]
%           [codedMessage, trellis] = convolutional(message)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [codedMessage, trellis] = E4(message)
    % Initialization
    D1 = 0;
    D2 = 0;
    D3 = 0;
    c1 = 0;
    codedMessage = [];
    trellis.numInputSymbols = 2^2;
    trellis.numOutputSymbols = 2^3;
    trellis.numStates  = 2^3;
    trellis.nextStates = zeros(trellis.numStates,trellis.numInputSymbols);
    trellis.outputs = zeros(trellis.numStates,trellis.numInputSymbols);
    
    for ii=1:2:length(message)
        u1 = message(ii);
        u2 = message(ii+1);
        symbol = bi2de([u1 u2],'left-msb');
        currentState = bi2de([D1 D2 D3],'left-msb');
        c1 = D3;
        D3 = mod(D2+u1,2);
        D2 = mod(D1+u2,2);
        D1 = c1;
        codedMessage = [codedMessage c1 u1 u2];
        nextState = bi2de([D1 D2 D3],'left-msb');
        trellis.nextStates(currentState+1,symbol+1) = nextState;
        trellis.outputs(currentState+1,symbol+1) = bi2de([c1 u1 u2],'left-msb');
    end
end