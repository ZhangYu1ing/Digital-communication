%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function converts a decimal number to an octal
%
% Parameter : decimal : the decimal number to convert
% 
% Return : octal : the result in octal
%
% Example : octal = dec2oct(23)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function octal = dec2oct(decimal)
    oct = [];
    while (decimal ~= 0)
        oct = [oct mod(decimal,8)];
        decimal = floor(decimal/8);
    end
    oct = fliplr(oct);
    octal = 0;
    for ii=1:length(oct)
        octal = octal + oct(ii)*10^(length(oct)-ii);
    end
end