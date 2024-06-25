function[t] = calc_relative_delay(sig1 , sig2 , fs)
% Calculate the relative delay between two signal using circulative shift and 
% ZNCC. The reltative delay is positive if sig2 lags behind sig1
    if ~(length(sig1) == length(sig2))
        error('signal length not the same!')
    end
    siglength = length(sig1);
    shiftvec = round(-siglength / 4 : siglength / 4);
    nc = 1;
    zncc_value = zeros(1 , siglength);
    for n_shift = shiftvec
        sig2_shift = circshift(sig2 , n_shift);
        sig1_centered = sig1 - mean(sig1);
        sig2_centered = sig2_shift - mean(sig2_shift);
        numerator = sum(sig1_centered .* sig2_centered);
        denominator = sqrt(sum(sig1_centered .^ 2) * sum(sig2_centered .^ 2));
        zncc_value(nc) = numerator / denominator;
        nc = nc + 1;
    end
    [~ , I] = max(zncc_value);
    t = shiftvec(I) / fs;
end