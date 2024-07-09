function [rf_fsa] = decode_encoded_beams(rf, transmit_delays, alpha, varargin)
    % Initialize variables for optional parameters
    center_freq = [];
    bandwidth = 0.3; % Set a default value for bandwidth
    
    % Check for additional parameters in varargin
    if nargin > 3
        for k = 1:2:length(varargin)
            if strcmp(varargin{k}, 'center_freq')
                center_freq = varargin{k+1};
            elseif strcmp(varargin{k}, 'bandwidth')
                bandwidth = varargin{k+1};
            end
        end
    end
    
    [~, n_receives, n_samples] = size(rf);
    [~, n_elements] = size(transmit_delays);
    
    % Frequency domain implementation
    RF = fft(rf, [], 3);
    frequency = (0:n_samples-1)/n_samples;
    
    % Apply bandpass filter if 'center_freq' is provided

    freq_ind = 1 : ceil(n_samples / 2); % Compute only half, assume symmetry
    freq_weight = ones(1 , ceil(n_samples / 2)); % Frequency domain weight
    if ~isempty(center_freq)
        freq_vec = linspace(0 , pi , ceil(n_samples / 2));
        freq_ind = freq_ind(abs((freq_vec - center_freq) ./ center_freq) < bandwidth);
        
        % 定义滚降系数 alpha，可以根据需要调整
        alpha = 0.2;
        
        % 计算升余弦窗
        n = linspace(-1 , 1 , length(freq_ind)); % 创建样本索引
        win = 1 .* (abs(n) <= (1-alpha) / 2) + ...
            0.5 .* (1 + cos((pi / alpha) .* (abs(n) - (1 - alpha)/2))) ...
            .* ((abs(n) > (1 - alpha)/2) & (abs(n) <= (1 + alpha)/2));
        freq_weight = win;
    end
    
    % Apply decoding matrix at each frequency
    RF_adj = zeros(n_receives, n_elements, n_samples, 'single');
    
    wc = 1;
    for i = freq_ind
        omega = 2*pi*frequency(i);
        H = exp(-1j*omega*transmit_delays);
        % Apply the decoding matrix
        Hinv=1\(H'*H+alpha*norm(squeeze(RF(:,:,i)))*eye(n_elements))*H';
        RF_adj(:,:,i) = Hinv*RF(:,:,i) .* freq_weight(wc);
        wc = wc + 1;
    end
    
    % Inverse FFT for real signal
    rf_fsa = ifft(RF_adj, [], 3, 'symmetric');
end