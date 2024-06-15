filename = './phantom/phantom_info_r3mm_is1450_bg1540_depth15.mat';
load(filename);
tone_burst_cycles = 4;
source.p_mask = mk_linear_array(arrayx, N_elements, pitch, Nx, Ny, dy);
sensor.mask = mk_linear_array(arrayx, N_elements, pitch, Nx, Ny, dy);

axiscor = ((1 : N_elements) - mean(1 :N_elements)) * pitch;

%% RUN THE SIMULATION


% set the input settings
input_args = {...
    'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false, ...
    'DisplayMask', source.p_mask, 'PlotScale', 'auto'};

% run the simulation if set to true, otherwise, load previous results

for elecount = 1 : N_elements

    % update the command line status
    disp('');
    disp(['Computing scan line ' num2str(elecount) ' of ' num2str(N_elements)]);

    % update the current steering angle
    ypos = round((elecount-(N_elements+1)/2)*pitch/dy+Ny/2);
    % create the tone burst signals
    input_signal = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles);
    [np , siglength] = size(input_signal);
    source.p_mask = zeros(Nx , Ny);
    source.p_mask(arrayx , ypos) = 1;
    source.p = input_signal;

    % run the simulation
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

    rf_frame.t_array = kgrid.t_array;
    rf_frame.rf_data = sensor_data;
    
    filename = ['./rf_signal/r3mm_is1450_bg1540_depth15_line_' num2str(elecount) '.mat'];
    save(filename, 'rf_frame')

end
    % save the scan lines to disk
    
