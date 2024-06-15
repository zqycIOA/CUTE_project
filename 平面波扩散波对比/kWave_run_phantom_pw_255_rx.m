filename = './phantom/phantom_info_r3mm_is1450_bg1540_depth15_x-5.mat';
% filename = './phantom/phantom_info_homo_with_target.mat';
load(filename);
tone_burst_cycles = 4;
source.p_mask = mk_linear_array(arrayx, N_elements*4-1, pitch/4, Nx, Ny, dy);
sensor.mask = mk_linear_array(arrayx, N_elements*2-1, pitch/2, Nx, Ny, dy);

xsize = dx * Nx;
ysize = dy * Ny;
xvec = linspace(-1,1,Nx) * xsize / 2;
yvec = linspace(-1,1,Ny) * ysize / 2;
[ys, xs] = meshgrid(1:Ny,1:Nx);
[yposes, xposes] = meshgrid(yvec, xvec);
ypos = yposes(source.p_mask==1);

axiscor = ((1 : N_elements) - mean(1 :N_elements)) * pitch;

%% RUN THE SIMULATION


% set the input settings
input_args = {...
    'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false, ...
    'DisplayMask', source.p_mask, 'PlotScale', 'auto'};

% run the simulation if set to true, otherwise, load previous results
phivec = (-5:0.5:5) * pi / 180;
pc = 1;
for phi = phivec

    % update the command line status
    disp('');
    disp(['Computing scan line ' num2str(pc) ' of ' num2str(length(phivec))]);

    % update the current steering angle
    tone_burst_offset = 1500 + ypos * sin(phi) ./ (1540 * kgrid.dt);

    % create the tone burst signals
    input_sig = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles, ...
        'SignalOffset', tone_burst_offset);
    input_sig = input_sig .* hanning(511);
    source.p = input_sig;

    % run the simulation
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

    rf_frame.t_array = kgrid.t_array - 1500 * kgrid.dt;
    rf_frame.rf_data = sensor_data;
    rf_frame.phi = phi;
    
    filename = ['./rf_signal/r3mm_is1450_bg1540_depth15_x-5_pw_255rx_line_' num2str(pc) '.mat'];
%     filename = ['./rf_signal/homo_with_target_pw_line_' num2str(pc) '.mat'];
    save(filename, 'rf_frame')
    pc = pc + 1;

end
    % save the scan lines to disk
    
