filename = '../phantom/phantom_info_r3mm_is1450_bg1540_depth15_x0.mat';
% filename = '../phantom/phantom_info_homo.mat';
load(filename);
tone_burst_cycles = 4;
source.p_mask = mk_linear_array(arrayx, N_elements*4-1, pitch/4, Nx, Ny, dy);

sensor_xpos = xvec * 0; sensor_xpos(50 : 5 : 550) = 1;
sensor_xvec = xvec(sensor_xpos == 1);
sensor_zpos = zvec * 0; sensor_zpos(50 : 5 : 650) = 1;
sensor_zvec = zvec(sensor_zpos == 1);
sensor.mask = sensor_zpos * sensor_xpos';

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
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', true, ...
    'DisplayMask', source.p_mask, 'PlotScale', 'auto' , 'LogScale' , true};

% run the simulation if set to true, otherwise, load previous results

phi = 5 * pi / 180;

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

filename = ['../field_signal/r3mm_is1450_bg1540_depth15_x0_pw_255rx_phi_' num2str(phi / pi * 180) '.mat'];
% filename = ['../field_signal/homo_pw_255rx_phi_' num2str(phi / pi * 180) '.mat'];
save(filename, 'rf_frame' , 'sensor_xvec' , 'sensor_zvec' , 'sensor_xpos' , 'sensor_zpos')
% save the scan lines to disk
    
