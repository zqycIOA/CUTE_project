%% Produce plane wave RF echo signal in heterogeneous medium
clear;

% simulation settings
DATA_CAST       = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations
% DATA_CAST       = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION  = true;         % set to false to reload previous results instead of running simulation

%% DEFINE THE K-WAVE GRID & ARRAY PROPERTIES
% set array parameters
N_elements = 128;
f0 = 5e6;
c0 = 1540;                      % [m/s]
lambda = c0/f0;
pitch = 300e-6;

% set the size of the perfectly matched layer (PML)
pml_x_size = 74;                % [grid points]
pml_y_size = 68;                % [grid points]

% set total number of grid points not including the PML
Nx = 800;     % [grid points]
Ny = 600;     % [grid points]

% calculate the spacing between the grid points
dx = 75e-6;                    % [m]
dy = 75e-6;                    % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

xsize = dx * (Nx-1);
ysize = dy * (Ny-1);
xvec = kgrid.x_vec;
yvec = kgrid.y_vec;
% create the time array
kgrid.dt = 6.25e-9;
sampling_freq = 1 / kgrid.dt;
% kgrid.makeTime(c0, [], t_end);
kgrid.t_array=(0 : 12000-1) * kgrid.dt;

%% DEFINE THE MEDIUM PARAMETERS


%% DEFINE THE INPUT SIGNAL

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = f0;     % [Hz]
tone_burst_cycles = 4;

%% DEFINE THE MEDIUM PROPERTIES
% define a random distribution of scatterers for the medium


% define properties
arrayx = 10;
xvec = kgrid.y_vec; zvec = kgrid.x_vec; zvec = zvec - zvec(1) - arrayx * dx;
[xs,zs] = meshgrid(xvec , zvec);
inclusion_map = makeDisc(Nx, Ny, round(15e-3 / dx) + arrayx , round(-5e-3 / dy) + Ny/2 , round(3e-3 / dx));
% inclusion_map = (sqrt(0.2 * (xs + 5e-3).^2 + (zs - 15e-3).^2) <= 4e-3) |...
%     (sqrt(0.2 * (xs - 5e-3).^2 + (zs - 15e-3).^2) <= 4e-3); 
% inclusion_map = sqrt((xs / 2).^2 + (zs - 15e-3).^2) < 4e-3;
hannfilt = hanning(20) * hanning(20)';
hannfilt = hannfilt ./ sum(sum(hannfilt));
sound_speed_map = conv2(inclusion_map * (1450 - 1540) , hannfilt , 'same') + 1540;

% point_target_map = 0;
% for pointz = 10e-3 : 5e-3 : 30e-3
%     point_target_map = point_target_map + makeDisc(Nx , Ny , round(pointz / dx) + arrayx , Ny/2 , round(0.1e-3 / dx));
% end
% 
% sound_speed_map(point_target_map == 1) = 340;
% sound_speed_map = inclusion_map * (1530 - 1500) + (1 - inclusion_map) * 0;
% sound_speed_map = ones(Nx, Ny) .* c0;
% bg_map = (zs > 25e-3) .* 1510 + (zs < 20e-3) .* 1500 + (zs <= 25e-3) .* (zs >=20e-3) .* ...
%     ((zs - 20e-3) * ((1510 - 1500) / 5e-3) + 1500);
% sound_speed_map = conv2 (sound_speed_map , hannfilt , 'same') + bg_map;


density_map = 1000 + c0 * rand(Nx , Ny) * 0.005;


% assign to the medium inputs
medium.sound_speed = sound_speed_map;
medium.density = density_map;

figure; imagesc(xvec * 1e3 , zvec * 1e3 , sound_speed_map);
set(gca,'DataAspectRatio',[1 1 1])
filename = './phantom/phantom_info_r3mm_is1450_bg1540_depth15_x-5.mat';
% filename = './phantom/homo_with_target.mat';
save(filename);