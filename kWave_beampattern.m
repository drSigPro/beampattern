%-- Script for estimating pressure pattern in orthogonal planes (xy, xz and
% yz) for various array configurations (e.g. L, T or Cross shaped)
%-- Author: Mahesh R Panicker (mahesh.signalproc@gmail.com) and Gayathri M
%-- Version: 1.0
%-- Date: July, 2024
%-------------------------------------------------------------------------%

%% Clear Screen and Variables
clc;
clearvars;

%% Generic Parameters
GPU                 = false;        % default running on CPU
if (gpuDeviceCount)
    GPU             = true;         % set to true to simulate using GPU, to false to simulate in CPU
end
USE_CPU_C_CODE      = false;        % set to true to simulate using CPU cores in the absence of CUDA
SimPlot             = false;        % set to true to view wave propagation.Will be active only if USE_CPU_CORES == false
disp_source         = 0;            %0= don't display the source with medium
MASK_PLANE          = 'xy';         % set to 'xy' or 'xz' or 'yz' to generate the beam pattern in different planes z is the depth direction, x is probe direction

%% Control Parameters
endDepth = 0.15;                    % Maximum depth in axial and lateral direction [m]
latDepth = 0.15;
source_freq = 725e3;                % Transmit frequency [Hz]
max_freq = 1.5*source_freq;         % Maximum supported frequency (~1.2 x desired transmit)
c_ave = 1500;                       % Average Speed of Sound
lamda = c_ave/max_freq;             % Estimate Wavelength

%% DEFINE THE K-WAVE GRID
%Important note: In k-wave x is the depth direction, y is the sensor (or lateral) direction and z is the elevation
% Calculate the spacing between the grid points
dx = lamda/2;                       % Assuming spacing[m]
dy = dx;
dz = dx;

% Number of grid points
Nx = ceil(1.2*endDepth/dx);
Ny = ceil(1.2*latDepth/dy);
Nz = ceil(1.2*latDepth/dz);

% create the k-space grid (2D grid here)
kgrid = kWaveGrid(Nx, dx, Ny,dy, Nz, dz);

PML = getOptimalPMLSize(kgrid);
PML_X_size = PML(1);
PML_Y_size = PML(2);
PML_Z_size = PML(3);


%% DEFINE THE MEDIUM PARAMETERS
background_map_mean = 1;    % To ensure mean speed of ~1540 m/s             
background_map_std = 0.001; % Change the standard deviation if you want a very hetrogenous medium
background_map = background_map_mean + background_map_std * randn([Nx, Ny, Nz]);
medium.sound_speed = 1540.*ones(Nx, Ny, Nz).*background_map;      % [m/s]
medium.density = 1000.*ones(Nx, Ny, Nz).*background_map;          % [kg/m^3]
medium.alpha_coeff = 0.05;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.05;
% medium.BonA = 6;


%% Create the time vector for wave evaluation
t_end = (Nx * dx) * 2.2 / c_ave;   % [s]
CFL = [];
kgrid.t_array = makeTime(kgrid, medium.sound_speed, CFL, t_end);


%% Creating the Transducer array
%% Transmit Array with single element transmit
% physical properties of the transducer
tr_element_width = 12*dx;  % width of each element [grid points]
tr_element_length = 12*dx;
tr_rotation        = [0,0,0];

% Create empty kWaveArray
karray_tr = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

element_spacing = 1; %Change spacing to have the array, '1' will be a continous array
for i = round(Nx/4):element_spacing:round(3*Nx/4) %Change this to change the start and end location of arrays
    % L shaped
    pos_element_x= [kgrid.x_vec(1), kgrid.y_vec(round(Ny/4)), kgrid.z_vec(i)];
    pos_element_y= [kgrid.x_vec(1), kgrid.y_vec(i), kgrid.z_vec(round(Nz/4))];

    % T shaped
    % pos_element_x= [kgrid.x_vec(1), kgrid.y_vec(round(Nx/2)), kgrid.z_vec(i)];
    % pos_element_y= [kgrid.x_vec(1), kgrid.y_vec(i), kgrid.z_vec(round(Nz/4))];
    
    % Cross shaped
    % pos_element_x= [kgrid.x_vec(round(0.05*Nx)), kgrid.y_vec(round(Ny/2)), kgrid.z_vec(i)];
    % pos_element_y= [kgrid.x_vec(round(0.05*Nx)), kgrid.y_vec(i), kgrid.z_vec(round(Nz/2))];
    karray_tr.addRectElement(pos_element_x,tr_element_width,tr_element_width, [0,0,0]);
    karray_tr.addRectElement(pos_element_y,tr_element_width,tr_element_width, [45,0,0]);
end

karray_tr.setArrayPosition([0, 0, 0], [0,0,0]);


%% Creating source from the transmit element
source.p_mask=zeros(Nx, Ny, Nz);
source.p_mask = karray_tr.getArrayBinaryMask(kgrid);
probe_geometry_tr = karray_tr.getElementPositions;

% Display source
if (disp_source)
    F_M=(source.p_mask *1000 + medium.sound_speed);
    volumeViewer(F_M);
end

%% DEFINE THE INPUT SIGNAL
% define properties of the input signal
source_strength     = 4e6;      % [Pa], Maximum is 4e6
source_cycles       = 5;        % number of tone burst cycles
source.p = source_strength*toneBurst(1/kgrid.dt, source_freq, source_cycles); %Create a toneBurst signal
source.p_mode = 'dirichlet';

%% DEFINE SENSOR MASK
% define a sensor mask through the central plane
sensor.mask = zeros(Nx, Ny, Nz);
switch MASK_PLANE
    case 'xy'

        % define mask
        %sensor.mask = karray_sensor_xy.getArrayBinaryMask(kgrid);
        sensor.mask(:, :, round(Nz/2)) = 1;

        % store y axis properties
        Ni = Nx;
        i_label = 'x';
        i_vec = kgrid.x_vec;
        Nj = Ny;
        j_vec = kgrid.y_vec;
        j_label = 'y';

    case 'xz'

        % define mask
        %sensor.mask = karray_sensor_xy.getArrayBinaryMask(kgrid);
        sensor.mask(:, round(Ny/2), :) = 1;

        % store z axis properties
        Ni = Nx;
        i_label = 'x';
        i_vec = kgrid.x_vec;
        Nj = Nz;
        j_vec = kgrid.z_vec;
        j_label = 'z';

    case 'yz'

        % define mask
        %sensor.mask = karray_sensor_xy.getArrayBinaryMask(kgrid);
        sensor.mask(round(Nx/2),: , :) = 1;

        % store z axis properties
        Ni = Ny;
        i_vec = kgrid.y_vec;
        i_label = 'y';
        Nj = Nz;
        j_vec = kgrid.z_vec;
        j_label = 'z';

end

% set the record mode such that only the rms and peak values are stored
sensor.record = {'p_rms','p_final'};
display_mask = source.p_mask;

%% RUN THE SIMULATION
% set the input settings
if(GPU)
    DATA_CAST       = 'gpuArray-single';
else
    DATA_CAST       = 'single';
end
input_args = {'DisplayMask', display_mask, ...
    'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_size, PML_Y_size, PML_Z_size], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotScale', [-1, 1] * source_strength};
% run the simulation
sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});

%% COMPUTE THE BEAM PATTERN USING SIMULATION STATISTICS
% reshape the returned rms and max fields to their original position
sensor_data.p_rms = reshape(sensor_data.p_rms, [Ni, Nj]);
p_rms = sensor_data.p_rms;
start_point = 1;

% plot the beam pattern using the pressure rms
figure;
imagesc(j_vec * 1e3, (i_vec(start_point:end) - min(i_vec(:))) * 1e3, p_rms(start_point:end,:)* 1e-6);
xlabel([j_label '-position [mm]']);
ylabel([i_label '-position [mm]']);
title('Total Beam Pattern Using RMS Of Recorded Pressure');
colormap(jet(256));
c = colorbar;
ylabel(c, 'Pressure [MPa]');
axis image;

