filename1 = '../field_signal/r3mm_is1450_bg1540_depth15_x0_pw_255rx_phi_0.mat';
% filename1 = '../field_signal/homo_pw_255rx_phi_0.mat';
load(filename1);
rf_data1 = rf_frame.rf_data;

filename2 = '../field_signal/r3mm_is1450_bg1540_depth15_x0_pw_255rx_phi_5.mat';
% filename2 = '../field_signal/homo_pw_255rx_phi_5.mat';
load(filename2);
rf_data2 = rf_frame.rf_data;

phi1 = 0 * pi / 180;
phi2 = 5 * pi / 180;

t_array = rf_frame.t_array;
fs = 1 / (t_array(2) - t_array(1));

c = 1540;
[sensor_xgrid , senzor_zgrid] = meshgrid(sensor_xvec , sensor_zvec);
sensor_xs = sensor_xgrid(:);
sensor_zs = senzor_zgrid(:);

totaldelay = zeros(length(sensor_xs) , 1);
tof_est1 = (sensor_xs .* sin(phi1) + sensor_zs .* cos(phi1)) / c;
tof_est2 = (sensor_xs .* sin(phi2) + sensor_zs .* cos(phi2)) / c;

for sc = 1 : length(sensor_zs)
    tof_temp1 = tof_est1(sc);
    tof_temp2 = tof_est2(sc);
    seg_length = 50 / 5e6; seg_vec = round(-seg_length * fs) : round(seg_length * fs);
    [~ , sign1] = min(abs(t_array - tof_temp1)); [~ , sign2] = min(abs(t_array - tof_temp2));
    sig_section1 = t_array * 0; sig_section1(sign1 + seg_vec) = 1;
    sig_section2 = t_array * 0; sig_section2(sign2 + seg_vec) = 1;
    sig1 = rf_data1(sc , :); sig1 = sig1(sig_section1 == 1);
    sig2 = rf_data2(sc , :); sig2 = sig2(sig_section2 == 1);
    totaldelay(sc) = calc_relative_delay(sig1 , sig2 , fs);
end

totaldelay1 = reshape(totaldelay , size(sensor_xgrid));

figure; imagesc(sensor_xvec * 1e3 , sensor_zvec * 1e3 , totaldelay1);