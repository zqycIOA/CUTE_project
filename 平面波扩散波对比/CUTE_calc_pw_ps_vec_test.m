clear;
% close all;
fs = 50e6;
sigmaxvec = (-20 : 0.5 : 20) * 1e-3; sigmazvec = (5 : 0.5 : 40) * 1e-3;
[sigmaxs , sigmazs] = meshgrid(sigmaxvec , sigmazvec);
pitch = 0.300e-3;
N_elements = 128;
axiscor = linspace(-0.5, 0.5, N_elements)*pitch*(N_elements-1);

axiscortx = (1 : N_elements) * pitch;
axiscortx = axiscortx - mean(axiscortx);
arrayparam.axiscortx = axiscortx - mean(axiscortx);
arrayparam.axiscorrx = axiscortx(1):pitch/2:axiscortx(end);

ROIparam.xvec = sigmaxvec;
ROIparam.zvec = sigmazvec;
ROIparam.c = 1540;

psparam.CMAvec = [-15 0 15] * pi / 180;
psparam.psstart = -5 * pi / 180;
psparam.psend = 5 * pi / 180;
psparam.psgap = 5 * pi / 180;
psparam.pslength = 1 * pi / 180;

linevec = 1:21;
lc = 1;
for linecount = linevec
    filename = ['./rf_signal/pw_r3mm_x0z15_is1450_bg1540_255rx/line_' num2str(linecount) '.mat'];
    load(filename);
    rf_frame.phi = rf_frame.phi;
    totalframes(lc) = rf_frame;
    clear rf_frame; 
    lc = lc + 1;
end

[totalps , totalmask] = calc_phase_shift_vec_pw_CMA(totalframes , psparam , ROIparam , arrayparam);
% totalps = calc_phase_shift_vec_pw_CMA(totalframes , psparam , ROIparam , arrayparam);

inclusion_map = sqrt(sigmaxs.^2 + (sigmazs - 15e-3).^2) < 3e-3;
ROIparam.Cmap0 = 1540 * (1 - inclusion_map) + 1450 * inclusion_map;
psref = calc_pw_ps_prediction(ROIparam , psparam , arrayparam);

figure;
imagesc(sigmaxvec * 1e3 , sigmazvec * 1e3 , totalps{2,2});
caxis([-2e-7 2e-7]); colorbar
axis([-20 20 0 40])

figure;
imagesc(sigmaxvec * 1e3 , sigmazvec * 1e3 , psref{2,2}.*totalmask{2,2});
% imagesc(sigmaxvec * 1e3 , sigmazvec * 1e3 , psref{2});
caxis([-2e-7 2e-7]); colorbar
axis([-20 20 0 40])

filename = './ps_info/ps_info_pw_r3mm_is1450_bg1540_depth15.mat';
save(filename , 'totalmask' , 'totalps' , 'psref' , 'ROIparam' , 'psparam' , 'arrayparam');