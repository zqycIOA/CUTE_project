function[PS] = calc_phase_shift_Hermitian_dw(xvec , zvec , c , rf_frame1 , ...
    rf_frame2 , axiscortx , axiscorrx)
rf_signal1 = rf_frame1.rf_data;
t_array = rf_frame1.t_array;
emitx1 = axiscortx(rf_frame1.emitele);
I1 = mk_dw_rf_image(xvec , zvec , c , emitx1 , rf_signal1 , ...
    t_array, axiscorrx);
nanind = isnan(I1);
I1(nanind) = 0;
inan = isnan(I1);
I1(inan) = 0;

rf_signal2 = rf_frame2.rf_data;
emitx2 = axiscortx(rf_frame2.emitele);
I2 = mk_dw_rf_image(xvec , zvec , c , emitx2 , rf_signal2 , ...
    t_array, axiscorrx);
nanind = isnan(I2);
I2(nanind) = 0;
inan = isnan(I2);
I2(inan) = 0;

DeltaTheta = I1.* conj(I2);
PS = medfilt2(angle(DeltaTheta)) / (2 * pi * 5e6);





