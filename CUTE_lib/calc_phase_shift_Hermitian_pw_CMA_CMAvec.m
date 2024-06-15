function[PS , mask] = calc_phase_shift_Hermitian_pw_CMA_CMAvec(xvec , zvec , c , rf_frame1 , ...
    rf_frame2 , axiscorrx , CMAvec)
rf_signal1 = rf_frame1.rf_data;
t_array = rf_frame1.t_array;
phi1 = rf_frame1.phi;
[I1 , mask1] = mk_pw_rf_CMA_image_CMAvec(xvec , zvec , c , phi1 , rf_signal1 , ...
    t_array - 2/5e6, axiscorrx , CMAvec);
nanind = isnan(I1);
I1(nanind) = 0;
inan = isnan(I1);
I1(inan) = 0;

rf_signal2 = rf_frame2.rf_data;
phi2 = rf_frame2.phi;
[I2 , mask2] = mk_pw_rf_CMA_image_CMAvec(xvec , zvec , c , phi2 , rf_signal2 , ...
    t_array - 2/5e6, axiscorrx , CMAvec);
nanind = isnan(I2);
I2(nanind) = 0;
inan = isnan(I2);
I2(inan) = 0;

DeltaTheta = I1.* conj(I2);
PS = angle(DeltaTheta) / (2 * pi * 5e6);
for cc = 1 : length(CMAvec)
    PS(cc , : , :) = medfilt2(squeeze(PS(cc , : , :)));
end
mask = mask1 | mask2;





