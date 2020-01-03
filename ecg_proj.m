function ecg_proj
    
    %ECG data file path
    file_ecg = '/Users/pproctor/Documents/PSU/ECE522/ecg_project_loc/custom_ecg.mat';
    file_t = '/Users/pproctor/Documents/PSU/ECE522/ecg_project_loc/custom_ecg_x.mat';
    imag_p = '/Users/pproctor/Documents/MATLAB/mathworks-2/matlab_codes/ECE_522_project/report_images';

    file_ecg_ex = '/ecg_example';
    type = 0; %1 default, 2 fast, 0 custom
    
    %Prepare the data 
    [ecg_data,x_data] = init_data(file_ecg,file_t,type,imag_p,file_ecg_ex);
    
    %Component values from circuit schematic
    r1 = 47e3;
    r2 = 1e6;
    c1 = .33e-6;
    c2 = 22e-9;
    comps = [r1 r2 c1 c2];
    Fs = 950;
    bins = 1024;
    title_d = 'Raw filter';

    %Freq. Spectrum with DC offset removed
    dft((ecg_data - mean(ecg_data)),Fs, bins);

    %Convert CTS transfer func. to DT transfer func. and perform filtering
    [Hd,num_d,den_d] = init_filt(comps,Fs,imag_p); 
    [filter_d,filtfilt_d] = filt_op(num_d,den_d,ecg_data,x_data,title_d,'/raw_filt',imag_p);
    [ecg_data_a,filter_d_a,D] = alignsignals(ecg_data,filter_d);
    
    %Calculate the RMSE and DFT for uncorrected band-pass filter
    H_i_filt_rmse = immse(ecg_data,filter_d);
    H_i_filt_rmse_a = immse(ecg_data_a(1:(end-D)),filter_d_a);
    H_i_filt_filt_rmse = immse(ecg_data,filtfilt_d);
    dft((filter_d-mean(filter_d)),Fs,bins);

    %Set the resolution and range of pole zero grid in polar coordinates
    idx_r_lb = .05;
    idx_r_ub = 1;
    r_step = .01;
    idx_phi_lb = 0;
    idx_phi_ub = pi;
    phi_step = .05;
    idx_t = 0;
    run = 0;

    [gd_i,f_i] = group_delay(num_d,den_d,Fs);
     
   close all;
    
   %Set number of iterations and threshold for PZ placement
   thresh = .05;
   iter = 1;
   tot_iter = 1;
   He = 1;
   gd_mat = zeros(tot_iter, length(f_i));
   gd_mean_track = zeros(tot_iter,1);
   
   %Optimization loop
   while (run < tot_iter && iter)
       %Col are phi, row are r
        gp_delay_st = zeros(ceil((idx_r_ub-idx_r_lb)/r_step),ceil((idx_phi_ub-idx_phi_lb)/phi_step));
        r_n = 0;
        
        idx_t = idx_t + 1;
        for r = idx_r_lb:r_step:idx_r_ub
            r_n = r_n + 1;
            phi_n = 1;
            for phi = idx_phi_lb:phi_step:idx_phi_ub
                H_ap = gen_AP(r,phi,Fs);
                H_temp = H_ap * Hd;
                [num_t,den_t] = tfdata(H_temp,'v');
                [delay,w_out] = group_delay(num_t,den_t,Fs);
                delay_diff = diff(delay(1: end));
                delay_mean = harmmean(abs(delay_diff(~isinf(delay_diff))));
                gp_delay_st(r_n,phi_n) = delay_mean;
                phi_n = phi_n + 1;

            end

        end

        %Determine all-pass that minimizes group delay slope
        [H_ap_min,num_min,den_min] = AP_min(r_step,phi_step,gp_delay_st,Fs);
        H_cas = Hd * H_ap_min;
        [num_cas,den_cas] = tfdata(H_cas,'v');
        [gd_cas,f_cas] = group_delay(num_cas,den_cas,Fs);
        gd_cas_diff = diff(gd_cas(1:end));
        gd_cas_mean = harmmean(abs(gd_cas_diff(~isinf(gd_cas_diff))));
        disp('Gd_cas_mean: ') 
        disp(gd_cas_mean);
    
        %Only add the AP to the system if it reduces group delay slope
        if gd_cas_mean <= thresh
            gd_mat(run+1,:) = gd_cas;
            gd_mean_track(run+1) = gd_cas_mean;
            Hd = Hd * H_ap_min;
            He = He * H_ap_min;
            thresh = gd_cas_mean;
        end
        
        run = run + 1;
        disp(run);
   end
    
    %Group delay for AP system(He) and AP cascaded with BP(Hd)
    [num_e, den_e] = tfdata(He,'v');
    [num_2,den_2] = tfdata(Hd,'v');
    [gd_e,f_win_e] = group_delay(num_e,den_e,Fs);
    [gd,f_win] = group_delay(num_2,den_2,Fs);
    gd_2nd = diff(gd);

    [row,col] = size(gd_mat);
    
    %Filter with AP cascaded with BP and calc. RMSE
    x_data=0.01:0.01:1.2;
    ecg_data = [ecg_data zeros(1,length(x_data) - length(ecg_data))];
    [filter_d_cas,filtfilt_cas_d] = filt_op(num_2,den_2,ecg_data,x_data,title_d,'/hd_filt',imag_p);
    [ecg_data_a,filter_d_cas_a,D_a] = alignsignals(ecg_data,filter_d_cas);
    H_d_filt_rmse = immse(ecg_data,filter_d_cas);
    H_d_filt_rmse_a = immse(ecg_data_a(1:(end-D_a)),filter_d_cas_a);
    H_d_filt_filt_rmse = immse(ecg_data,filtfilt_cas_d);
    



end

function [Hd,num_d,den_d] = init_filt(params,Fs,imag_p)
    %Initialize the CTS band-pass filter
    
    r1 = params(1);
    r2 = params(2);
    c1 = params(3);
    c2 = params(4);
    
    num = [c1*r2 0];
    den = [c1*c2*r1*r2 c1*r2+c2*r1-c1*r1 1];
    
    [gd_o,f_o] = group_delay(num,den,Fs);
    
    [Hd,num_d,den_d] = initialize(num,den,Fs,imag_p);
    
end

function [gd,f]= group_delay(num,den,fs)
    %Calculate the group delay over a specific freq. range
    
    fin = 5:.0475:(195-.0475);
    [gd,f] = grpdelay(num,den,fin,fs);
end

function [init_H,num_i,den_i] = initialize(num,den,fs,file_p)
    %Convert the CTS. transfer function to DT via ZOH.
    
    Hc=tf(num,den);
    figure;
    impulse(Hc);
    init_H = c2d(Hc,1/fs,'zoh');
    [num_i,den_i] = tfdata(init_H,'v');
    bode_config(init_H,fs,file_p,'/bode_mag'); 
    figure;
    pzmap(init_H)
    
end

function [H_temp] = gen_AP(m,phi,fs)
    %Generate all pass poles and zeros over specified range
    
    syms k
    z = m*exp(1j*phi);
    z_con = conj(z);
    p = m*exp(1j*phi);
    p_con = conj(p);
    
    num_coeff = sym2poly(expand((1-z*k)*(1-z_con*k)));
    den_coeff = sym2poly(expand((k-p)*(k - p_con)));
    
    H_temp = tf(num_coeff, den_coeff,1/fs);
end

function [ecg_load,x_load] = init_data(file_name_1,file_name_2,type,file_p,file_n)
    %Load the initial ECG data depending on the desired heart rate
    
    data_ecg = load(file_name_1);
    x_ecg = load(file_name_2);
    
    if type == 2
        ecg_load = data_ecg.ecg(33:102);
        x_load = x_ecg.x(33:102);
    elseif type == 1
        ecg_load = data_ecg.ecg(51:141);
        x_load = x_ecg.x(51:141);
    elseif type == 0
        ecg_load = data_ecg.ecg(51:141);
        x_load = x_ecg.x(51:141);
    end
    
    figure;
    %plot(x_load,ecg_load);
    title('Simulated ECG Data')
    xlabel('Time[s]');
    ylabel('Amplitude [mV]');
    hold on;
    ecg_load = ecg_load - ecg_load(1);
    plot(x_load,ecg_load);
    
    hold off;
end

function dft(input_d,Fs,bins)
    %Zero pad the input signal and take the FFT
    
    zero_pad = 4^ceil(log2(length(input_d))) - length(input_d);
    input_d = [input_d zeros(1,zero_pad)];
    f_d = fft(input_d,bins);
    f_ax = (1:1:bins/2)*Fs/bins;
    figure;
    plot(f_ax,abs(f_d(1:(bins/2))));
    xlabel('Frequency [Hz]');
    ylabel('|H(f)|')
    
    
end

function [filt_data,filt_2_data] = filt_op(num,den,data,x,title_n,file_n,file_p)
    %Filtering the data with two different methods
    
    filt_data = filter(num,den,data);
    filt_2_data = filtfilt(num,den,data); 
    hold off;
    
end

function [H_min_temp,num_min,den_min] = AP_min(r_st,phi_st,results_mat,fs)
    %Determine the all pass pole zero placement that
    %makes the group delay more constant
    
    [row_min,col_min] = find(results_mat == min(results_mat(:)));
    r_min = r_st*row_min;
    phi_min = phi_st*col_min;
    H_min_temp = gen_AP(r_min, phi_min,fs);
    [num_min,den_min] = tfdata(H_min_temp,'v'); 
    
end

function bode_config(H_f,fs,file_p,file_n)
    %Configuring Bode Plots

    figure;
    opts = bodeoptions('cstprefs');
    opts.FreqUnits = 'Hz';
    opts.FreqScale = 'Linear';
    opts.magVisible = 'off';
    opts.Title.String = 'Phase Resp.';
    opts.Xlim = [.1,700];
    h = bodeplot(H_f,opts);
        
end







