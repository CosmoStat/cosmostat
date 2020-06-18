i=double(imread('lenna512.pgm'));
imax = max(max(i));

% checking energy...
disp(sprintf('Initial energy = %.10g', sum(sum(i.^2)) ));

all_filter_types = {'spline', '7-9'};

for iter_filter_type = 1:2
    filter_type = all_filter_types{iter_filter_type};
    for nb_scales = 1:3;

        btree = cwpt2_btree(nb_scales,1);       % 1 = full
        twp=cwpt2_interface(i, 'forward', filter_type, 'quad', btree);
        rwp=cwpt2_interface(twp, 'inverse', filter_type, 'quad', btree);

        % checking reconstruction...
        disp(sprintf('FULL @ %d relative E_rcn: %.10g', nb_scales, mean(mean(abs(i-rwp)))/imax   ));

        % checking energy...
        if nb_scales == 1
            disp(sprintf('Energy = %.10g',   sum(sum(sum(twp.^2)))/4   ));
                 % correct for nb_scales = 1 only
        end;

        tw=cwpt2_interface(i, 'forward', filter_type, 'quad', nb_scales);
        rw=cwpt2_interface(tw, 'inverse', filter_type, 'quad', 0);

        % checking reconstruction...
        disp(sprintf('WAVELET @ %d relative E_rcn: %.10g', nb_scales, mean(mean(abs(i-rw)))/imax    ));

        % checking energy...
        if nb_scales==1
            disp(sprintf('Energy = %.10g',   sum(sum(sum(tw.^2)))/4   ));
                 % correct for nb_scales = 1 only
        end;

        tw3=cwpt2_interface(i, 'forward', filter_type, 'tri', nb_scales);
        rw3=cwpt2_interface(tw3, 'inverse', filter_type, 'tri', 0);

        % checking reconstruction...
        disp(sprintf('3-WAY @ %d relative E_rcn: %.10g', nb_scales, mean(mean(abs(i-rw3)))/imax    ));

        % checking energy...
        if nb_scales==1
            disp(sprintf('Energy (should be diff) = %.10g',   sum(sum(sum(tw3.^2)))/4   ));
                 % correct for nb_scales = 1 only
        end;

    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    % 1D wavelet transform
    dt1d = cwt(d1d,nb_scales);
    plotzoom(dt1d);
    dr1d = cwti(dt1d);

    disp(max(abs(dr1d-d1d')));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checking energy distribution

if 0

    btree = cwpt2_btree(3,1);       % 1 = full

    for node = double(cwpt2_get_leaf_order(btree, 4))';     % 4 = quad-tree
        figure(1);
        btree_freqzone_fft(node, 1);
        sp = btree_freqzone_fft(node, 2);
        sample_signal = fftshift(real(ifft2(sp)));

        dt = cwpt2(sample_signal,btree);
        dr = cwpt2i(dt,btree);

        disp(max(max(abs(dr-sample_signal))));

        cwpt2_show_energy(dt,btree);
        pause;
    end;

end;


% 2D Continuous Wavelet Packet Transform package
% (c) 2002-2005 Let It Wave, all rights reserved
