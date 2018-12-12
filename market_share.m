    clc
	clearvars
    close all
   
    %%  Determine input
%   Select input file and sheet
    filn        =   [pwd '/AE4423_Datasheets.xlsx'];
    filn2       =   [pwd '/Group8_results.xlsx'];
    
    frequencies = xlsread(filn2,'Group8-data','A1:X24');
    frequencies_c = xlsread(filn,'Group 8','C89:Z112');
    
    demand_low = xlsread(filn,'Group 8','C63:Z86');
    demand_high = xlsread(filn,'Group 8','C37:Z60');
    
    a = 1.0;
    b = 1.7;
    
    %% Market Share
    demand_low_new = zeros(length(frequencies),length(frequencies'));
    demand_high_new = zeros(length(frequencies),length(frequencies'));
    for i = 1:length(frequencies)
        for j = 1:length(frequencies')
            freq_d = frequencies(i,j);
            freq_i = min(frequencies(3,j),frequencies(i,3));
            freq_c = frequencies_c(i,j);
            ms = (freq_d^a + freq_i^b)/(freq_d^a + freq_i^b + freq_c^a + 1e-15);
            d_lo = ms*demand_low(i,j);
            d_hi = ms*demand_high(i,j);
            demand_low_new(i,j) = d_lo;
            demand_high_new(i,j) = d_hi;
        end
    end
    demand_new = [demand_low_new;demand_high_new];
    xlswrite(filn2,demand_new,'new_demands')