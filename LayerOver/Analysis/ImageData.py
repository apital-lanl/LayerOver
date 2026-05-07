"""
Copyright 2026. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2026-05-04
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Module for image analysis functions and standards.

"""

import matplotlib.pyplot as plt
import numpy as np
import os



def dynamic_threshold(array, num_bins = 100,
                        show_threshold_graph = True,
                        name = '', threshold = True,
                        calculation_range = 'below',
                        calculation_type = 'outside_CLT',
                        fix_bins = False):
        
    ''' 
    Description: General-purpose function for taking a bin::cnts histogram, assuming central limit theorem, and 
        and returning associated values. If thresholding and calculation values are set, do those and return.
        Taken from SEAM library 2026-05-04.

    INPUT:
        lorem_ipsum     lorem
    ACTION:

    OUTPUT:
        return_dict     dict; contains the following keys:
            'name'                      description
            'cnts'                      list of counts for each bin    
            'bins'                      list of bin edges
            'max_middle_cnt'            maximum count in the middle bins (i.e. not edge bins)
            'bin_size'                  size of the bins (assumes consistent bin sizes)
            'lower_FWHM_bin_idx'        index of the lower bin closest to the FWHM value
            'upper_FWHM_bin_idx'        index of the upper bin closest to the FWHM value
            'FWHM'                      count value at the FWHM (i.e. the minimum of the two closest bins to the half-max value)
            'calculation'               string describing the calculation performed (e.g. 'outside_CLT;below')
            'calculation_cnt_sum'       sum of counts that meet the criteria of the calculation (e.g. sum of counts below the lower threshold)
            'upper_FWHM_pix_value'      pixel value corresponding to the upper FWHM bin
            'peak_pix_value'            pixel value corresponding to the peak bin
            'lower_FWHM_pix_value'      pixel value corresponding to the lower FWHM bin
            'max_threshold_pix_value'   pixel value corresponding to the upper threshold bin
            'min_threshold_pix_value'   pixel value corresponding to the lower threshold bin
    '''

    #Define and initialize variables
    averaging_window = 7

      #Only pull volumes greater than 0 for obvious reasons
    array= np.array(array[array>0])
    average_array = np.convolve(array, np.ones(averaging_window), 'valid') / averaging_window
    array = average_array

        # define pixel integers if 256 is given (strong assumption)
    if num_bins == 256:
        cnts, bins = np.histogram(array, bins= list(range(0,256)))
    else:    
        cnts, bins = np.histogram(array, bins = num_bins)
          
        #check to see if 100 bins is overkill; only applicable in small systems (vignettes from images, for example)
    zero_count = 0
    for cnt in cnts:
        #flag negative values as well for future normalization
        if cnt == 0:
            zero_count += 1
        #reduce the number of bins by the number of zero_counts
        #NOTE: hard-coded 10% below can be adjusted
    if (zero_count >= (0.1*num_bins)) and (not fix_bins):
        zero_int = int(num_bins-zero_count)  #force 'int' in case 'num_bins' isn't
        print(f"Too many zero-cnts with {num_bins}.")
        print(f"\t Reducting to {zero_int} histogram bins.")
        num_bins = zero_int
        cnts, bins = np.histogram(array, bins = num_bins)
            
        
        #get summary variables
    bin_middle_idx = num_bins//2
    bin_size = bins[bin_middle_idx+1]-bins[bin_middle_idx]   #Pick the middle bins to accomodate any edge-bin wierdness; should be consistent bin sizes
    global_max_cnt = max(cnts)
    cnt_max = max(cnts[2:bins.shape[-0]-3])   #skip first and last 2 bins to avoid edge effects (i.e. lots of zeros or )
        
    # check if 'name' is a filename and pull basename if so
    if name != '':
        if os.path.isfile(name):
            name = os.path.basename(name)
    else:
        name = "GenericName"
        
    #ASSUME CLT and get FWHM
    max_idx = np.where(cnts== cnt_max)[0] +1
    if type(max_idx)!=int:
        max_idx = max_idx[0]
    half_max = cnt_max//2
    diff_array = np.array([abs(value-half_max) for value in cnts])
    sorted_diff_array = np.sort(diff_array)
        #TODO: verify CLT by checking Gaussian shape; if not CLT, do something

    #Look at the first 10 values closest to the half-max value and try to pull the indices
    upper_halfmax_idx = num_bins-2   #start at the top end of the histogram
    upper_diff = abs(upper_halfmax_idx-max_idx)  #measure idx distance to upper end of array
    lower_halfmax_idx = -1   #start at bottom end of the histogram; offset by 1 to account for bin indexing
    lower_diff = abs(lower_halfmax_idx-max_idx)  #measure idx distance to lower end of array
        
    #Step through closest values to FWHM expectations and save them
    for value in sorted_diff_array[0:10]:
        #if np.where() returns more than one value, only take the first
        this_idx = np.where(diff_array==value)[0]
            #make sure 'this_idx' isn't an array
        if type(this_idx)!=int:
            this_idx = this_idx[0]
        this_diff = abs(this_idx-max_idx)
            
        #Check if closer upper bound
        if (this_idx > max_idx):
            if this_diff < upper_diff:
                upper_diff = this_diff
                upper_halfmax_idx = this_idx
                    
        #Check if closer lower bound
        if (this_idx < max_idx):
            if this_diff < lower_diff:
                lower_diff = this_diff
                lower_halfmax_idx = this_idx
                    
    #Make sure bounds are in range
    if lower_halfmax_idx < 0:
        lower_halfmax_idx = 0 
    if upper_halfmax_idx > (len(cnts)-2):
        upper_halfmax_idx = len(cnts)-2
        
    #Clean up and generate some summary variables
    upper_halfmax_idx +=1
    lower_halfmax_idx +=1
        #assume FWHM ~= to 2.4 std. devs. and try to find the upper and lower points 
    bin_2stddev_idx_diff = 2* max( abs(max_idx- lower_halfmax_idx), abs(upper_halfmax_idx-max_idx))
    lower_threshold_idx = int(max_idx-bin_2stddev_idx_diff)
    if lower_threshold_idx < 0:
        lower_threshold_idx = 0
    lower_threshold = bins[lower_threshold_idx]
    upper_threshold_idx = int(max_idx+bin_2stddev_idx_diff)
    if upper_threshold_idx > (len(cnts)-1):
        upper_threshold_idx = len(cnts)-1
    upper_threshold = bins[upper_threshold_idx]
    intial_lower_threshold = lower_threshold
        
    # print(f"lower_halfmax_idx: {lower_halfmax_idx}")
    # print(f"upper_halfmax_idx: {upper_halfmax_idx}")
        
    cnts_FWHM = max(cnts[lower_halfmax_idx], cnts[upper_halfmax_idx])- \
        abs(cnts[lower_halfmax_idx] - cnts[upper_halfmax_idx])
    bin_FWHM = max(bins[lower_halfmax_idx], bins[upper_halfmax_idx])- \
        abs(bins[lower_halfmax_idx] - bins[upper_halfmax_idx])
            
        #get variance of 'cnts' values, assume Gaussian, and set 5 std. dev. threshold
    if calculation_type == 'outside_CLT':
        #     #set minimum threshold as 1/8 of the FWHM; 
        #     #TODO: make this less weird and hard-coded; pick a more rational minimum
        # max_cnt_threshold = round(cnt_max/16)
        #Hard-code a threshold for cnts; choosing 25 based on testing
        max_cnt_threshold = 25

        if calculation_range == 'below':
                
            #Make sure the threshold is below about 1000 cnts
            if cnts[lower_threshold_idx] > max_cnt_threshold:
                iteration_cnt = 0
                max_iterations = 10
                while (iteration_cnt < max_iterations) and (cnts[lower_threshold_idx] > max_cnt_threshold):
                    lower_threshold_idx -= 1
                    lower_threshold = bins[lower_threshold_idx]
            #Count the thresholds      
            calc_sum = 0
            for this_bin, this_cnt in zip(bins[1:lower_threshold_idx], cnts[0:lower_threshold_idx]):
                if this_cnt <= max_cnt_threshold:
                    calc_sum += this_cnt
            
        if calculation_range == 'above':
            #Make sure the threshold is below about 1000 cnts
            if cnts[upper_threshold_idx] > max_cnt_threshold:
                iteration_cnt = 0
                max_iterations = 10
                while (iteration_cnt < max_iterations) and (cnts[upper_threshold_idx] > max_cnt_threshold):
                    upper_threshold_idx += 1
                    upper_threshold = bins[upper_threshold_idx]
            #Count the thresholds      
            calc_sum = 0
            for this_bin, this_cnt in zip(bins[upper_threshold_idx:-1], cnts[(upper_threshold_idx-1):-1]):
                if this_cnt <= max_cnt_threshold:
                    calc_sum += this_cnt 
        
    #Plot everything if called for
    if show_threshold_graph:
            #plot the peak max and pseudo FWHM
        plt.plot([bins[lower_halfmax_idx], bins[upper_halfmax_idx]],[half_max, half_max], color = 'm', linewidth = 4)
        plt.plot([bins[max_idx], bins[max_idx]], [0, cnt_max], color = 'm', linewidth = 1)
            #show the actual closest-bin level to the FWHM that the psuedo is calculated
        plt.plot([bins[lower_halfmax_idx], bins[upper_halfmax_idx]],[cnts_FWHM, cnts_FWHM], color = 'purple', linewidth = 5)
            #plot the upper and lower 5-sigma threshholds
        plt.plot([intial_lower_threshold, intial_lower_threshold], [0, half_max], color = 'gray', linewidth = 2)
        plt.plot([upper_threshold, upper_threshold], [0, half_max], color = 'lime', linewidth = 2, alpha = 0.6)
        plt.plot([lower_threshold, lower_threshold], [0, half_max], color = 'red', linewidth = 2, alpha = 0.6)
            #finish plotting the histogram values
        plt.plot(bins[1::], cnts, linewidth = 1)
        plt.scatter(bins[1::], cnts)
        plt.title(f"Height hist- {name}")
        plt.legend(['FWHM', 'Peak', 'Closest FWHM Bin', 'Initial Lower Threshold', 'Upper Threshold', 'Lower Threshold'])
        plt.show()
            
    return_dict = {
        'name': name,
        'cnts': cnts.tolist(),
        'bins': bins.tolist(),
        'max_middle_cnt': cnt_max,
        'bin_size': bin_size,
        'lower_FWHM_bin_idx': lower_halfmax_idx,
        'upper_FWHM_bin_idx': upper_halfmax_idx,
        'FWHM': cnts_FWHM,
        'calculation': calculation_type + ';' + calculation_range,
        'calculation_cnt_sum': calc_sum,
        'upper_FWHM_pix_value': bins[upper_halfmax_idx],
        'peak_pix_value': min(bins[lower_halfmax_idx], bins[upper_halfmax_idx]) + abs(bins[upper_halfmax_idx]-bins[lower_halfmax_idx]),
        'lower_FWHM_pix_value': bins[lower_halfmax_idx],
        'max_threshold_pix_value': upper_threshold,
        'min_threshold_pix_value': lower_threshold
        }
        
    return return_dict