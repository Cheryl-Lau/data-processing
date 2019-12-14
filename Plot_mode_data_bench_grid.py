'''
 #################################################################
    Comparing Period Spacings of ADIPLS/GYRE grid and benchmark
 #################################################################

  Reads both ADIPLS/GYRE benchmark and ADIPLS/GYRE grid
  Filter best-matching models with gradient and intercept uncertainties 
  of benchmark
  Extract parameter set of best-matching models
  Plot all best-matching models on a graph with benchmark 

  Evaluate chi2 for all grid models 
  Plot chi2 variation with each parameter combination in the set

'''

import numpy as np
import matplotlib.pyplot as plt 
import os, fnmatch
import pandas as pd



# Set benchmark model directory 
benchmark_work_dir = 'work_bench_adipls'

# Set grid models directory 
grid_work_dir = 'work_gGbA_fine_2' 



# Settings for filtering adipls output -
# Enter azimuthal order m
m = 0  
# Frequency range (uHz) 
nu1 = 3.8
nu2 = 38.6


# Controls for fine grids testing stage -
# used for obtaining parameters which are sufficiently close to benchmark  
# when there is not enough best-matching models 

bench_err_filter = False     # if True, use the original grad_err and intercept_err as filter 
grad_filter = 0.004         # if False, set a new grad_filter
intercept_filter = 0.00578     # and a new intercept_filter

best_matching_filter = False    # if True, require at least 4 models 
num_best_matching = 1          # if False, set a new required number 


# fill this the smallest chi2 value to scale opacity of plots 
plot_scale_factor = 2190.131487




def main(benchmark_work_dir,grid_work_dir):
    
    ### Processing Benchmark - 
    
    if fnmatch.fnmatch(benchmark_work_dir,'work_bench_gyre'):
        
        # call subroutine gyre_p_spacing_import
        gyre_output_filename = fnmatch.filter(os.listdir(benchmark_work_dir), 'gyre_output_*')[0]
        print('---- Reading GYRE benchmark <'+gyre_output_filename+'> ----')
        period_bench,p_spacing_bench = gyre_p_spacing_import(benchmark_work_dir,gyre_output_filename, \
            'profile_puls_summary.txt')
            
        # call subroutine gradient
        grad_bench,grad_err_bench,intercept_bench,intercept_err_bench = gradient(period_bench,p_spacing_bench) 
        
        # extract benchmark label from filename
        label_bench = gyre_output_filename.replace('gyre_output_','')
        
    elif fnmatch.fnmatch(benchmark_work_dir,'work_bench_adipls'):
        
        adipls_filename = fnmatch.filter(os.listdir(benchmark_work_dir), 'save_*_mode.data')[0]
        print('---- Reading ADIPLS benchmark <'+adipls_filename+'> ----')
        period_bench,p_spacing_bench = adipls_p_spacing_import(benchmark_work_dir,adipls_filename)
        grad_bench,grad_err_bench,intercept_bench,intercept_err_bench = gradient(period_bench,p_spacing_bench)
        
        label = adipls_filename.replace('save_','')
        label_bench = label.replace('_mode.data','')
     
    # when adjusting fine grids, aim to get user-defined filters to within this error 
    print('Benchmark gradient:           ',grad_bench)
    print('Benchmark grad error:         ',grad_err_bench)
    print('Benchmark intercept error:    ',intercept_err_bench)
     
    ### Processing Grid -
       
    if fnmatch.fnmatch(grid_work_dir,'work_grid_gyre_*') or fnmatch.fnmatch(grid_work_dir,'work_gG*'):
        # call subroutine to loop over whole grid, get param_set with chi2 for the whole grid
        # and the best-matching models 
        print('---- Begin GYRE grid processing ----')
        chi2_list,M_list,O_list,H_list,R_list,chi2_list_best,M_list_best,O_list_best, \
            H_list_best,R_list_best = gyre_grid_scan(grid_work_dir,period_bench, \
            p_spacing_bench,grad_bench,grad_err_bench,intercept_bench,intercept_err_bench) 
              
    elif fnmatch.fnmatch(grid_work_dir,'work_grid_adipls_*') or fnmatch.fnmatch(grid_work_dir,'work_gA*'):
        print('---- Begin ADIPLS grid processing ----')
        chi2_list,M_list,O_list,H_list,R_list,chi2_list_best,M_list_best,O_list_best, \
            H_list_best,R_list_best = adipls_grid_scan(grid_work_dir,period_bench, \
            p_spacing_bench,grad_bench,grad_err_bench,intercept_bench,intercept_err_bench)     
            
            
    # call subroutine param_correlation to generate Contour plots for the whole grid 
    param_correlation(chi2_list,M_list,O_list,H_list,R_list,label_bench)
    
    # sort outputs in chi2 ascending order 
    param_set = zip(chi2_list,M_list,O_list,H_list,R_list)
    sort_param = sorted(param_set)
    chi2_list,M_list,O_list,H_list,R_list = zip(*sort_param)
    
    # if best-matching outputs are not empty 
    if chi2_list_best:
        best_set = zip(chi2_list_best,M_list_best,O_list_best,H_list_best,R_list_best)
        sort_best = sorted(best_set)
        chi2_list_best,M_list_best,O_list_best,H_list_best,R_list_best = zip(*sort_best) 
        
        # Proceed to output best-matching models results if at least 4 models have passed through filter
        if best_matching_filter:
            if len(chi2_list_best) >= 4: 
                print('Parameter sets of the best-matching models -')
                outputs = pd.DataFrame({'M': M_list_best,'O': O_list_best,'H': H_list_best, \
                                        'R': R_list_best,'chi2': chi2_list_best})
                print(outputs[['M','O','H','R','chi2']])
            else:
                print('Insufficient best-matching models')
        
        # For temporarily changing the filters to make fine grids 
        elif not best_matching_filter and num_best_matching:
            if len(chi2_list_best) >= num_best_matching: 
                print('Parameter sets of the models within the given gradient and intercept range -')
                outputs = pd.DataFrame({'M': M_list_best,'O': O_list_best,'H': H_list_best, \
                                        'R': R_list_best,'chi2': chi2_list_best})
                print(outputs[['M','O','H','R','chi2']])
            else:
                print('Insufficient models') 
                
        else:
            print('Please enter a best-matching filtering condition.')
    
    else:
            print('No output models')
        
    return
    
    

def gyre_p_spacing_import(gyre_work_dir,gyre_output_filename,gyre_filename):
    # Process individual GYRE model frequency outputs 
    
    with open(gyre_work_dir+"/"+gyre_output_filename+"/"+gyre_filename) as gyre_f:
        data = gyre_f.readlines()
        gyre_data = []
        for i,line in enumerate(data):
            if i >= 6:  # skip file header 
                data = line.split()                
                data = [float(d) for d in data]
                gyre_data.append(data)
      
    gyre_order = [g[0] for g in gyre_data]   
    gyre_freq = [g[1] for g in gyre_data]
    
    # convert uHz to cycles per day 
    gyre_freq = np.array(gyre_freq)*10**(-6)*86400
    # Periods in days
    gyre_pd = 1./np.array(gyre_freq)


    # Difference between periods of consecutive radial order n only
    gyre_p_spacing = []
    index_list = []
    for i in range(len(gyre_pd)-1):
        if abs(gyre_order[i+1] - gyre_order[i]) == 1:
            p_spacing = abs(gyre_pd[i+1] - gyre_pd[i]) 
            gyre_p_spacing.append(p_spacing)
            index_list.append(i)  # those to keep
            
    # use index list to make new period set which does not include those with skipped modes 
    gyre_period = [gyre_pd[i] for i in index_list] 
    
    return gyre_period,gyre_p_spacing  
    
    
    
def adipls_p_spacing_import(adipls_work_dir,adipls_filename):
    # Process individual ADIPLS model frequency outputs

    with open(adipls_work_dir+"/"+adipls_filename) as adipls_f:
        data = adipls_f.readlines()
        adipls_data = []
        for i,line in enumerate(data):
            data = line.split()                
            data = [float(d) for d in data]
            adipls_data.append(data)
            
    # only g-modes (n<0) and m=-1        
    adipls_order = [a[1] for a in adipls_data if a[1] < 0 and a[2] == m] 
    
    # only g-modes (n<0) and user-specified m; only freqs within defined freq range
    adipls_splitted_freq = [a[9] for a in adipls_data if a[1] < 0 and a[2] == m \
        and a[9] < nu2 and a[9] > nu1] 
    adipls_splitted_freq.sort()  
    
    # convert uHz to cycles per day 
    adipls_splitted_freq = np.array(adipls_splitted_freq)*10**(-6)*86400
    # Periods in days
    adipls_pd = 1./np.array(adipls_splitted_freq)
    

    # Difference between periods of consecutive radial order n only 
    adipls_p_spacing = []
    index_list = []
    for i in range(len(adipls_pd)-1):
        if abs(adipls_order[i+1] - adipls_order[i]) == 1:
            p_spacing = abs(adipls_pd[i+1] - adipls_pd[i]) 
            adipls_p_spacing.append(p_spacing)
            index_list.append(i)  
            
    adipls_period = [adipls_pd[i] for i in index_list]

    return adipls_period,adipls_p_spacing
        


def gyre_grid_scan(grid_work_dir,period_bench,p_spacing_bench,grad_bench,grad_err_bench, \
    intercept_bench,intercept_err_bench):
    # Loop over all GYRE output files inside grid directory 
    # main subroutine for controlling further computation

    # Create a single figure object for ploting all best-matching models
    fig1 = plt.figure(figsize=(14,8))
    ax1 = fig1.add_subplot(111)
    # Initialize param set for storing all grid models and the best-matching models 
    M_list = []
    O_list = []
    H_list = []
    R_list = []
    chi2_list = []
    M_list_best = []
    O_list_best = []
    H_list_best = []
    R_list_best = []
    chi2_list_best = []
    
    gyre_output_dir = fnmatch.filter(os.listdir(grid_work_dir), 'gyre_output_*')
    
    for gyre_output_filename in gyre_output_dir:
        # extract label from filename
        label = gyre_output_filename.replace('gyre_output_','')
        M,M_value, O,O_value, H,H_value, R,R_value = label.split('_')
        M_list.append(float(M_value))
        O_list.append(float(O_value))
        H_list.append(float(H_value))
        R_list.append(float(R_value))

        # go inside file to get the profile_puls_summary.txt
        gyre_puls_summary = os.listdir(grid_work_dir+"/"+gyre_output_filename)
        
        for gyre_filename in gyre_puls_summary:
            # call subroutine gyre_p_spacing_import
            period_grid,p_spacing_grid = gyre_p_spacing_import(grid_work_dir,gyre_output_filename,gyre_filename)
            # call subroutine gradient
            grad_grid,grad_err_grid,intercept_grid,intercept_err_grid = gradient(period_grid,p_spacing_grid)
            
            # Compute chi2 between this grid model and the benchmark - 
            # call subroutine data_interpolation 
            bench_P_output_list,grid_P_output_list = data_interpolation(period_grid,period_bench)
            # call subroutine compute_chi2
            chi2 = compute_chi2(bench_P_output_list,grid_P_output_list)
            chi2_list.append(chi2)
            
            # Filter best-matching models by gradient and intercept error of benchmark 
            grad_diff = abs(grad_grid - grad_bench)
            intercept_diff = abs(intercept_grid - intercept_bench)
            
            if bench_err_filter: 
                
                if grad_diff <= grad_err_bench and intercept_diff <= intercept_err_bench:
                    # call subroutine to plot p_spacings of best-matching models 
                    plot_p_spacings(ax1,period_grid,p_spacing_grid,period_bench,p_spacing_bench,chi2)                    
                    # store param set for the best-matching models 
                    M_list_best.append(float(M_value))
                    O_list_best.append(float(O_value))
                    H_list_best.append(float(H_value))
                    R_list_best.append(float(R_value))
                    chi2_list_best.append(chi2)
            
            # Temporarily changing filters to get more good-fit models       
            elif not bench_err_filter and grad_filter and intercept_filter:
                
                if grad_diff <= grad_filter and intercept_diff <= intercept_filter:
                    # call subroutine to plot p_spacings of best-matching models 
                    plot_p_spacings(ax1,period_grid,p_spacing_grid,period_bench,p_spacing_bench,chi2)                   
                    # store param set for the best-matching models 
                    M_list_best.append(float(M_value))
                    O_list_best.append(float(O_value))
                    H_list_best.append(float(H_value))
                    R_list_best.append(float(R_value))
                    chi2_list_best.append(chi2)    
                    
            else:
                print('Please enter a gradient filtering condition.')
                            
    return chi2_list,M_list,O_list,H_list,R_list,chi2_list_best,M_list_best,O_list_best, \
            H_list_best,R_list_best 



def adipls_grid_scan(grid_work_dir,period_bench,p_spacing_bench,grad_bench,grad_err_bench, \
    intercept_bench,intercept_err_bench):
    # Loop over all ADIPLS output files inside grid directory 
    # main subroutine for controlling further computation
    
    # Create figure object for ploting the best-matching models
    fig1 = plt.figure(figsize=(14,8))
    ax1 = fig1.add_subplot(111)
    # Initialize param set for all grid models and the best-matching models 
    M_list = []
    O_list = []
    H_list = []
    R_list = []
    chi2_list = []
    M_list_best = []
    O_list_best = []
    H_list_best = []
    R_list_best = []
    chi2_list_best = []
    
    adipls_output = fnmatch.filter(os.listdir(grid_work_dir), 'save_*_mode.data')
    
    for adipls_filename in adipls_output: 
        # extract label from filename
        label = adipls_filename.replace('save_','')
        label = label.replace('_mode.data','')   
        M,M_value, O,O_value, H,H_value, R,R_value = label.split('_') 
        M_list.append(float(M_value))
        O_list.append(float(O_value))
        H_list.append(float(H_value))
        R_list.append(float(R_value))       
    
        # call subroutine adipls_p_spacing_import
        period_grid,p_spacing_grid = adipls_p_spacing_import(grid_work_dir,adipls_filename)
        # call subroutine gradient 
        grad_grid,grad_err_grid,intercept_grid,intercept_err_grid = gradient(period_grid,p_spacing_grid)
        
        # Compute chi2 between this grid model and the benchmark 
        # call subroutine data_interpolation 
        P_output_list_1,P_output_list_2 = data_interpolation(period_grid,period_bench)
        # call subroutine compute_chi2
        chi2 = compute_chi2(P_output_list_1,P_output_list_2)
        chi2_list.append(chi2)  
        
        # Extract best-matching models with gradient filter 
        grad_diff = abs(grad_grid - grad_bench)
        intercept_diff = abs(intercept_grid - intercept_bench)
        
        if bench_err_filter: 
            
            if grad_diff <= grad_err_bench and intercept_diff <= intercept_err_bench:
                # call subroutine to plot p_spacings of best-matching models 
                plot_p_spacings(ax1,period_grid,p_spacing_grid,period_bench,p_spacing_bench,chi2)                    
                # store param set for the best-matching models 
                M_list_best.append(float(M_value))
                O_list_best.append(float(O_value))
                H_list_best.append(float(H_value))
                R_list_best.append(float(R_value))
                chi2_list_best.append(chi2)
        
        # Temporarily changing filters to get more good-fit models       
        elif not bench_err_filter and grad_filter and intercept_filter:
            
            if grad_diff <= grad_filter and intercept_diff <= intercept_filter:
                # call subroutine to plot p_spacings of best-matching models 
                plot_p_spacings(ax1,period_grid,p_spacing_grid,period_bench,p_spacing_bench,chi2)                   
                # store param set for the best-matching models 
                M_list_best.append(float(M_value))
                O_list_best.append(float(O_value))
                H_list_best.append(float(H_value))
                R_list_best.append(float(R_value))
                chi2_list_best.append(chi2)    
                
        else:
            print('Please enter a gradient filtering condition.')                   
           
    return chi2_list,M_list,O_list,H_list,R_list,chi2_list_best,M_list_best,O_list_best, \
            H_list_best,R_list_best 
          
                    

def gradient(period,p_spacing):   
    # Finds gradient (and error) of any given period spacing pattern 
    
    # Do not consider the short period range with large oscillations P < 1.5
    # and remove the corresponding p_spacing values
    period_cropped = []
    p_spacing_cropped = []
    
    for i,pd_set in enumerate(zip(period,p_spacing)):
        # proceed if period has values beyond 1.5 (some ADIPLS don't)
        if max(period) > 1.5:
            if pd_set[0] >= 1.5:
                period_cropped.append(pd_set[0])
                p_spacing_cropped.append(pd_set[1])
        else:        
            if pd_set[0] >= 1.0: # otherwise set filter at 1.0
                period_cropped.append(pd_set[0])
                p_spacing_cropped.append(pd_set[1])        
    
    # Defining cut-off with gradient of benchmark 
    poly_coeff,covar_matrix = np.polyfit(period_cropped,p_spacing_cropped,1,cov=True)
    # diagonal of covar matrix are variance for each coefficient    
    grad = poly_coeff[0]
    grad_err = np.sqrt(covar_matrix[0][0])
    intercept = poly_coeff[1]
    intercept_err = np.sqrt(covar_matrix[1][1])
    
    return grad,grad_err,intercept,intercept_err

                    
                                                            
def data_interpolation(period_grid,period_bench):
    # Interpolation - Make same number of points in period_grid and period_bench
    # before calculating chi2
    
    bench_P_output_list = []
    grid_P_output_list = []

    # Use the period list with greater number of points to define the cut boundaries
    if len(period_grid) >= len(period_bench):
        period_min = period_bench[-1]
        period_max = period_bench[0]
    elif len(period_bench) > len(period_grid):
        period_min = period_grid[-1]
        period_max = period_grid[0]
    
    # Create a set of specific_P reference points, then assign a closest period from the 
    # grid and one from the bench to each point
    for specific_P in np.arange(period_min,period_max,0.02):  
            
        # Interpolating ADIPLS results  
        bench_P_output_list.append(min(period_bench, key=lambda x:abs(x-specific_P)))
        
        # Interpolating GYRE results 
        grid_P_output_list.append(min(period_grid, key=lambda x:abs(x-specific_P)))
        
    return bench_P_output_list,grid_P_output_list



def compute_chi2(P_output_list_1,P_output_list_2):
   # Calculate chi2 between pulsation freqs from P_output_list_1 and P_output_list_2
   
    i_max = len(P_output_list_1)   
    chi2_sum = 0
    for i in range(0,i_max): 
        chi2_sum += (((1/P_output_list_1[i]) - (1/P_output_list_2[i]))/0.00068)**2    
        # 0.00068 - Rayleigh limit of Kepler; frequency resolution
    chi2_sum = chi2_sum*(1/(i_max-4.))     
              
    return chi2_sum
                                                        
         
                                                                                                               
def plot_p_spacings(ax1,period_grid,p_spacing_grid,period_bench,p_spacing_bench,chi2):
    
    opacity = 1./chi2 * plot_scale_factor
    
    ax1.plot(period_grid,p_spacing_grid,'-',color='seagreen',markersize=1,label='Grid',alpha=opacity)
    ax1.plot(period_bench,p_spacing_bench,'-',color='red',markersize=1,label='Benchmark')
    
    ax1.set_xlabel('Period (days)',fontsize=27,labelpad=10)
    ax1.set_ylabel('Period Spacing (days)',fontsize=27,labelpad=22) 
    ax1.set_xbound(0.0,3.3)
    ax1.set_ybound(0.015,0.045)
    ax1.tick_params(labelsize=17)
    plt.grid(True)
    plt.show()
    
    return                                                                     
                                                                                    


def param_correlation(chi2_list,M_list,O_list,H_list,R_list,label_bench):
    # Loop over each combination pair in the parameter set w the while keeping the oher two 
    # parameters constant at around the benchmark param
    # Plot the 2D variation of chi2 values with each pair of parameters 
    
    # extract benchmark param_set 
    M_b,M_value_b, O_b,O_value_b, H_b,H_value_b, R_b,R_value_b = label_bench.split('_') 
    bench_param_val_set = [float(M_value_b),float(O_value_b),float(H_value_b),float(R_value_b)]
    bench_param_set = ['M','O','H','R']
    
    # call subroutine gen_combination to generate all grid parameter combinations (and labels)
    param_full_list = M_list,O_list,H_list,R_list
    param_full_label = 'M','O','H','R'
    
    comb_list = gen_combination(*param_full_list)  # * to unpack tuple
    comb_label = gen_combination(*param_full_label)
    
    # Loop over each combination in the param_set
    for param_comb_set in zip(comb_list,comb_label):  
        
        param_comb_list = param_comb_set[0]
        param_comb_label = param_comb_set[1]
        
        # the two varied parameters 
        vary_param_list1 = param_comb_list[0]   # e.g. M_list
        vary_param_list2 = param_comb_list[1]   # e.g. O_list
        vary_param1 = param_comb_label[0]   # e.g. M
        vary_param2 = param_comb_label[1]   # e.g. O
                
        # the other two contant parameters
        alt_comb_set = []
        for param in zip(param_full_list,param_full_label):
            if param[0] != vary_param_list1 and param[0] != vary_param_list2:
                alt_comb_set.append(param)
                
        alt_param_comb_list = [alt_set[0] for alt_set in alt_comb_set] # e.g. [H_list,R_list]
        alt_param_comb_label = [alt_set[1] for alt_set in alt_comb_set]
                        
        const_param_list1 = alt_param_comb_list[0]  # e.g. H_list
        const_param_list2 = alt_param_comb_list[1]  # e.g. R_list
        const_param1 = alt_param_comb_label[0]   # e.g. H
        const_param2 = alt_param_comb_label[1]   # e.g. R


        # Within each const_param_list, obtain the value closest to corresponding param of benchmark 
        for param_set in zip(bench_param_val_set,bench_param_set):
            if param_set[1] == const_param1:   # At e.g. both are 'H' strings
                const_param1_closest = min(const_param_list1, key=lambda x:abs(x-param_set[0]))
            if param_set[1] == const_param2:
                const_param2_closest = min(const_param_list2, key=lambda x:abs(x-param_set[0]))
        
        const_param_index_list = []
        extracted_vary_param_list1 = []
        extracted_vary_param_list2 = []
        extracted_chi2_list = []
        
        # find the index list of models of which both const_param are kept constant at values close to benchmark
        for i,(const_param1_val,const_param2_val) in enumerate(zip(const_param_list1,const_param_list2)):
            if const_param1_val == const_param1_closest and const_param2_val == const_param2_closest:
                const_param_index_list.append(i)
                # using the model index, extract the required varying parameter list along with chi2 
                extracted_vary_param_list1.append(vary_param_list1[i])
                extracted_vary_param_list2.append(vary_param_list2[i])
                extracted_chi2_list.append(chi2_list[i])
                
        
        # call subroutine plot_contour to create 2D surface plot of chi2 for each param_comb
        plot_contour(vary_param1,vary_param2,extracted_vary_param_list1,extracted_vary_param_list2,extracted_chi2_list)
        
    return 
        
        
        
def gen_combination(list1,list2,list3,list4):
    # Generate all combinations of param_list; order not important
    
    param_full_list1 = [list1,list2,list3,list4]
    param_full_list2 = [list1,list2,list3,list4]
    comb_list = []
    for param1 in param_full_list1:
        for param2 in param_full_list2:
            if param1 != param2:
                if ([param1,param2]) not in comb_list and ([param2,param1]) not in comb_list:
                    comb_list.append([param1,param2])  
                                     
    return comb_list   # in matrix form ([param1,param2],[param1,param2],...)
        

            
def plot_contour(param1_label,param2_label,param_list1,param_list2,chi2_list):
    # subroutine for 2D contour plot of a given pair of param1 and param2
    
    # Sort param_list values in ascending order, sort the corresponding chi2 values 
    list_set = zip(param_list1,param_list2,chi2_list)
    l = sorted(list_set)
    param_list1,param_list2,chi2_list = zip(*l)    
    # NB: order of param_list1 and param_list2 do not affect final results 
    
    # param_list1 - number of columns, x-axis
    # param_list2 - number of rows, y-axis
    # chi2_matrix shape - (number of rows, number of columns)
 
    compact_param_list2 = [] # without duplicated values 
    chi2_matrix = []
    
    for i,param_2 in enumerate(param_list2): # loop through each y value        
        if param_2 not in compact_param_list2:  
            compact_param_list2.append(param_2)
                
            # define a new row list at each looped y-value
            row_list = []
            chi2_row = []
            for param1,param2,chi2 in zip(param_list1,param_list2,chi2_list):
                # at where param2 is equal to the param_2 (y value) being processed at this iteration
                if param2 == param_2:                    
                    row_list.append(param1)
                    # add in the chi2 value to each x on this row
                    chi2_row.append(chi2)
                    
        # do not put in row for duplicated values in param_list2
        if chi2_row not in chi2_matrix:        
            chi2_matrix.append(chi2_row)

    compact_param_list1 = []
    for param_1 in param_list1:
        if param_1 not in compact_param_list1:
            compact_param_list1.append(param_1)
    
    # Plot contour with the two axis and 2D matrix 
    
    label = ['M','O','H','R']
    label_name = ['Stellar Mass','Convective core overshoot','Hydrogen abundance','Rotation rate']
    for l,l_name in zip(label,label_name):
        if l == param1_label:
            param1_label_n = l_name
        elif l == param2_label:
            param2_label_n = l_name
    
    fig2 = plt.figure(figsize=(10.5,9))
    ax2 = fig2.add_subplot(111)
    levels = np.arange(min(chi2_list),max(chi2_list),(max(chi2_list)-min(chi2_list))/50.)
    cmap = plt.get_cmap('gnuplot2')
    CS = ax2.contour(compact_param_list1,compact_param_list2,chi2_matrix,levels=levels,cmap=cmap)
    plt.colorbar(CS)
#    ax2.clabel(CS, fontsize=10, inline=1,color='black',fmt='%6.0f')
    ax2.set_xlabel(param1_label_n,fontsize=27,labelpad=10)
    ax2.set_ylabel(param2_label_n,fontsize=27,labelpad=10)
    ax2.tick_params(labelsize=13)
    xstart,xend = ax2.get_xlim()
    ystart,yend = ax2.get_ylim()
    xstepsize = ((xend-xstart)/5.)
    ystepsize = ((yend-ystart)/5.)
    plt.xticks(np.arange(xstart,xend+xstepsize,step=xstepsize))
    plt.yticks(np.arange(ystart,yend+ystepsize,step=ystepsize))
    plt.show()
    
    return
            
            
            
main(benchmark_work_dir,grid_work_dir)
                            
                
                                











 