##################################################################
# Plot Period spacing pattern 1st local Minimum against Rotation
##################################################################

# For models only varying in rotation rate, other parameters kept constant
# No other model outputs in work dir

import numpy as np
import matplotlib.pyplot as plt 
import os, fnmatch

# Set work directories
adipls_work_dir = 'work_adipls_11' 
gyre_work_dir = 'work_gyre_11'

# Set azimuthal order m
m = 0  

# Frequency range (uHz) for filtering adipls output
nu1 = 3.8
nu2 = 38.6



def main_program(adipls_work_dir,gyre_work_dir):
    
    adipls_rotation_rates, adipls_first_min_list = adipls_file_scan(adipls_work_dir)
    gyre_rotation_rates, gyre_first_min_list = gyre_file_scan(gyre_work_dir)
    
    # plot results on same graph 
    fig = plt.figure(figsize=(13,9))
    ax = fig.add_subplot(111)
    ax.plot(adipls_rotation_rates,adipls_first_min_list,'o',color='blue',markersize=8,label='2nd order')
    ax.plot(gyre_rotation_rates,gyre_first_min_list,'o',color='orange',markersize=8,label='TAR')
    
#    ax.set_title('Location of 1st dip against Rotation',fontsize=20)
    ax.set_xlabel('Rotation rate (critical)',fontsize=27,labelpad=22)
    ax.set_ylabel('Period of 1st dip (days)',fontsize=27,labelpad=22) 
    ax.set_xbound(0.0,0.17)
    ax.set_ybound(0.45,0.57)
    ax.tick_params(labelsize=17)
    ax.legend(prop={'size': 20})
    plt.grid(True)
    plt.show()
    
    return 
    
    

def adipls_file_scan(adipls_work_dir):
    # loop over all adipls output files 
    adipls_rotation_rates = []
    adipls_first_min_list = []
    
    adipls_output = fnmatch.filter(os.listdir(adipls_work_dir), 'save_*_mode.data')

    for adipls_filename in adipls_output: 
        # extract label from filename
        adipls_label = adipls_filename.replace('save_','')
        adipls_label = adipls_label.replace('_mode.data','')
        # store rotation value
        M,M_value, O,O_value, H,H_value, R,R_value = adipls_label.split('_')
        adipls_rotation_rates.append(float(R_value))  
        
        # call subroutine adipls_p_list
        adipls_period,adipls_p_spacing = adipls_p_list(adipls_work_dir,adipls_filename)
        # call subroutine crop_data
        adipls_period_cropped,adipls_p_spacing_cropped = crop_data(adipls_period,adipls_p_spacing)
        # call subroutine first_min
        adipls_first_min = first_min(adipls_period_cropped,adipls_p_spacing_cropped)
        adipls_first_min_list.append(adipls_first_min)
    
    return adipls_rotation_rates, adipls_first_min_list



def gyre_file_scan(gyre_work_dir):
    # loop over all gyre output files 
    gyre_rotation_rates = []
    gyre_first_min_list = []

    gyre_output_dir = fnmatch.filter(os.listdir(gyre_work_dir), 'gyre_output_*')
    
    for gyre_output_filename in gyre_output_dir:
        # extract label from filename
        gyre_label = gyre_output_filename.replace('gyre_output_','')
        # store rotation value
        M,M_value, O,O_value, H,H_value, R,R_value = gyre_label.split('_')
        gyre_rotation_rates.append(float(R_value)) 

        # go inside file to get the profile_puls_summary.txt
        gyre_puls_summary = os.listdir(gyre_work_dir+"/"+gyre_output_filename)
        
        for gyre_filename in gyre_puls_summary:
            # call subroutine gyre_p_list
            gyre_period,gyre_p_spacing = gyre_p_list(gyre_work_dir,gyre_output_filename,gyre_filename)
            # call subroutine crop_data
            gyre_period_cropped,gyre_p_spacing_cropped = crop_data(gyre_period,gyre_p_spacing)
            # call subroutine first_min
            gyre_first_min = first_min(gyre_period_cropped,gyre_p_spacing_cropped)
            gyre_first_min_list.append(gyre_first_min)
                         
    return gyre_rotation_rates, gyre_first_min_list



def adipls_p_list(adipls_work_dir,adipls_filename):
    # process individual model freq outputs 
    adipls_data = []
    with open(adipls_work_dir+"/"+adipls_filename) as adipls_f:
        data = adipls_f.readlines()
        for i in range(len(data)):
            line = data[i].split()
            adipls_data.append(line)

    # only freqs within defined freq range
    adipls_splitted_freq = [float(a[9]) for a in adipls_data if float(a[1]) < 0 \
        and float(a[2]) == m and float(a[9]) < nu2 and float(a[9]) > nu1]         
    adipls_splitted_freq.sort()  # ascending order            
    # convert uHz to cycles per day 
    adipls_splitted_freq = np.array(adipls_splitted_freq)*10**(-6)*86400
    # Periods in days
    adipls_period = 1./np.array(adipls_splitted_freq)

    # Difference between consecutive periods 
    adipls_p_spacing = [abs(adipls_period[i+1] - adipls_period[i]) for i in range(len(adipls_period)-1)]        
    adipls_p_spacing.reverse()

    # remove last element in period - same dimension as p_spacing
    adipls_period = adipls_period[:-1]             
    adipls_period = adipls_period.tolist()  
    adipls_period.reverse()
    
    return adipls_period,adipls_p_spacing 
            
            

def gyre_p_list(gyre_work_dir,gyre_output_filename,gyre_filename):
    # process individual model freq outputs 
    gyre_data = []
    with open(gyre_work_dir+"/"+gyre_output_filename+"/"+gyre_filename) as gyre_f:
        data = gyre_f.readlines()
        for i in range(len(data)):
            if i >= 6:
                line = data[i].split()
                gyre_data.append(line)
        
    gyre_freq = [float(g[1]) for g in gyre_data]
    # convert uHz to cycles per day 
    gyre_freq = np.array(gyre_freq)*10**(-6)*86400
    # Periods in days
    gyre_period = 1./np.array(gyre_freq)

    # Difference between consecutive periods 
    gyre_p_spacing = [abs(gyre_period[i+1] - gyre_period[i]) for i in range(len(gyre_period)-1)]

    # remove last element in period - same dimension as p_spacing
    gyre_period = gyre_period[:-1]    
    
    return gyre_period,gyre_p_spacing    



def crop_data(period,p_spacing):
    # Remove data with period >0.6 days 
    period_cropped = [pd for pd in period if pd <= 0.6]
    
    # Remove the corresponding values in p_spacing 
    # find index of those remained in period_cropped
    p_spacing_cropped = []
    for pd_c in period_cropped:
        for i,j in enumerate(period):
            if j == pd_c:
                p_spacing_cropped.append(p_spacing[i])  
    return period_cropped,p_spacing_cropped 
   
   
       
def first_min(period,p_spacing):
    current_p_spacing = p_spacing[0] # start from first one
    turning_p_spacing = p_spacing[0] # initialize variable
    turning_p_spacing_list = []
    
    for p_spacing_element in p_spacing:
        if p_spacing_element <= current_p_spacing:
            current_p_spacing = p_spacing_element # store into current_p_spacing in each iteration
        elif p_spacing_element > current_p_spacing: #passed turning point
            turning_p_spacing = p_spacing_element
            turning_p_spacing_list.append(turning_p_spacing)
            
    # get index of this turning point
    for i,j in enumerate(p_spacing):
        if j == turning_p_spacing_list[0]: 
            turning_index = i-1  # real min point at index -1
    return period[turning_index]                  



main_program(adipls_work_dir,gyre_work_dir)






