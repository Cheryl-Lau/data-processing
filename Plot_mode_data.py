import numpy as np 
import matplotlib.pyplot as plt 

gyre_work_dir = 'work_gyre_R0'
gyre_output_filename = 'gyre_output_M_1.55_O_0.03_H_0.5_R_0'
gyre_filename = 'profile_puls_summary.txt'

adipls_work_dir = 'work_adipls_R0'
adipls_filename = 'save_M_1.55_O_0.03_H_0.5_R_0_mode.data'


# Settings for filtering adipls output -
# Enter azimuthal order m
m = 0  
# Frequency range (uHz) 
nu1 = 3.8
nu2 = 38.6



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
            index_list.append(i)  # those to keep
            
    # use index list to make new period set which does not include those with skipped modes 
    adipls_period = [adipls_pd[i] for i in index_list]

    return adipls_period,adipls_p_spacing
    
    
gyre_period,gyre_p_spacing = gyre_p_spacing_import(gyre_work_dir,gyre_output_filename,gyre_filename)
adipls_period,adipls_p_spacing = adipls_p_spacing_import(adipls_work_dir,adipls_filename)

fig1 = plt.figure(figsize=(13,9.2))
ax1 = fig1.add_subplot(111)
ax1.plot(adipls_period,adipls_p_spacing,'-',color='blue',markersize=1,label='2nd Pert.')
ax1.plot(gyre_period,gyre_p_spacing,'-',color='orange',markersize=2,label='TAR')


ax1.set_xlabel('Period (days)',fontsize=27,labelpad=22)
ax1.set_ylabel('Period Spacing (days)',fontsize=27,labelpad=22) 
ax1.set_xbound(0.0,3.3)
ax1.set_ybound(0.0,0.06)
ax1.tick_params(labelsize=17)
plt.grid(True)
plt.show()