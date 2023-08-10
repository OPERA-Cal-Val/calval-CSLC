import papermill as pm
import warnings
import os
warnings.filterwarnings('ignore')

# Start runtime evaluation
import timeit
start = timeit.default_timer()

f = open('./Rosamond/t064_135523_iw2_datelist_all.txt')   #reading text file with a single column of dates
datels = f.read().splitlines()

#Parameters
data_dir = 's3://opera-provisional-products/CSLC/pst_adt_common/gamma_v.0.3/Rosamond/Ascending'
save_dir = './Rosamond/A064_run4'
burst_id = 't064_135523_iw2'
snr_threshold = 15
solidtide = 'True'
cr_network = 'Rosamond'

# Create folders
os.makedirs(f'{save_dir}/pngs', exist_ok=True)
os.makedirs(f'{save_dir}/ipynbs', exist_ok=True)
os.makedirs(f'{save_dir}/summary', exist_ok=True)
os.makedirs(f'{save_dir}/crdata', exist_ok=True)

#Check if file exist
outcsv1 = 'ALE_{cr_network}_{burst_id}_allDates.csv'
outcsv2 = 'ALE_{cr_network}_allCRs.csv'

if os.path.isfile(outcsv1):
   os.remove(outcsv1)

if os.path.isfile(outcsv2):
   os.remove(outcsv2)

# Run the ALE for each date
for d in datels:
    print(f'Processing date: {d}')
    pm.execute_notebook('ALE_COMPASS_Stream.ipynb',
                        f'{save_dir}/ipynbs/ALE_COMPASS_{burst_id}_{d}.ipynb',
                        parameters={'data_dir':data_dir,
                                    'save_dir':save_dir,
                                    'burst_id':burst_id,
                                    'date':d, 
                                    'snr_threshold':snr_threshold, 
                                    'solidtide':solidtide,
                                    'cr_network':cr_network},
                        kernel_name='calval-CSLC')
    
# End runtime evaluation
stop = timeit.default_timer()
print(f'Time: ', (stop - start)/60, 'min.')