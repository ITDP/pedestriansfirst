

import dask.delayed
import random
from dask import compute, delayed
from dask.distributed import Client
pd.options.display.float_format = '{:.2f}'.format

@dask.delayed
def run_analysis(hdc):
  return regional_analysis(
    hdc,
    to_test=['pnrt'],
    analyse=True,
    summarize=False,
    jurisdictions=False,
    current_year=2023,
    prep=False,
  )

# hdcs = pd.read_csv('input_data/hdc_to_run.csv', dtype={'hdc_new' : str})
hdcs = pd.read_csv('input_data/pbf/pbf_hdc_country.csv', dtype={'hdc' : str})
hdcs = hdcs.sort_values('pop', ascending=True)

hdcs_all = hdcs["hdc"].tolist()



# set up cluster with 20 workers. Each worker uses 1 thread and has a 64GB memory limit.
client = Client(n_workers=6,
                threads_per_worker=1,
                memory_limit='24GB')

# have a look at your workers
client


# # Initialize the Dask client to enable parallel execution
# client = Client()

# from 1 to 700: 20-22 workers
# from 700 to 900: 14-16 workers
# from 900 to 1000: 8-10 workers
# from 900 to 1114: 4-6 workers

start = timeit.default_timer()
# Create a list of delayed tasks
results = [run_analysis(hdc) for hdc in hdcs_all[951:1114]]
# Compute the results (this will trigger the parallel execution)
compute(*results)
end = timeit.default_timer()

http://localhost:8787/


# Close the client after computation is done
client.close()


#########################################################################################
# SEE THE ONES WE ARE MISSING!
#########################################################################################

# Get list of files in "logs" directory
hdcs_ran = os.listdir("logs/")

# Filter files that contain the word "running"
hdcs_ran = [file for file in hdcs_ran if "finished" in file]

# Extract the first 5-digit number from each file name
hdcs_ran = [re.search(r'\d{5}', file).group() for file in hdcs_ran if re.search(r'\d{5}', file)]

# Get the difference between all IDs and the running ones (hdcs_left)
hdcs_left = list(set(hdcs_all) - set(hdcs_ran))
len(hdcs_left)

# Close the client after computation is done
client.close()

client = Client(n_workers=8,
                threads_per_worker=1,
                memory_limit='24GB')

# Create a list of delayed tasks
results = [run_analysis(hdc) for hdc in hdcs_left]
# Compute the results (this will trigger the parallel execution)
compute(*results)



#########################################################################################
# RUN ADDITIONAL CITIES
#########################################################################################


# 03105 - guadalajara, run with buffer = 40000
 hdcs_additional = ['02006' , '02249' , '03048', '99999'] #ITDP "Cycling Cities" below 500k
 # 02249 already run
 
 for hdc in hdcs_additional:
  regional_analysis("02006", 
                    to_test =  [
                         'healthcare',
                         'schools',
                         'hs',
                         'bikeshare',
                         'carfree',
                         'blocks',
                         'density',
                         'pnft', # VOLTAR PRO PADRAO!!!
                         'pnrt',
                         'pnpb', #protected bikeways
                         'pnab', #all bikeways
                         'pnst', #combo transit + bike
                         'highways',
                         ],
                    analyse=True, 
                    summarize=True, 
                    jurisdictions=False, 
                    prep = False,
                    current_year=2023)
                    
                    
 regional_analysis("03105", 
                    to_test =  [
                         'healthcare',
                         'schools',
                         'hs',
                         'bikeshare',
                         'carfree',
                         'blocks',
                         'density',
                         'pnft', # VOLTAR PRO PADRAO!!!
                         'pnrt',
                         'pnpb', #protected bikeways
                         'pnab', #all bikeways
                         'pnst', #combo transit + bike
                         'highways',
                         ],
                    analyse=True, 
                    summarize=True, 
                    jurisdictions=False, 
                    prep = False,
                    current_year=2023)


#########################################################################################
# calculate the indicators! (summarize = True)
#########################################################################################

@dask.delayed
def run_calculations(hdc):
  return regional_analysis(
    hdc,
    to_test=[],
    analyse=False,
    summarize=True,
    jurisdictions=False,
    current_year=2023,
    prep=False
  )


# Close the client after computation is done
client.close()

client = Client(n_workers=24,
                threads_per_worker=1,
                memory_limit='24GB')


# Create a list of delayed tasks
results = [run_calculations(hdc) for hdc in hdcs_all]
# Compute the results (this will trigger the parallel execution)
compute(*results)


start = timeit.default_timer()
run_analysis("00178")
end = timeit.default_timer()



#########################################################################################
# Run PNFT for all




#########################################################################################

# test for jakarta
regional_analysis(
    "07910",
    to_test=["pnft"],
    analyse=True,
    summarize=True,
    jurisdictions=False,
    current_year=2023,
    prep=False
  )
  
  
  regional_analysis(
    "01156",
    to_test=["pnft"],
    analyse=True,
    summarize=True,
    jurisdictions=False,
    current_year=2023,
    prep=False
  )
