

import dask.delayed
import random
from dask import compute, delayed
from dask.distributed import Client
pd.options.display.float_format = '{:.2f}'.format


# hdcs = pd.read_csv('input_data/hdc_to_run.csv', dtype={'hdc_new' : str})
hdcs = pd.read_csv('input_data/pbf/pbf_hdc_country.csv', dtype={'hdc' : str})
hdcs = hdcs.sort_values('pop', ascending=True)

hdcs_all = hdcs["hdc"].tolist()







#########################################################################################
# RUN STEP BY STEP
#########################################################################################

# first, run only juris and prep (not parallel, since it already runs parallel)

 for hdc in hdcs_all:
  regional_analysis(hdc, 
                    to_test =  [],
                    analyse=False, 
                    summarize=False, 
                    jurisdictions=True, 
                    prep = False,
                    current_year=2023)
                    
                    
# now, parallel
# @dask.delayed
def run_analysis(hdc):
  return regional_analysis(
    hdc,
    to_test=['pnpb', #protected bikeways
            'pnab', #all bikeways
            # 'healthcare',
                # 'schools',
                # 'hs',
                # 'bikeshare',
                # 'carfree',
                'blocks',
                # 'density',
                'pnft', # VOLTAR PRO PADRAO!!!
                'pnrt',

                'pnst', #combo transit + bike
                'highways'],
    summarize=False,
    jurisdictions=False,
    current_year=2024,
    prep=False,
    analyse=True
    )

# set up cluster with 20 workers. Each worker uses 1 thread and has a 64GB memory limit.
client = Client(n_workers=10,threads_per_worker=1,memory_limit='64GB')

# have a look at your workers
client


# # Initialize the Dask client to enable parallel execution
# client = Client()

# from 1 to 700: 20-22 workers
# from 700 to 900: 14-16 workers
# from 900 to 1000: 8-10 workers
# from 900 to 1114: 4-6 workers

@dask.delayed
def run_analysis_safe(hdc):
    try:
        return run_analysis(hdc)  # Your original function
    except Exception as e:
        return f"Error in {hdc}: {e}"  # Return error message instead of crashing

# Create a list of delayed tasks
results = [run_analysis_safe(hdc) for hdc in hdcs_all[501:900]]
# Compute the results (this will trigger the parallel execution)
results1 = compute(*results)

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

client = Client(n_workers=3,threads_per_worker=1,memory_limit='48GB')

# Create a list of delayed tasks
results = [run_analysis_safe(hdc) for hdc in hdcs_left[0:50]]
# Compute the results (this will trigger the parallel execution)
results1 = compute(*results)


 for hdc in hdcs_left:
   run_analysis(hdc)


#########################################################################################
# RUN ONLY PNFT!
#########################################################################################v

# List all files in the directory
files = os.listdir("input_data/gtfs/2024")

# Extract numeric patterns from filenames
hdc_gtfs = {re.search(r"(?<=gtfs_)\d{5}", f).group() for f in files if re.search(r"(?<=gtfs_)\d{5}", f)}

# Convert the set to a sorted list (optional)
hdc_gtfs = sorted(hdc_gtfs)
len(hdc_gtfs)


# now, parallel
# @dask.delayed
def run_analysis(hdc):
  return regional_analysis(
    hdc,
    to_test=['pnab', 'pnpb'],
    summarize=False,
    jurisdictions=False,
    current_year=2024,
    prep=False,
    analyse=True
    )

# set up cluster with 20 workers. Each worker uses 1 thread and has a 64GB memory limit.
client = Client(n_workers=20, threads_per_worker=1,memory_limit='64GB')

# have a look at your workers
client



@dask.delayed
def run_analysis_safe(hdc):
    try:
        return run_analysis(hdc)  # Your original function
    except Exception as e:
        return f"Error in {hdc}: {e}"  # Return error message instead of crashing
      
# filter only the ones with GTFS

# where is 08154
hdcs_all.index('08154')

# Create a list of delayed tasks
results = [run_analysis_safe(hdc) for hdc in hdcs_all[968:969]]
# Compute the results (this will trigger the parallel execution)
results1 = compute(*results)

out = hdcs_all[1004:1117]
# remove the ones without the gtfs
out1 =  list(set(out) - set(hdc_gtfs))
len(out1)
# Create a list of delayed tasks
results = [run_analysis_safe(hdc) for hdc in out1]
# Compute the results (this will trigger the parallel execution)
results1 = compute(*results)

  #########################################################################################
# calculate the indicators! (summarize = True)
#########################################################################################

def run_indicators(hdc):
  return regional_analysis(
    hdc,
    to_test=[],
    analyse=False,
    summarize=True,
    jurisdictions=False,
    current_year=2024,
    prep=False
  )

@dask.delayed
def run_indicators_safe(hdc):
    try:
        return run_indicators(hdc)  # Your original function
    except Exception as e:
        return f"Error in {hdc}: {e}"  # Return error message instead of crashing

# Close the client after computation is done
client.close()

client = Client(n_workers=20,
                threads_per_worker=1,
                memory_limit='48GB')


# Create a list of delayed tasks
results = [run_indicators_safe(hdc) for hdc in hdcs_all[1001:1117]]
# Compute the results (this will trigger the parallel execution)
results1 = compute(*results)


Marcstart = timeit.default_timer()
run_analysis("00178")
end = timeit.default_timer()



#########################################################################################
# country calculations
#########################################################################################


calculate_country_indicators(current_year = 2023)
calculate_country_indicators(current_year = 2024)








#########################################################################################
# new calculations - may 2025
#########################################################################################



# first, run only juris and prep (not parallel, since it already runs parallel)
def run_analysis(hdc):
  return regional_analysis(
    hdc,
    to_test=[],
    summarize=False,
    jurisdictions=True,
    current_year=2024,
    prep=False,
    analyse=False
    )

client = Client(n_workers=4,threads_per_worker=1,memory_limit='64GB')



@dask.delayed
def run_analysis_safe(hdc):
    try:
        return run_analysis(hdc)  # Your original function
    except Exception as e:
        return f"Error in {hdc}: {e}"  # Return error message instead of crashing

# Create a list of delayed tasks
results = [run_analysis_safe(hdc) for hdc in hdcs_all]
# Compute the results (this will trigger the parallel execution)
results1 = compute(*results)



# only pnrt
# now, parallel
# @dask.delayed
def run_analysis(hdc):
  return regional_analysis(
    hdc,
    to_test=['pnrt'],
    summarize=False,
    jurisdictions=False,
    current_year=2024,
    prep=False,
    analyse=True
    )

# set up cluster with 20 workers. Each worker uses 1 thread and has a 64GB memory limit.
client = Client(n_workers=6,threads_per_worker=1,memory_limit='64GB')

# have a look at your workers
client



@dask.delayed
def run_analysis_safe(hdc):
    try:
        return run_analysis(hdc)  # Your original function
    except Exception as e:
        return f"Error in {hdc}: {e}"  # Return error message instead of crashing

# Create a list of delayed tasks
results = [run_analysis_safe(hdc) for hdc in hdcs_left]
# Compute the results (this will trigger the parallel execution)
results1 = compute(*results)

http://localhost:8787/


# Close the client after computation is done
client.close()               






# run anaylises for everythiong
def run_analysis(hdc):
  return regional_analysis(
    hdc,
    to_test=['pnpb', #protected bikeways
            'pnab', #all bikeways
            'healthcare',
                'schools',
                'hs',
                'bikeshare',
                'carfree',
                'blocks',
                'density',
                'pnft', # VOLTAR PRO PADRAO!!!
                'pnrt',
                'pnst', #combo transit + bike
                'highways'],
    summarize=True,
    jurisdictions=False,
    prep=False,
    analyse=False,
    current_year=2024
    )

# set up cluster with 20 workers. Each worker uses 1 thread and has a 64GB memory limit.
client = Client(n_workers=23,threads_per_worker=1,memory_limit='64GB')


@dask.delayed
def run_analysis_safe(hdc):
    try:
        return run_analysis(hdc)  # Your original function
    except Exception as e:
        return f"Error in {hdc}: {e}"  # Return error message instead of crashing

# Create a list of delayed tasks
results = [run_analysis_safe(hdc) for hdc in hdcs_all]
# Compute the results (this will trigger the parallel execution)
results1 = compute(*results)

http://localhost:8787/


# Close the client after computation is done
client.close()               



#########################################################################################
# country calculations
#########################################################################################


calculate_country_indicators(current_year = 2023, rt_and_pop_years = [2023])
calculate_country_indicators(current_year = 2024, rt_and_pop_years = [1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2024, 2025])
