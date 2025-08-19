    
def download_ghsl(proj='mw', resolution='1000'):
    """Download population data from the Global Human Settlement Layer
    
    Can get data in mollewiede or WGS84 projections
    Hard-coded to use 2023 release
    """
    
    proj_code = {'mw':54009,'ll':4326}[proj]
    if not os.path.exists('input_data/'):
        os.mkdir('input_data/')
    if not os.path.exists(f'input_data/ghsl_data_{resolution}m_{proj}/'):
        os.mkdir('input_data/ghsl_data_{resolution}m_mw/')
    for year in tqdm(range(1975, 2031, 5)): 
        letter='E'
        if proj == 'mw':
            name = f'GHS_POP_{letter}{year}_GLOBE_R2023A_{proj_code}_{resolution}'
        if not os.path.exists(f'input_data/ghsl_data_{resolution}m_{proj}/{name}/{name}_V1_0.tif'):
            zippath, _ = urllib.request.urlretrieve(f'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/{name}/V1-0/{name}_V1_0.zip')
            with zipfile.ZipFile(zippath, "r") as f:
                f.extractall(f'input_data/ghsl_data_{resolution}m_{proj}/{name}')
