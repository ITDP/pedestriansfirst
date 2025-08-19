def get_number_jurisdictions(input_folder_prefix = 'cities_out/'):
    """return the total number of jursidictions in all calculated agglomerations
    """
    folders = os.listdir(input_folder_prefix)
    total = 0
    for folder in tqdm(folders):
        data = pd.read_csv(f'{input_folder_prefix}/{folder}/indicator_values.csv')
        if 'admin_level' in data.columns:
            total += len(data[~np.isnan(data.admin_level)])
    print(total)
    return total
