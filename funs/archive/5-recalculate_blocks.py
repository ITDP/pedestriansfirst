def recalculate_blocks(input_folder_prefix = 'cities_out/',current_year=2024,block_patch_length=1000):
    """re-calculate block density 'patches' from existing block geodata
    """
    
    folders = os.listdir(input_folder_prefix)
    for folder in tqdm(folders):
        print(folder)
        blocks_latlon = gpd.read_file(f"{input_folder_prefix}/{folder}/geodata/blocks/blocks_latlon_{current_year}.geojson")
        blocks_utm = ox.project_gdf(blocks_latlon)
        
        block_patches_latlon, block_unbuffered_patches_latlon = pedestriansfirst.make_patches(
            gpd.GeoDataFrame(geometry=[shapely.geometry.box(*blocks_latlon.total_bounds)], crs=4326), 
            blocks_utm.crs, 
            patch_length=block_patch_length
            )
        centroids = blocks_latlon.centroid
        block_unbuf_patches_utm = block_unbuffered_patches_latlon.to_crs(blocks_utm.crs)
        patch_densities = block_unbuffered_patches_latlon
        for patch_idx  in list(patch_densities.index):
            try:
                patch_densities.loc[patch_idx,'block_count'] = centroids.intersects(patch_densities.loc[patch_idx,'geometry']).value_counts()[True]
                #import pdb; pdb.set_trace()
            except KeyError:
                patch_densities.loc[patch_idx,'block_count'] = 0 
        patch_densities_utm = patch_densities.to_crs(blocks_utm.crs)
        patch_densities_utm['density'] = patch_densities_utm.block_count / (patch_densities_utm.area / 1000000)
        patch_densities_latlon = patch_densities_utm.to_crs(epsg=4326)
        try:
            os.rename(f"{input_folder_prefix}/{folder}/geodata/blocks/block_densities_latlon_{current_year}.geojson",
                  f"{input_folder_prefix}/{folder}/geodata/blocks/OLD_block_densities_latlon_{current_year}.geojson")
        except FileNotFoundError:
            pass
        patch_densities_latlon.to_file(f"{input_folder_prefix}/{folder}/geodata/blocks/block_densities_latlon_{current_year}.geojson", driver='GeoJSON')

def make_block_only_folder(input_folder_prefix = 'cities_out/', output_folder_prefix='cities_out_block_data_only/'):
    if not os.path.exists(output_folder_prefix):
        os.makedirs(output_folder_prefix)
    
    folders = os.listdir(input_folder_prefix)
    for folder in tqdm(folders):
        os.makedirs(f'{output_folder_prefix}/{folder}/geodata/blocks')
        shutil.copy(f"{input_folder_prefix}/{folder}/geodata/blocks/block_densities_latlon_2024.geojson",f"{output_folder_prefix}/{folder}/geodata/blocks/block_densities_latlon_2024.geojson")
        shutil.copy(f"{input_folder_prefix}/{folder}/geodata/blocks/blocks_latlon_2024.geojson",f"{output_folder_prefix}/{folder}/geodata/blocks/blocks_latlon_2024.geojson")
