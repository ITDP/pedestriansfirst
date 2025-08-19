if __name__ == '__main__':
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
    ox.utils.config(log_console = False)
    ucdb = gpd.read_file('input_data/ghsl/SMOD_V1s6_opr_P2023_v1_2020_labelUC_DB_release.gpkg')
    ucdb.index =  ucdb['ID_UC_G0']
    hdcs_to_test = [11480 , 461 , 576 , 8265 , 4494] #ITDP "Cycling Cities" below 500k
    hdcs_to_test += list(ucdb[(int(sys.argv[2]) < ucdb['POP_2020'])&(ucdb['POP_2020'] < int(sys.argv[1]))].sort_values('POP_2020', ascending=False).ID_UC_G0)
    for hdc in hdcs_to_test:
        hdc = int(hdc)
        #if len(sys.argv) == 1:
        #divide_by = 1
        #remainder = 0
        # else: 
        divide_by = int(sys.argv[3])
        remainder = int(sys.argv[4])
        print (f"{hdc}%{divide_by}={hdc % divide_by}, compare to {remainder}, {ucdb.loc[hdc,'NAME_MAIN']}, pop {ucdb.loc[hdc,'POP_2020']}")
        if hdc % divide_by == remainder and ucdb.loc[hdc,'NAME_MAIN'] != 'N/A':
            if not os.path.exists(f'cities_out/ghsl_region_{hdc}/indicator_values.csv'):
                if os.path.exists(f'cities_out/ghsl_region_{hdc}/geodata/blocks/blocks_latlon_2024.geojson'):
                    regional_analysis(hdc)#, analyze=False)
                else:
                    regional_analysis(hdc)
    calculate_country_indicators()
