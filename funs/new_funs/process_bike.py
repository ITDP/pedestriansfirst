 def process_bike(G_allroads, quilt_allbike, quilt_protectedbike, center_nodes):
 
 # convert graph to data.frame
        ways_gdf = ox.graph_to_gdfs(G_allroads, nodes=False)
        # ways_gdf.to_csv('teste_kaue/Output.csv', index = False)
        print("Finished building gdfs graph")
        
        # remove columns?
        ways_gdf = ways_gdf.drop(columns=['name', 'width', 'oneway', 'reversed', 
                                          'length', 'lanes', 'ref', 'maxspeed', 
                                          'access', 'bridge', 'service', 'junction', 'area'])
        
        for col in ['highway','cycleway','bicycle','cycleway:left','cycleway:right','cycleway:both']:
            if not col in ways_gdf.columns:
                ways_gdf[col] = ''
                
        tagged_cycleways = ways_gdf[(ways_gdf['highway'] == 'cycleway')]
        cycle_paths = ways_gdf[(ways_gdf['highway'] == 'path') & (ways_gdf['bicycle'] == 'designated')]
        on_street_tracks = ways_gdf[(ways_gdf['cycleway:right'] == 'track') |
                                       (ways_gdf['cycleway:left'] == 'track') |
                                       (ways_gdf['cycleway:both'] == 'track') |
                                       (ways_gdf['cycleway'] == 'track') |
                                       (ways_gdf['cycleway:right'] == 'opposite_track') |
                                        (ways_gdf['cycleway:left'] == 'opposite_track') |
                                        (ways_gdf['cycleway:both'] == 'opposite_track') |
                                        (ways_gdf['cycleway'] == 'opposite_track')
                                       ]
        on_street_lanes = ways_gdf[(ways_gdf['cycleway:right'] == 'lane') |
                                       (ways_gdf['cycleway:left'] == 'lane') |
                                       (ways_gdf['cycleway:both'] == 'lane') |
                                       (ways_gdf['cycleway'] == 'lane')]
        
        total_protectedbike = pd.concat([tagged_cycleways, cycle_paths, on_street_tracks])
        total_allbike = pd.concat([total_protectedbike, on_street_lanes])
        
        print("Finished filtering")
        
                
        del ways_gdf
        gc.collect()
        
        print("Finished cleaning memory")
        
        # exclude tiny, unconnected segments that aren't near larger ones
        max_jump = 300
        min_radius = 1000
        
        #remove isolated small lanes from protected network
        total_protectedbike['in_real_network'] = "unknown"
        for idx in total_protectedbike.index:
            already_identified = total_protectedbike.loc[idx,'in_real_network'] == "unknown"
            if type(already_identified).__name__ == 'Series':
                already_identified = already_identified.any()
            if already_identified:
                connected_indices = [idx]
                if shapely.minimum_bounding_radius(total_protectedbike.loc[connected_indices,'geometry'].unary_union) > min_radius:
                    total_protectedbike.loc[connected_indices,'in_real_network'] = "yes"
                else:
                    for i in range(0,10): #just so we don't end up in an infinite loop somehow
                        connected_network = total_protectedbike.loc[connected_indices,'geometry'].unary_union
                        nearby = total_protectedbike[total_protectedbike.distance(connected_network) < max_jump]
                        if 'yes' in nearby.in_real_network.unique():
                            total_protectedbike.loc[connected_indices,'in_real_network'] = "yes"
                            break
                        if set(connected_indices) == set(nearby.index):
                            if shapely.minimum_bounding_radius(total_protectedbike.loc[connected_indices,'geometry'].unary_union) > min_radius:
                                total_protectedbike.loc[connected_indices,'in_real_network'] = "yes"
                            else:
                                total_protectedbike.loc[connected_indices,'in_real_network'] = "no"
                            break
                        else:
                            connected_indices = list(nearby.index)
                    if shapely.minimum_bounding_radius(total_protectedbike.loc[connected_indices,'geometry'].unary_union) > min_radius:
                        total_protectedbike.loc[connected_indices,'in_real_network'] = "yes"
        total_protectedbike = total_protectedbike[total_protectedbike.in_real_network == "yes"]
        
        print("Finish remove isolated small lanes from protected network")
        
        #remove isolated small lanes from allbike network
        total_allbike['in_real_network'] = "unknown"
        for idx in total_allbike.index:
            already_identified = total_allbike.loc[idx,'in_real_network'] == "unknown"
            if type(already_identified).__name__ == 'Series':
                already_identified = already_identified.any()
            if already_identified:
                connected_indices = [idx]
                if shapely.minimum_bounding_radius(total_allbike.loc[connected_indices,'geometry'].unary_union) > min_radius:
                    total_allbike.loc[connected_indices,'in_real_network'] = "yes"
                else:
                    for i in range(0,10): #just so we don't end up in an infinite loop somehow
                        connected_network = total_allbike.loc[connected_indices,'geometry'].unary_union
                        nearby = total_allbike[total_allbike.distance(connected_network) < max_jump]
                        if 'yes' in nearby.in_real_network.unique():
                            total_allbike.loc[connected_indices,'in_real_network'] = "yes"
                            break
                        if set(connected_indices) == set(nearby.index):
                            if shapely.minimum_bounding_radius(total_allbike.loc[connected_indices,'geometry'].unary_union) > min_radius:
                                total_allbike.loc[connected_indices,'in_real_network'] = "yes"
                            else:
                                total_allbike.loc[connected_indices,'in_real_network'] = "no"
                            break
                        else:
                            connected_indices = list(nearby.index)
                    if shapely.minimum_bounding_radius(total_allbike.loc[connected_indices,'geometry'].unary_union) > min_radius:
                        total_allbike.loc[connected_indices,'in_real_network'] = "yes"
        total_allbike = total_allbike[total_allbike.in_real_network == "yes"]
        
        print("Finish remove isolated small lanes from allbike network")
                
        quilt_allbike = pd.concat([quilt_allbike, total_allbike])
        quilt_protectedbike = pd.concat([quilt_protectedbike, total_protectedbike])
        
        center_nodes['pnpb'] = set()
        center_nodes['pnab'] = set()
        for edge in total_protectedbike.index:
            center_nodes['pnpb'].add(edge[0])
            center_nodes['pnpb'].add(edge[1])
        for edge in total_allbike.index:
            center_nodes['pnab'].add(edge[0])
            center_nodes['pnab'].add(edge[1])

            
        return quilt_allbike, quilt_protectedbike, center_nodes
