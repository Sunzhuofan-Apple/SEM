import cellxgene_census

with cellxgene_census.open_soma() as census:
    
    cell_metadata = census["census_data"]["homo_sapiens"].obs.read(
        value_filter = "sex == 'female' and cell_type in ['microglial cell', 'neuron']",
        column_names = ["assay", "cell_type", "tissue", "tissue_general", 
"suspension_type", "disease"]
    )

    cell_metadata = cell_metadata.concat()

    cell_metadata = cell_metadata.to_pandas()

    print(cell_metadata)

