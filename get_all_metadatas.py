import GEOparse
import pandas
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import numpy as np
from math import log10

import geo_retrieve

def get_all_metadatas(gse_file: str, output_file: str):
    """
    This function creates a list of all possible test designs from the
    GSEs provided in the file "lung_cancer_data_files.xlsx" and saves
    this in a new file. All metadatas will also be downloaded to the
    folder "all_metadatas".

    :param gse_file: location of file "lung_cancer_data_files.xlsx"
    :param output_file: name of output file.
    """
    all_lung_cancer_data = pandas.read_excel(gse_file)
    alldesignspandas = pandas.DataFrame(columns=['GSE', 'gse_title', 'Characteristic', 'Design', 'Number', ])
    for index, row in all_lung_cancer_data.iterrows():
        soft_file_string = geo_retrieve.get_geo_file(row["GEO ID"].strip(), "getmap")
        gse = geo_retrieve.get_gse_obj(soft_file_string)
        meta_pd = geo_retrieve.get_metadata_pandas(gse, "all_metadatas/" + row["GEO ID"].strip() + "metadata")
        design_dict = geo_retrieve.get_design_dict(meta_pd)
        for x in design_dict:
            for y in design_dict[x]:
                alldesignspandas = pandas.concat([alldesignspandas,
                                                  pandas.DataFrame({'GSE': [row["GEO ID"].strip()],
                                                                    'gse_title': [row["Title"]],
                                                                    'Characteristic': [x],
                                                                    'Design': [y],
                                                                    'Number': [design_dict[x][y]]})],
                                                 ignore_index=True)
    alldesignspandas.to_excel(output_file)



if __name__ == '__main__':
    get_all_metadatas("lung_cancer_data_files.xlsx", "all_ros1_GSE_designs.xlsx")
