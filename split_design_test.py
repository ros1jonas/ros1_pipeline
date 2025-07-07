import GEOparse
import pandas
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import numpy as np
from math import log10

import geo_retrieve

# GSE255958


if __name__ == '__main__':
    soft_file_string = geo_retrieve.get_geo_file("GSE255958", "getmap")
    gse = geo_retrieve.get_gse_obj(soft_file_string)
    meta_pd = geo_retrieve.get_metadata_pandas(gse, "pre_split_metadata")
    geo_retrieve.get_metadata_pandas(gse)
    design_dict = geo_retrieve.get_design_dict(meta_pd)
    # print(design_dict)
    # geo_retrieve.get_combined_design_meta(meta_pd, ["characteristics_ch1.4.treatment",
    #                                                 "characteristics_ch1.3.treatment"],
    #                                       "characteristics_combined", "combined_pre_split_meta")
    combine_meta, new_collumnname = geo_retrieve.get_combined_design_meta(meta_pd, ["characteristics_ch1.4.treatment",
                                                    "characteristics_ch1.3.treatment"],
                                          "characteristics_combined")
    # geo_retrieve.get_split_design_meta(combine_meta, "characteristics_combined")
    combined_meta, simplified_col_name = geo_retrieve.create_simplified_design(combine_meta, "characteristics_combined",
                                          [["untreated"],["Osimertinib", "9d"]], excelstring="simplified_meta")
    # geo_retrieve.get_count_pd_from_tar("handmatige_counts_download/GSE255958_RAW.tar", combined_meta)
    organoid_meta_pd = geo_retrieve.get_meta_pd_from_xlsx("handmatig_aangepaste_meta/organoid_255958_meta.xlsx")
    organoid_counts_pd = geo_retrieve.get_count_pd_from_tar("handmatige_counts_download/GSE255958_RAW.tar", organoid_meta_pd)
    print(organoid_counts_pd)
    # designstring, control_name, test_list = geo_retrieve.get_designstring_controlname_testlist(organoid_meta_pd)
    dds = geo_retrieve.get_deseq2dataset(organoid_meta_pd, organoid_counts_pd, "characteristics_ch1.3.treatment")
    dds.deseq2()
    lfc = geo_retrieve.getfoldchange(dds, "foldchange")
    multidict = geo_retrieve.get_multidict(dds, "characteristics_ch1.3.treatment", "Control",
                                     ["Osimertinib"])
    dict_met_volcs = geo_retrieve.get_multi_volcano(multidict, 20, 5, True)