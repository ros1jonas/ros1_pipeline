import GEOparse
import pandas
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import numpy as np
from math import log10

import geo_retrieve

if __name__ == '__main__':
    soft_file_string = geo_retrieve.get_geo_file("GSE64027", "getmap")
    gse = geo_retrieve.get_gse_obj(soft_file_string)
    meta_pd = geo_retrieve.get_metadata_pandas(gse, "alternate_metadata")
    geo_retrieve.get_metadata_pandas(gse)
    rawcounts_pd = (geo_retrieve.get_counts_pandas_from_raw_counts
                    ("handmatige_counts_download/GSE64027_kassie05NNK.txt", meta_pd))
    designstring, control_name, test_list = geo_retrieve.get_designstring_controlname_testlist(meta_pd)
    dds = geo_retrieve.get_deseq2dataset(meta_pd, rawcounts_pd, designstring)
    lfc = geo_retrieve.getfoldchange(dds, "foldchange")
    multidict = geo_retrieve.get_multidict(dds, designstring, control_name,
                                     test_list)
    dict_met_volcs = geo_retrieve.get_multi_volcano(multidict, 20, 5, True)