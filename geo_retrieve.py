import GEOparse
import pandas
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import numpy as np
from math import log10
import json
import tarfile
import os
import gzip

import Volcano_obj


def get_count_pd_from_tar(tar_path : str, meta_pd : pandas):
    """
    Creates a pandas object containing the counts data from a tar file
    and a pandas metadata object. If certain samples are omitted from
    the metadata, then they will also be omitted in the counts pandas
    object. If the "supplementary_file_1" column in the metadata is
    not fully filled out, the counts pandas object will not be complete.

    :param str tar_path: File path to tar file.
    :param pandas meta_pd: pandas dataframe containing the metadata.
    :return: dataframe containing the raw counts.
    :rtype: pandas
    """
    data = {}
    gene_index = []
    index_done = False
    with tarfile.open(tar_path, "r") as tf:
        for x in tf.getmembers():
            tf.extract(member=x.name)
            with gzip.open(x.name, 'rt') as f:
                file_content = f.read()
            meta_title = ""
            for index, row in meta_pd.iterrows():
               if x.name in row["supplementary_file_1"]:
                   meta_title = index
            if meta_title != "":
                templist = []
                counts_loc = None
                for y in file_content.split("\n"):
                    if counts_loc == None:
                        teller = 0
                        for z in y.split("\t"):
                            if z in ["expected_count"]:
                                counts_loc = teller
                            teller += 1
                    elif y != "":
                        templist.append(round(float(y.split("\t")[counts_loc])))
                        if index_done == False:
                            gene_index.append((y.split("\t")[0]))
                index_done = True
                data[meta_title] = templist
            os.remove(x.name)
    counts_df = pandas.DataFrame(data, index=gene_index)
    return counts_df.T


def get_meta_pd_from_xlsx(metastring : str):
    """
    Create a metadata pandas dataframe from a metadata Excel file.

    :param str metastring: File path to metadata xlsx file.
    :return: dataframe containing the raw counts.
    :rtype: pandas
    """
    meta_pd = pandas.read_excel(metastring)
    return meta_pd.set_index("title")


def get_geo_file(geo_code: str, dest_map: str=""):
    """
    Downloads a soft file of an GSE number. The soft file is either saved
    to a user specified folder (if the parameter dest_map is filled), or
    directly saved in the current directory
    (if the parameter dest_map is not filled). This function will always
    return a string containing the full file path to the created
    soft file.

    A soft file can contain data tables and accompanying descriptive
    information of GSE records.
    This soft file will be further used to create a GSE object.

    :param str geo_code: GEO database identifier of GSE file.
    :param str dest_map: Name of folder where soft file is downloaded to.
    :return: file path to downloaded soft file.
    :rtype: int
    """
    if dest_map != "":
        GEOparse.get_GEO(geo=geo_code, destdir="./" + dest_map + "/")
        return "./" + dest_map + "/"+geo_code+"_family.soft.gz"
    else:
        GEOparse.get_GEO(geo=geo_code)
        return geo_code+"_family.soft.gz"
    # GEOparse.get_GEO_file(geo=geo_code, destdir="./getmap/")
    # gse = GEOparse.get_GEO(geo="GSE1563", destdir="./")


def get_gse_obj(soft_file_string: str)->GEOparse.GSE:
    """
    Creates an GSE object with a soft file. This function receives
    a single parameter, soft_file_string, a string which
    contains the file path to the soft file.
    This function returns the GSE object.

    The GSE object will be used to create a metadata table.
    the metadata is necessary to create a dds (DeseqDataSet) object, which will
    be further used to create ds (DeseqStats) objects.
    The metadata is also necessary to create the raw counts pandas table
    which is also required to create the dds object.

    :param str soft_file_string: File path to soft file to retrieve.
    :return: GSE object
    :rtype: gse
    """
    gse = GEOparse.get_GEO(filepath =soft_file_string)
    return (gse)


def get_metadata_pandas(gse, excelstring=None)->pandas.DataFrame:
    """
    Creates a pandas table containing the sequencing metadata from a GSE object.
    Optionialy, the user may opt to save the metadata in an Excel file
    by filling out the excelstring parameter. The metadata is returned
    as a pandas table.

    The metadata contains, among others, the various test designs.
    This metadata is necessary to create a dds (DeseqDataSet) object, which will
    be further used to create ds (DeseqStats) objects.
    The metadata is also necessary to create the raw counts pandas table
    which is also required to create the dds object.

    :param GEOparse.GSE gse: GSE object from which metadata is retrieved.
    :param str excelstring: Name of file where metadata will be saved,
        nothing is saved if empty.
    :return: Table containing the sequencing metadata.
    :rtype: pandas
    """
    meta_pd = gse.phenotype_data
    indexed_pd = meta_pd.set_index("title")
    if excelstring != None:
        indexed_pd.to_excel(excelstring + ".xlsx")
    return indexed_pd


def get_counts_pandas_from_raw_counts(rawcounts, meta_pd, excelstring=None)->pandas.DataFrame:
    """
    Creates a pandas table containing all the gene counts
    per sample from a raw counts file. This function also
    requires the metadata pandas table to get a list of sample names.
    The counts table is arranged in such a way that it can be
    utilized by DeseqDataSet without further modification.
    The user may opt to save the metadata in an Excel file by
    filling out the excelstring parameter.

    The counts table is necessary to create a dds (DeseqDataSet) object,
    which will be further used to create ds (DeseqStats) objects.

    :param str rawcounts: File path to a .txt file containing the reads.
    :param pandas meta_pd: Table containing the sequencing metadata.
    :param str excelstring: Name of file where count data will be saved,
        nothing is saved if empty.
    :return: Table containing all the gene counts per sample.
    :rtype: pandas
    """
    rawcounts_pd = pandas.read_table(rawcounts, index_col=0)
    testlist = meta_pd.index.to_list()
    rawcounts_dropped = rawcounts_pd[testlist].copy()
    rawcounts_dropped.reindex_like(rawcounts_dropped)
    rawcounts_pd_transposed = rawcounts_dropped.T
    if excelstring != None:
        rawcounts_dropped.to_excel(excelstring + ".xlsx")
    return rawcounts_pd_transposed


def filterrawcounts(rawcounts_pd):
    """
    Unused
    """
    print(rawcounts_pd)
    genes_to_keep = rawcounts_pd.columns[rawcounts_pd.sum(axis=0) >= 10]
    trimmed_counts = rawcounts_pd[genes_to_keep]
    print (trimmed_counts)
    # trimmed_counts.to_excel("trimmed_counts.xlsx")
    return trimmed_counts


def get_deseq2dataset(meta_pd, rawcounts_pd, designstring)->DeseqDataSet:
    """
    Creates a Deseq2 dataset object containing all the test cases.
    As parameters, this function requires a pandas table containing
    the metadata, a pandas table containing the counts data and the
    designstring, a string which contains the name of the column
    containing the test design names.

    The dds (DeseqDataSet) object is used to create ds (DeseqStats)
    objects. ds objects contain, among others, tables with the P-values
    and fold changes of the different test cases, these will be
    used to fill multidicts which are used to create volcano plots.

    :param pandas meta_pd: Table containing the sequencing metadata
    :param pandas rawcounts_pd: Table containing all the gene counts per sample.
    :param str designstring: Name of column containing different test cases.
    :return: Deseq2 dataset object containing all the test cases.
    :rtype: DeseqDataSet
    """
    dds = DeseqDataSet(
    counts=rawcounts_pd,
    metadata=meta_pd,
    design='~'+ designstring)
    dds.deseq2()
    return dds


def get_deseqstats(dds, designstring, tested_lvl: str, control_name: str)->DeseqStats:
    """
    Creates a Deseq2 stats object containing p-value estimations
    and fold changes of a specific test case. This function
    requires as parameters: a dds object, a string which
    contains the name of the column containing the test design
    names (designstring), the name of the test design (tested_lvl) and
    the name of the control test (control_name).

    Ds objects contain, among others, tables with the P-values
    and fold changes of the different test cases, these will be
    used to fill multidicts which are used to create volcano plots.

    :param pydeseq2.DeseqStats dds: Deseq2 dataset object containing all the test cases.
    :param str designstring: Name of column containing the test design names.
    :param str control_name: Name of the control test case.
    :param str tested_lvl: Name of the test case of which P-values will be estimated.
    :return: Deseq2 stats object containing p-value estimations for a specific test case.
    :rtype: DeseqStats
    """
    ds = DeseqStats(
        dds,
        contrast=[designstring, tested_lvl, control_name],
        alpha=0.05,
        cooks_filter=True,
        independent_filter=True,
    )
    ds.run_wald_test()
    return ds


def make_p_value_dict(dds, control_name: str, tested_list: list):
    """
    unused
    """
    ds_list = {}
    for tested_lvl in tested_list:
        stats = get_deseqstats(dds, tested_lvl, control_name)
        ds_list[tested_lvl] = stats.p_values
        stats.p_values.to_excel("p_valuestest_" + tested_lvl + ".xlsx")
        normaltest = pandas.DataFrame(stats.base_mean)
        normaltest.to_excel("basemeans" + tested_lvl + ".xlsx")
    return ds_list


def get_multidict(dds, designstring, control_name: str, tested_list: list, excelstring=None)->dict:
    """
    Creates a dictionary containing the P values,
    fold change and base means per test case.
    A dds (deseqdataset) object is required as the first parameter,
    This dds is used to generate ds (deseqstats) objects for each
    test case, from which the data is retrieved.



    The user may opt to save the multidict in an Excel file by
    filling out the excelstring parameter. The Excel file will contain
    a sheet for each test design, with each sheet containing a column
    for the P values, fold changes and base means.



    :param pydeseq2.DeseqStats dds: Deseq2 dataset object containing all the test cases.
    :param str designstring: Name of column containing the test design names.
    :param str control_name: Name of the control test case.
    :param tested_list: List with names of all test cases.
    :param excelstring: Name of file where test case data will be saved,
        nothing is saved if empty.
    :return: Dictionary containing the P values,
        fold change and base means per test case.
    :rtype: Dict
    """
    multidict = {}
    for tested_lvl in tested_list:
        stats = get_deseqstats(dds, designstring, tested_lvl, control_name)
        for label, content in stats.LFC.items():
            if label not in ["Ensembl_ID", "Intercept"]:
                if label.split("[")[1] == "T." + tested_lvl + "]":
                    multidict[tested_lvl] = {
                        "p_values": stats.p_values,
                        "fold_change": content,
                        "base_means": pandas.Series(stats.base_mean)
                    }
    if excelstring != None:
        writer = pandas.ExcelWriter(excelstring+ '.xlsx', engine='xlsxwriter')
        names = list(multidict.keys())
        dataframes = []
        for x in multidict:
            temppandas = multidict[x]["p_values"].to_frame(name="P values")
            temppandas["fold change"] = multidict[x]["fold_change"]
            temppandas["base means"] = multidict[x]["base_means"].to_list()
            dataframes.append(temppandas)
        for i, frame in enumerate(dataframes):
            frame.to_excel(writer, sheet_name=names[i])
        writer.close()
    return multidict


def get_multi_volcano(multidict, p_limit: int, fold_limit: int, save=False)->dict:
    """
    Loops through a dictionary containing test cases with
    their P values and log folds. Adds a volcano
    plot object to each test case.

    :param dict multidict: Dictionary containing the P values,
        fold change and base means per test case.
    :param int p_limit: Maximum limit to which higher p-values wil be reduced.
    :param int fold_limit: Maximum limit to which higher and
        lower fold values wil be reduced.
    :param save: If true, every volcano plot will be saved to pdf file.
    :return: Dictionary containing the P values,
        fold change, base means and volcano plot per test case.
    :rtype: Dict
    """
    saveteller = 1
    for condition_name in multidict:
        if save == True:
            multidict[condition_name]["volcano"] = Volcano_obj.Volcano_obj(multidict[condition_name], condition_name,
                               p_limit, fold_limit, "volcano"+str(saveteller)+".pdf")
            saveteller+=1
        else:
            multidict[condition_name]["volcano"] = Volcano_obj.Volcano_obj(multidict[condition_name],
                                                                           condition_name, p_limit, fold_limit)
    return multidict


def get_single_volcano(conditiondict, condition_name, p_limit: int, fold_limit: int, save_string=None):
    """
    Unused
    """
    modified_p_values_list = []
    for x in conditiondict["p_values"]:
        if (log10(x) * -1) > p_limit:
            modified_p_values_list.append(p_limit)
        else:
            modified_p_values_list.append(log10(x) * -1)
    modified_p_values = pandas.Series(modified_p_values_list)
    numpy_fold = conditiondict["fold_change"].to_numpy()
    modified_fold_list = []
    for x in numpy_fold:
        if x > fold_limit:
            modified_fold_list.append(fold_limit)
        elif x < -fold_limit:
            modified_fold_list.append(-fold_limit)
        else:
            modified_fold_list.append(x)
    modified_fold = pandas.Series(modified_fold_list)
    plt.scatter(modified_fold, modified_p_values, marker='.')
    plt.xlabel('FoldChange(log2)')
    plt.ylabel('P-value(-10log)')
    plt.title(condition_name)
    if save_string != None:
        plt.savefig(save_string)
        plt.clf()


def getfoldchange(dds, excelstring=None):
    """
    Unused

    Create a pandas table containing the gene expression fold changes
    per condition compared to the control group.

    :param pydeseq2.dds dds: Deseq2 dataset.
    :return: Table containing the gene expression fold changes.
    :rtype: pandas
    """
    dds.fit_LFC()
    lfc = dds.varm["LFC"]
    if excelstring != None:
        lfc.to_excel(excelstring + ".xlsx")
    return lfc



def volcano_plot(p_values, lfc, p_limit: int, fold_limit: int):
    """
    Unused

    Creates volcano plots based on provided p_values and foldchanges.
    
    :param pandas p_values: 
    :param pandas lfc:
    :param int p_limit: Maximum limit to which higher p-values wil be reduced.
    :param int fold_limit: Maximum limit to which higher and lower p-values wil be reduced.
    """
    numpy_p_values = p_values.to_numpy()
    modified_p_values_list = []
    for x in numpy_p_values:
        if (log10(x)*-1) > p_limit:
            modified_p_values_list.append(p_limit)
        else:
            modified_p_values_list.append(log10(x)*-1)
    modified_p_values = pandas.Series(modified_p_values_list)
    # colors = ['tab:blue', 'tab:orange', 'tab:green']
    color_count = 0
    for label, content in lfc.items():
        if label != "Intercept":
            numpy_fold = content.to_numpy()
            modified_fold_list = []
            for x in numpy_fold:
                if x > fold_limit:
                    modified_fold_list.append(fold_limit)
                elif x < -fold_limit:
                    modified_fold_list.append(-fold_limit)
                else:
                    modified_fold_list.append(x)
            modified_fold = pandas.Series(modified_fold_list)
            plt.scatter(modified_fold, modified_p_values, marker='.')
            color_count+=1
            # plt.xlim(-5, 5)
            # plt.ylim(0, 20)
            plt.xlabel('FoldChange(log2)')
            plt.ylabel('P-value(-10log)')
            plt.savefig( "volcano" + str(color_count) + ".pdf")
            plt.clf()


def get_designstring_controlname_testlist(meta_pd)->tuple[str, str, list]:
    """
    Retrieves from the metadata the name of the column containing
    the names of test designs, the name of the control test case and
    a list containing all the test case names.

    :param pandas meta_pd: Table containing the sequencing metadata.
    :return: Name of the column containing the names of test designs.
    :rtype: str
    :return: The name of the control test case.
    :rtype: str
    :return: List containing all the test case names.
    :rtype: list
    """
    designstring = ""
    control_name = ""
    test_list = []
    for (columnName, columnData) in meta_pd.items():
        if "_" in columnName:
            if columnName.split("_")[0] == "characteristics":
                if columnName.split(".")[2] == "treatment":
                    designstring = columnName
                    for testcase in columnData:
                        if control_name == "":
                            control_name = testcase
                        if control_name != testcase:
                            if testcase not in test_list:
                                test_list.append(testcase)
    if designstring == "":
        raise Exception("designstring couldn't be found")
    else:
        return designstring, control_name, test_list


def get_design_dict(meta_pd : pandas):
    """
    Creates from a metadata pandas dataframe a dictionary containing
    all characteristics columns containing more than one design.
    The values for each column contain further dictionaries containing
    the different designs as keys and the number of times these designs
    appear as values.

    These dicts can be used to create an overview of different possible
    designs of multiple GSE series.

    :param meta_pd: Pandas dataframe containing the metadata.
    :return: dictionary containing
        all characteristics columns containing more than one design.
    :rtype: dict
    """
    design_dict = {}
    for (columnName, columnData) in meta_pd.items():
        if columnName.split("_")[0] == "characteristics":
            counter = 0
            testdict = {}
            for index, value in columnData.items():
                if value not in testdict.keys():
                    testdict[value] = 1
                else:
                    testdict[value] += 1
                counter+=1
            if len(testdict) > 1:
                design_dict[columnName.split("_")[1]] = testdict
    return design_dict


def get_combined_design_meta(meta_pd : pandas, designlist : list, combined_col_name : str = None,  excelstring=None):
    """
    This function receives a metadata pandas dataframe and a list
    containing two metadata columns to be combined into a new
    third column, which is added to the metadata dataframe and
    returned. by filling out the parameter "combined_col_name" a new
    custom name can be given to this combined collumn, else the function
    will make one by combining the names of the provided collumns.
    The new metadata is also saved as a xlsx file if "excelstring"
    is filled out.


    :param meta_pd: metadata pandas dataframe
    :param designlist: list containing two metadata columns to be
        combined into a new third column
    :param combined_col_name: Name of new column, optional.
    :param excelstring: Name of file where new metadata will be saved,
        nothing is saved if empty.
    :return: New metadata pandas dataframe with added combined column.
    :rtype: pandas dataframe
    :return: Name of new column
    :rtype: str
    """

    new_collumnname = ""
    for x in designlist:
        new_collumnname += x
    if combined_col_name is not None:
        new_collumnname = combined_col_name
    serieslist = []
    for index, row in meta_pd.iterrows():
        combined_content = ""
        for x in designlist:
            if combined_content != "":
                combined_content += " "
            if str(row.loc[x]) != "nan":
                combined_content += str(row.loc[x])
        serieslist.append(row.add(pandas.Series([combined_content, index], index=[new_collumnname, "title"]), fill_value=""))
    combined_meta = pandas.concat(serieslist, axis = 1).T.set_index("title")
    if excelstring != None:
        combined_meta.to_excel(excelstring + ".xlsx")
    return combined_meta, new_collumnname


def create_simplified_design(
        meta_pd : pandas,
        targetdesign : str,
        simplify_list : list,
        simplified_col_name : str = "char_simple",
        excelstring=None):
    """
    Adds a new collumn to a metadata pandas dataframe containing
    the labels "target", "control" or "other", depending on certain
    keywords provided by the parameter simplify_list and whether these
    appear in the column "targetdesign".

    simplify_list contains two elements in following order:
    elements that signify control (list with strings)
    elements that signify target (list with strings)

    currently ALL elements have to be present


    :param meta_pd: pandas dataframe containing metadata.
    :param targetdesign: Column from which simplified column
        will be created
    :param simplify_list: list containg two list containing keywords
        that define "control" and "target".
    :param simplified_col_name: Name of new column, optional.
    :param excelstring: Name of file where new metadata will be saved,
        nothing is saved if empty.
    :return: New metadata pandas dataframe with added simple column.
    :rtype: pandas dataframe
    :return: Name of new column.
    :rtype: str
    """
    serieslist = []
    for index, row in meta_pd.iterrows():
        new_content = "other"
        for x in simplify_list[0]:
            if x in str(row.loc[targetdesign]):
                new_content = "control"
            else:
                new_content = "other"
                break
        if new_content != "control":
            for x in simplify_list[1]:
                if x in str(row.loc[targetdesign]):
                    new_content = "target"
                else:
                    new_content = "other"
                    break
        serieslist.append(
            row.add(pandas.Series([new_content, index], index=[simplified_col_name, "title"]), fill_value=""))
    combined_meta = pandas.concat(serieslist, axis=1).T.set_index("title")
    if excelstring != None:
        combined_meta.to_excel(excelstring + ".xlsx")
    return combined_meta, simplified_col_name

# def get_split_design_meta(
#         meta_pd : pandas,
#         targetdesign : str,
#         split_intervals : list,
#         split_on : list = [" "],
#         split_column_names : list = []
# ):
#     """
#     Will probably be unused.
#
#     :param meta_pd:
#     :param targetdesign:
#     :param split_intervals:
#     :param split_on:
#     :param split_column_names:
#     :return:
#     """
#     pass


if __name__ == '__main__':
    soft_file_string = get_geo_file("GSE229195", "getmap")
    gse = get_gse_obj(soft_file_string)
    meta_pd = get_metadata_pandas(gse)
    # get_metadata_pandas(gse, "metadata")
    rawcounts_pd = (get_counts_pandas_from_raw_counts
                    ("handmatige_counts_download/GSE229195_Raw_gene_counts_matrix.txt", meta_pd, "rawcountstest"))
    # filterrawcounts(rawcounts_pd)
    designstring, control_name, test_list = get_designstring_controlname_testlist(meta_pd)
    dds = get_deseq2dataset(meta_pd, rawcounts_pd, designstring)
    # lfc = getfoldchange(dds, "foldchange")
    # variance_stabilizing_transform(dds)
    # ds = get_deseqstats(dds, "IN10018, 24h", "Control, 24h")
    # p_value_dict = make_p_value_dict(dds, control_name,
    #                                  test_list)
    # volcano_plot(ds.p_values, lfc, 20, 5)
    multidict = get_multidict(dds, designstring, control_name,
                                     test_list, "multidict")
    dict_met_volcs = get_multi_volcano(multidict, 20, 5, True)
    # print(dict_met_volcs)
    dict_met_volcs["IN10018, 24h"]["volcano"].show_plot()
