import GEOparse
import pandas
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import numpy as np
from math import log10


class Volcano_obj:
    """
    Object used to store,show and save volcano plots.
    :param dict  conditiondict: Dictionary containing the P values,
        log fold changes and base means of a test case.
    :param str condition_name: Name of the test case.
    :param int p_limit: Maximum limit to which higher p-values wil be reduced.
    :param int fold_limit: Maximum limit to which higher and
        lower fold values wil be reduced.
    :param str save_string: Name of file where volcano plot will be saved,
        nothing is saved if empty.
    """

    def __init__(
            self,
            conditiondict,
            condition_name,
            p_limit: int,
            fold_limit: int,
            save_string=None
    ):
        self.conditiondict = conditiondict
        self.condition_name = condition_name
        self.p_limit = p_limit
        self.fold_limit = fold_limit
        self.transformed_p_values = self.transform_p_values()
        self.transformed_folds = self.transform_folds()
        if save_string is not None:
            self.save_string = save_string
            self.save_plot()


    def transform_p_values(self)->list:
        """
        Creates a modified P values list that
        does not exceed the P value limit.

        :return: Modified P values list that does not exceed the P value limit.
        :rtype: List
        """
        modified_p_values_list = []
        for x in self.conditiondict["p_values"]:
            if (log10(x) * -1) > self.p_limit:
                modified_p_values_list.append(self.p_limit)
            else:
                modified_p_values_list.append(log10(x) * -1)
        return modified_p_values_list

    def transform_folds(self)->pandas.Series:
        """
        Creates a modified log fold change list that
        does not exceed the log fold value limit.

        :return: Modified log fold change list that does not exceed
        the log fold value limit.
        :rtype: List
        """
        numpy_fold = self.conditiondict["fold_change"].to_numpy()
        modified_fold_list = []
        for x in numpy_fold:
            if x > self.fold_limit:
                modified_fold_list.append(self.fold_limit)
            elif x < -self.fold_limit:
                modified_fold_list.append(-self.fold_limit)
            else:
                modified_fold_list.append(x)
        modified_fold = pandas.Series(modified_fold_list)
        return modified_fold


    def show_plot(self):
        """
        Show the volcano plot in a pop-up window.
        """
        plt.scatter(self.transformed_folds, self.transformed_p_values, marker='.')
        plt.xlabel('FoldChange(log2)')
        plt.ylabel('P-value(-10log)')
        plt.title(self.condition_name)
        plt.show()
        plt.clf()


    def save_plot(self, new_save_string = None):
        """
        Saves the volcano plot to a file. save_string is used as
        file name if no new file name is provided.

        :param str new_save_string: New file name to overwrite the
            old file name.
        """
        plt.scatter(self.transformed_folds, self.transformed_p_values, marker='.')
        plt.xlabel('FoldChange(log2)')
        plt.ylabel('P-value(-10log)')
        plt.title(self.condition_name)
        if new_save_string is not None:
            self.save_string = new_save_string
        plt.savefig(self.save_string)
        plt.clf()


    def __str__(self)->str:
        return ("Volcano plot: \n"
                "Condition name: %s, p_limit: %s,"
                " fold_limit: %s."
                % (self.condition_name, self.p_limit,
                   self.fold_limit))