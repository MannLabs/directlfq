import directlfq.config as config
import logging
config.setup_logging()

LOGGER = logging.getLogger(__name__)

a4_dims = (11.7, 8.27)
a4_width_no_margin = 10.5

# %% ../nbdev_nbs/05_visualizations.ipynb 2
import seaborn as sns
import matplotlib
import matplotlib.cm
import matplotlib.colors
import matplotlib.pyplot as plt

class AlphaPeptColorMap():
    def __init__(self):

        #colorlist = ["#3FC5F0", "#42DEE1", "#7BEDC5", "#FFD479", "#16212B"]
        colorlist = ["#3FC5F0","#16212B", "#FFD479", "#42DEE1", "#7BEDC5" ]
        self.colorlist = [matplotlib.colors.to_rgba(x) for x in colorlist]
        self.colorlist_hex = colorlist
        self.colormap_linear = matplotlib.colors.LinearSegmentedColormap.from_list("alphapept",self.colorlist)
        self.colormap_discrete = matplotlib.colors.LinearSegmentedColormap.from_list("alphapept",self.colorlist, N=5)
        self.seaborn_mapname_linear = None
        self.seaborn_mapname_linear_cut = None
    

class CmapRegistrator():
    def __init__(self):
        self._existing_matplotlib_cmaps = None
        self._define_existing_cmaps()
    
    def _define_existing_cmaps(self):
        self._existing_matplotlib_cmaps = [name for name in plt.colormaps() if not name.endswith('_r')]

    def register_colormap(self, name, colorlist):
        linmap = matplotlib.colors.LinearSegmentedColormap.from_list(name, colorlist)
        matplotlib.cm.register_cmap(name, linmap)
    

# %% ../nbdev_nbs/05_visualizations.ipynb 3
import pandas as pd
import directlfq.normalization as lfq_norm
import directlfq.protein_intensity_estimation as lfq_protint
import matplotlib.pyplot as plt

class IonTraceCompararisonPlotter():
    def __init__(self, proteome_df, selected_protein, axis_unnormed, axis_normed):
        self._proteome_df = proteome_df
        self._selected_protein = selected_protein
        self._protein_df_before_norm = None
        self._protein_df_after_norm = None
        
        self.axis_unnormed = axis_unnormed
        self.axis_normed = axis_normed

        self._prepare_data_and_plot_ion_traces_before_and_after_normalization()

    def _prepare_data_and_plot_ion_traces_before_and_after_normalization(self):
        self._define_protein_dataframes()
        self._plot_before_norm()
        self._plot_after_norm()

    def _define_protein_dataframes(self):
        self._define_protein_df_before_norm()
        self._define_protein_df_after_norm()

    def _define_protein_df_before_norm(self):
        self._protein_df_before_norm = pd.DataFrame(self._proteome_df.loc[self._selected_protein])
        self._protein_df_before_norm = self._protein_df_before_norm.dropna(axis='columns', how='all')
    
    def _define_protein_df_after_norm(self):
        self._protein_df_after_norm = lfq_norm.NormalizationManagerProtein(self._protein_df_before_norm.copy(), num_samples_quadratic = 10).complete_dataframe
        self._protein_df_after_norm = self._protein_df_after_norm.dropna(axis='columns', how='all')

    def _plot_before_norm(self):
        IonTraceVisualizer(self._protein_df_before_norm,ax= self.axis_unnormed)
    
    def _plot_after_norm(self):
        visualizer = IonTraceVisualizer(self._protein_df_after_norm, ax=self.axis_normed)
        median_list = lfq_protint.get_list_with_protein_value_for_each_sample(self._protein_df_after_norm, min_nonan=1)
        visualizer.add_median_trace(median_list)


class IonTraceCompararisonPlotterNoDirectLFQTrace(IonTraceCompararisonPlotter):
    def __init__(self, proteome_df, selected_protein, ax):
        self._proteome_df = proteome_df
        self._selected_protein = selected_protein
        self._protein_df_before_norm = None
        self._protein_df_after_norm = None
        
        self.axis_normed = ax
        
        self._prepare_data_and_plot_ion_traces_before_and_after_normalization()

    
    def _prepare_data_and_plot_ion_traces_before_and_after_normalization(self):
        self._define_protein_dataframes()
        self._plot_after_norm()

    def _plot_after_norm(self):
        visualizer = IonTraceVisualizer(self._protein_df_before_norm, ax=self.axis_normed)

        

# %% ../nbdev_nbs/05_visualizations.ipynb 4
import seaborn as sns
import matplotlib.cm
import numpy as np

class IonTraceVisualizer():
    def __init__(self, protein_df, ax):
        self._protein_df = protein_df
        self._plot_df = None
        self._num_samples = None
        self._ax = ax
        self._define_inputs_and_plot_ion_traces()
    
    def _define_inputs_and_plot_ion_traces(self):
        self._define_num_samples()
        self._define_prepared_dataframe()
        self._plot_ion_traces()

    def _define_num_samples(self):
        self._num_samples = len(self._protein_df.columns)

    def _plot_ion_traces(self):
        plot_values = self._plot_df.values #row contains intensity trace
        for idx in range(plot_values.shape[0]):
            x_values = np.array(range(plot_values.shape[1]))
            y_values = plot_values[idx]
            nan_mask = np.isfinite(y_values)
            self._ax.plot(x_values[nan_mask], y_values[nan_mask],color='grey', alpha=0.5)
            self._ax.scatter(x_values[nan_mask], y_values[nan_mask], color='grey', marker = 'o', s = 11)
        self._ax.set_xticks(range(self._num_samples))
        self._annotate_x_ticks(sample_names=self._protein_df.columns)

    def _define_prepared_dataframe(self):
        #drop all rows that contain less than 1 non nan value
        self._plot_df = self._protein_df.copy()
        self._plot_df = self._plot_df.dropna(axis='rows', thresh=1)
        self._plot_df.columns = range(self._num_samples)
        #self._plot_df = self._plot_df.T

    def add_median_trace(self, list_of_median_values):
        sns.lineplot(x = range(len(list_of_median_values)), y = list_of_median_values, ax=self._ax,color='black', linewidth=3)
    
    #function that annotates x ticks of an axis with the sample names
    def _annotate_x_ticks(self, sample_names):
        self._ax.set_xticklabels(sample_names, rotation=90)



# %% ../nbdev_nbs/05_visualizations.ipynb 5
import seaborn as sns

class MultiOrganismMultiMethodBoxPlot():
    def __init__(self, method_ratio_results_table, ax, organisms_to_plot, fcs_to_expect):
        self._method_ratio_results_table = method_ratio_results_table
        self._colorlist_hex =['#bad566', '#325e7a', '#ffd479'] + AlphaPeptColorMap().colorlist_hex
        self._fcs_to_expect = fcs_to_expect
        self._organisms_to_plot = organisms_to_plot
        
        self.ax = ax

        self.plot_boxplot()
        self._add_expected_fold_changes()

    def plot_boxplot(self):
        color_palette = sns.color_palette(self._colorlist_hex, n_colors=len(self._fcs_to_expect))
        sns.violinplot(data=self._method_ratio_results_table, x="method", y = "log2fc", hue= "organism", palette=color_palette, hue_order=self._organisms_to_plot, ax=self.ax)
    
    def _add_expected_fold_changes(self):
        if self._fcs_to_expect is not None:
            for idx, fc in enumerate(self._fcs_to_expect):
                color = self._colorlist_hex[idx]
                self.ax.axhline(fc, color = color)
    
    

# %% ../nbdev_nbs/05_visualizations.ipynb 6
import matplotlib.pyplot as plt
import itertools

def plot_withincond_fcs(normed_intensity_df, cut_extremes = True):
    """takes a normalized intensity dataframe and plots the fold change distribution between all samples. Column = sample, row = ion"""

    samplecombs = list(itertools.combinations(normed_intensity_df.columns, 2))

    for spair in samplecombs:#compare all pairs of samples
        s1 = spair[0]
        s2 = spair[1]
        diff_fcs = normed_intensity_df[s1].to_numpy() - normed_intensity_df[s2].to_numpy() #calculate fold changes by subtracting log2 intensities of both samples

        if cut_extremes:
            cutoff = max(abs(np.nanquantile(diff_fcs,0.025)), abs(np.nanquantile(diff_fcs, 0.975))) #determine 2.5% - 97.5% interval, i.e. remove extremes
            range = (-cutoff, cutoff)
        else:
            range = None
        plt.hist(diff_fcs,80,density=True, histtype='step',range=range) #set the cutoffs to focus the visualization
        plt.xlabel("log2 peptide fcs")

    plt.show()

# %% ../nbdev_nbs/05_visualizations.ipynb 7
import matplotlib.pyplot as plt
import itertools

def plot_relative_to_median_fcs(normed_intensity_df):

    median_intensities = normed_intensity_df.median(axis=1)
    median_intensities = median_intensities.to_numpy()
    
    diff_fcs = []
    for col in normed_intensity_df.columns:
        median_fcs = normed_intensity_df[col].to_numpy() - median_intensities
        diff_fcs.append(np.nanmedian(median_fcs))
    plt.hist(diff_fcs,80,density=True, histtype='step')
    plt.xlabel("log2 peptide fcs")
    plt.show()



import numpy as np
class MultiOrganismIntensityFCPlotter():
    def __init__(self, ax, resultstable_w_ratios, organisms_to_plot = None, fcs_to_expect = None, title = ""):
        LOGGER.info('init MultiOrganismIntensityFCPlotter')
        self.ax = ax
        self._color_list_hex = ['#ffd479', '#325e7a', '#bad566']
        self._resultstable_w_ratios = resultstable_w_ratios
        self._organism_column = resultstable_w_ratios.organism_column
        self._log2fc_column = resultstable_w_ratios.log2fc_column
        self._mean_intensity_column = resultstable_w_ratios.mean_intensity_column
        
        self._organisms_to_plot = self._get_organisms_to_plot(organisms_to_plot)
        self._fcs_to_expect = fcs_to_expect

        self._title = self._get_title(title)
        self._scatter_per_organism()
        self._add_expected_lines()

    def _get_organisms_to_plot(self, organisms_to_plot):
        if organisms_to_plot is not None:
            return organisms_to_plot
        else:
            return sorted(list(set(self._resultstable_w_ratios.formated_dataframe[self._organism_column].astype('str'))))
    
    def _get_title(self, title):
        if title !="":
            self._print_infos_about_data()
            return title
        return self._generate_title()

    def _print_infos_about_data(self):
        for organism in self._organisms_to_plot:
            subtable_organism = self._get_organism_subtable(organism)
            LOGGER.info(self._get_stats_of_organism(organism, subtable_organism))

    def _generate_title(self):
        for organism in self._organisms_to_plot:
            subtable_organism = self._get_organism_subtable(organism)
            title += self._get_stats_of_organism(organism, subtable_organism)
        return title

    def _scatter_per_organism(self):
        complete_table = self._resultstable_w_ratios.formated_dataframe.copy()
        complete_table[self._mean_intensity_column] = np.log2(complete_table[self._mean_intensity_column])
        complete_table = self._remove_omitted_organisms_from_table(complete_table)
        color_palette = sns.color_palette(self._color_list_hex, n_colors=len(self._organisms_to_plot))
        sns.scatterplot(data= complete_table, x =self._mean_intensity_column, y= self._log2fc_column, hue=self._organism_column, alpha=0.15, ax=self.ax, 
        hue_order=self._organisms_to_plot, palette=color_palette, size=0.2)
        self.ax.set_title(self._title)
    
    def _remove_omitted_organisms_from_table(self, complete_table):
        row_w_permitted_organism = [x in self._organisms_to_plot for x in complete_table["organism"]]
        return complete_table[row_w_permitted_organism]

    def _add_expected_lines(self):
        if self._fcs_to_expect is not None:
            for idx, fc in enumerate(self._fcs_to_expect):
                color = self._color_list_hex[idx]
                self.ax.axhline(fc, color = color)

    def _get_organism_subtable(self, organism):
        complete_table = self._resultstable_w_ratios.formated_dataframe
        return complete_table[complete_table[self._organism_column] == organism]
    
    def _get_stats_of_organism(self, organism, subtable_organism):
        fcs = subtable_organism[self._log2fc_column].to_numpy()
        fcs = fcs[np.isfinite(fcs)]
        median_fc = np.nanmedian(fcs)
        std_fc = np.nanstd(fcs)
        num_ratios = sum(~np.isnan(fcs))
        return f"{organism} num:{num_ratios} median_FC:{median_fc:.2} STD:{std_fc:.2}\n"