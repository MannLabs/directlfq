import os
import re
from io import StringIO
import itertools
import pandas as pd
import numpy as np
from scipy import stats
import directlfq.gui_textfields as gui_textfields
import panel.pane
import panel as pn
import directlfq.utils as lfq_utils



class BaseWidget(object):

    def __init__(self, name):
        self.name = name
        self.__update_event = pn.widgets.IntInput(value=0)
        self.depends = pn.depends(self.__update_event.param.value)
        self.active_depends = pn.depends(
            self.__update_event.param.value,
            watch=True
        )

    def trigger_dependancy(self):
        self.__update_event.value += 1


class HeaderWidget(object):
    """This class creates a layout for the header of the dashboard with the name of the tool and all links to the MPI website, the MPI Mann Lab page and the GitHub repo.

    Parameters
    ----------
    title : str
        The name of the tool.

    Attributes
    ----------
    header_title : pn.pane.Markdown
        A Panel Markdown pane that returns the title of the tool.
    mpi_biochem_logo : pn.pane.PNG
        A Panel PNG pane that embeds a png image file of the MPI Biochemisty logo and makes the image clickable with the link to the official website.
    mpi_logo : pn.pane.JPG
        A Panel JPG pane that embeds a jpg image file of the MPI Biochemisty logo and makes the image clickable with the link to the official website.
    github_logo : pn.pane.PNG
        A Panel PNG pane that embeds a png image file of the GitHub logo and makes the image clickable with the link to the GitHub repository of the project.

    """

    def __init__(
        self,
        title,
        img_folder_path,
        github_url
    ):
        self.header_title = pn.pane.Markdown(
            f'# {title}',
            sizing_mode='stretch_width',
        )
        self.biochem_logo_path = os.path.join(
            img_folder_path,
            "mpi_logo.png"
        )
        self.mpi_logo_path = os.path.join(
            img_folder_path,
            "max-planck-gesellschaft.jpg"
        )
        self.github_logo_path = os.path.join(
            img_folder_path,
            "github.png"
        )

        self.mpi_biochem_logo = pn.pane.PNG(
            self.biochem_logo_path,
            link_url='https://www.biochem.mpg.de/mann',
            width=60,
            height=60,
            align='start'
        )
        self.mpi_logo = pn.pane.JPG(
            self.mpi_logo_path,
            link_url='https://www.biochem.mpg.de/en',
            height=62,
            embed=True,
            width=62,
            margin=(5, 0, 0, 5),
            css_classes=['opt']
        )
        self.github_logo = pn.pane.PNG(
            self.github_logo_path,
            link_url=github_url,
            height=70,
            align='end'
        )

    def create(self):
        return pn.Row(
            self.mpi_biochem_logo,
            self.mpi_logo,
            self.header_title,
            self.github_logo,
            height=73,
            sizing_mode='stretch_width'
        )


class MainWidget(object):

    def __init__(
        self,
        description,
        manual_path
    ):
        self.project_description = pn.pane.Markdown(
            description,
            margin=(10, 0, 10, 0),
            css_classes=['main-part'],
            align='start'
        )
        self.manual = pn.widgets.FileDownload(
            file=manual_path,
            label='Download Manual',
            button_type='default',
            align='center',
            auto=True,
            height=31,
            width=200,
            margin=(0, 20, 0, 0)
        )

    def create(self):
        LAYOUT = pn.Row(
            self.project_description,
            pn.layout.HSpacer(width=500),
           # self.manual,
            background='#eaeaea',
            align='center',
            sizing_mode='stretch_width',
            height=190,
            margin=(10, 8, 10, 8),
            css_classes=['background']
        )
        return LAYOUT


class RunPipeline(BaseWidget):

    def __init__(self):
        super().__init__(name="Data")
        # DATA FILES
        self.path_analysis_file = pn.widgets.TextInput(
            name='Specify an analysis file:',
            placeholder='Enter the whole path to the MQ | Spectronaut | DIA-NN output file',
            width=900,
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.path_output_folder = pn.widgets.TextInput(
            name='Specify a path to the output folder:',
            placeholder='Enter the whole path to the output folder',
            width=900,
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        self.samplemap_title = pn.pane.Markdown(
            'Load an experiments-to-conditions file:',
            margin=(10,0,0,15)
        )
        self.samplemap = pn.widgets.FileInput(
            accept='.tsv,.csv,.txt',
            margin=(-10,0,5,15)
        )
        self.samplemap_table = pn.widgets.Tabulator(
            layout='fit_data_fill',
            height=300,
            show_index=False,
            width=570,
            align='center',
            margin=(15, 12, 10, 18)
        )
        self.assign_cond_pairs = pn.widgets.CrossSelector(
            width=870,
            height=300,
            align='center',
            margin=(15, 15, 15, 15)
        )
        # RUN PIPELINE
        self.run_pipeline_button = pn.widgets.Button(
            name='Run pipeline',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.run_pipeline_progress = pn.indicators.Progress(
            active=False,
            bar_color='light',
            width=250,
            align='center',
            margin=(-10, 0, 20, 0)
        )
        self.visualize_data_button = pn.widgets.Button(
            name='Visualize data',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(40, 0, 0, 0)
        )
        self.run_pipeline_error = pn.pane.Alert(
            width=600,
            alert_type="danger",
            # object='test warning message',
            margin=(-20, 10, -5, 16),
        )


    def create(self):
        self.layout = pn.Card(
            pn.Row(
                pn.Column(
                    self.path_analysis_file,
                    self.run_pipeline_error,
                    self.path_output_folder,
                    pn.Row(
                        pn.Column(
                            self.samplemap_title,
                            self.samplemap
                        ),
                        pn.Card(
                            self.samplemap_table,
                            header='Experiment to condition map',
                            collapsed=False,
                            margin=(20, 0, 20, 0),
                            width=601,
                        )
                    ),
                    pn.Card(
                        self.assign_cond_pairs,
                        header='Select condition pairs for the analysis',
                        collapsed=True,
                        margin=(5, 0, 20, 15),
                        width=901,
                    ),
                    margin=(20, 30, 10, 10),
                ),
                pn.Spacer(sizing_mode='stretch_width'),
                pn.Column(
                    gui_textfields.Descriptions.project_instruction,
                    gui_textfields.Cards.spectronaut,
                    gui_textfields.Cards.diann,
                    gui_textfields.Cards.alphapept,
                    gui_textfields.Cards.maxquant,
                    self.run_pipeline_button,
                    self.run_pipeline_progress,
                    self.visualize_data_button,
                    align='center',
                    margin=(100, 40, 0, 0),
                )
            ),
            title='Run Pipeline | Visualize data',
            collapsed=False,
            header_background='#eaeaea',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            #height=1000,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )
        self.path_analysis_file.param.watch(
            self.activate_after_analysis_file_upload,
            'value'
        )
        self.samplemap.param.watch(
            self.update_samplemap,
            'value'
        )
        self.samplemap_table.param.watch(
            self.add_conditions_for_assignment,
            'value'
        )
        self.run_pipeline_button.param.watch(
            self.run_pipeline,
            'clicks'
        )
        self.visualize_data_button.param.watch(
            self.visualize_data,
            'clicks'
        )
        return self.layout

    def activate_after_analysis_file_upload(self, *args):
        self.set_default_output_folder()
        self.import_exp_data()
        self.extract_sample_names()

    def natural_sort(self, l):
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(l, key = alphanum_key)

    def set_default_output_folder(self):
        if not self.path_output_folder.value:
            self.path_output_folder.value = os.path.join(
                os.path.dirname(self.path_analysis_file.value),
                'results'
            )

    def import_exp_data(self):
        self.data = aqutils.import_data(
            input_file=self.path_analysis_file.value
        )

    def extract_sample_names(self):
        sample_names = aqutils.get_samplenames_from_input_df(self.data)
        self.samplemap_table.value = pd.DataFrame(
            data={'sample': self.natural_sort(sample_names), 'condition': str()}
        )

    def update_samplemap(self, *args):
        file_ext = os.path.splitext(self.samplemap.filename)[-1]

        if file_ext=='.csv':
            sep=','
        else:
            sep='\t'

        self.samplemap_table.value = pd.read_csv(
            StringIO(str(self.samplemap.value, "utf-8")),
            sep=sep,dtype='str'
        )

    def add_conditions_for_assignment(self, *args):
        unique_condit = self.samplemap_table.value.condition.unique()
        comb_condit = ['_vs_'.join(comb) for comb in itertools.permutations(unique_condit, 2)]
        self.assign_cond_pairs.options = comb_condit

    def run_pipeline(self, *args):
        self.run_pipeline_progress.active = True

        data_processed, samplemap_df_processed = aqutils.prepare_loaded_tables(
            self.data,
            self.samplemap_table.value
        )

        if self.assign_cond_pairs.value:
            cond_combinations = [tuple(pair.split('_vs_')) for pair in self.assign_cond_pairs.value]
        else:
            cond_combinations = [tuple(pair.split('_vs_')) for pair in self.assign_cond_pairs.options]
        
        diffmgr.run_pipeline(input_file=self.path_analysis_file.value,
            samplemap_df=self.samplemap_table.value,
            results_dir=self.path_output_folder.value,
            condpair_combinations=cond_combinations
        )

        self.trigger_dependancy()
        self.run_pipeline_progress.active = False

    def visualize_data(self, *args):
        self.trigger_dependancy()


class Tabs(object):

    def __init__(self, pipeline):
        self.layout = None
        self.pipeline = pipeline

    def create(
        self,
        # tab_list=None
    ):
        # self.tabs = tab_list
        return self.pipeline.depends(self.create_layout)

    def create_layout(self, *args):
        if self.pipeline.path_output_folder.value is not None and self.pipeline.visualize_data_button.clicks != 0:
            self.pipeline.layout.collapsed = True
            self.layout = pn.Tabs(
                tabs_location='above',
                margin=(30, 10, 5, 8),
                sizing_mode='stretch_width',
            )
            # self.layout += self.tabs
            self.layout += [
                # ('Multiple Comparison', MultipleComparison(  #Commented out for pre-release
                #     self.pipeline.path_output_folder.value
                #     ).create()
                # ),
                ('Single Comparison', SingleComparison(
                    self.pipeline.path_output_folder.value,
                    self.pipeline.samplemap_table.value,
                    ).create()
                ),
            ]
            self.active = 0
            return self.layout


class MultipleComparison(object):

    def __init__(self, output_folder):
        self.condpairs_to_compare = pn.widgets.CrossSelector(
            width=870,
            height=300,
            align='center',
            margin=(15, 15, 15, 15)
        )
        self.output_folder = output_folder
        self.layout = None
        self.heatmap = None

    def create(self):
        self.extract_conditions_from_folder()
        self.condpairs_to_compare.param.watch(
            self.return_clustered_heatmap,
            'value'
        )
        self.layout = pn.Column(
            self.condpairs_to_compare,
            self.heatmap
        )
        return self.layout

    def extract_conditions_from_folder(self):
        self.condpairs_to_compare.options = [f.replace(".results.tsv", "").replace('VS', 'vs') for f in os.listdir(self.output_folder) if re.match(r'.*results.tsv', f)]

    def return_clustered_heatmap(self, *args):
        if self.condpairs_to_compare.value:
            cond_combinations = [tuple(pair.split('_vs_')) for pair in self.condpairs_to_compare.value]
        else:
            cond_combinations = [tuple(pair.split('_vs_')) for pair in self.condpairs_to_compare.options]
        
        if len(cond_combinations)>1:
            overview_dataframe = aqplot.get_sample_overview_dataframe(
                results_folder=self.output_folder,
                condpairs_to_compare=cond_combinations
            )

            clustered_dataframe = aqplot.get_clustered_dataframe(overview_dataframe)

            self.layout[1] = pn.Pane(
                self.plot_heatmap(
                    clustered_dataframe.T,
                    title='Significant proteins heatmap',
                    colormap='RdBu'
                )
            )

    def plot_heatmap(self, df, title, colormap):
        """
        This function plots a simple Heatmap.
        """
        print('inside plotting function')
        fig = go.Figure()
        fig.update_layout(
            title={
                'text': title,
                'y': 0.95,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            height=1000,
            width=700,
            annotations=[dict(xref='paper', yref='paper', showarrow=False, text='')],
            template='plotly_white'
        )
        fig.add_trace(
            go.Heatmap(
                z=df.values.tolist(),
                x = list(df.columns),
                y = list(df.index),
                colorscale=colormap
            )
        )
        return fig


class SingleComparison(object):

    def __init__(self, output_folder, sample_to_cond):
        self.condpairs_selector = pn.widgets.Select(
            name='Select a pair of conditions:',
            options=['No conditions'],
            width=870,
            align='center',
            margin=(15, 15, 15, 15)
        )
        self.protein = pn.widgets.AutocompleteInput(
            name='Select protein:',
            placeholder="Type the first letters of the protein's name",
            min_characters=2,
            disabled=True,
            margin=(0, 40, 40, 10)
        )
        
        self.output_folder = output_folder
        self.sample_to_cond = sample_to_cond
        self.layout = None
        self.volcano_plot = None

    def create(self):
        self.extract_conditions_from_folder()
        self.condpairs_selector.param.watch(
            self.run_after_pair_cond_selection,
            'value'
        )
        self.protein.param.watch(
            self.visualize_after_protein_selection,
            'value'
        )
        self.layout = pn.Column(
            self.condpairs_selector,
            None,
            self.volcano_plot,
            pn.layout.Divider(sizing_mode='stretch_width'),
            pn.Column(
                self.protein,
                None,
                None,
                sizing_mode = 'stretch_width'
                # pn.Row(
                #     None,
                #     None,
                #     sizing_mode='stretch_width'
                # )
            ),

        )

        return self.layout

    def extract_conditions_from_folder(self):
        self.condpairs_selector.options.extend([f.replace(".results.tsv", "").replace('VS', 'vs') for f in os.listdir(self.output_folder) if re.match(r'.*results.tsv', f)])


    def run_after_pair_cond_selection(self, *args):
        if self.condpairs_selector.value != 'No conditions':
            self.cond1, self.cond2 = self.condpairs_selector.value.split('_vs_')
            self.iontree_condpair = aqutils.read_condpair_tree(self.cond1, self.cond2, self.output_folder)
            self.result_df = aqplot.get_diffresult_dataframe(
                self.cond1,
                self.cond2,
                results_folder=self.output_folder
            )
            self.normalized_intensity_df = aqplot.get_normed_peptides_dataframe(
                self.cond1,
                self.cond2,
                results_folder=self.output_folder
            )

            quantinfo = aqplot.CondpairQuantificationInfo().init_from_loaded_tables(diffresults_df=self.result_df, normed_df=self.normalized_intensity_df, condpair_root_node= self.iontree_condpair,samplemap_df=self.sample_to_cond)
            self.protein_df_getter = aqplot.ProteinIntensityDataFrameGetter(quantinfo)

            self.sample2cond = dict(zip(self.sample_to_cond["sample"], self.sample_to_cond["condition"]))

            c1_normed = aqplot.subset_normed_peptides_df_to_condition(
                self.cond1,
                self.sample_to_cond,
                self.normalized_intensity_df
            )
            c2_normed = aqplot.subset_normed_peptides_df_to_condition(
                self.cond2,
                self.sample_to_cond,
                self.normalized_intensity_df
            )

            self.protein.options = self.result_df.protein.values.tolist()
            self.protein.disabled = False

            self.volcano_plot = pn.pane.Plotly(
                aqplot.plot_volcano_plotly(
                    self.result_df
                ),
                align='center',
                config={'responsive': True}
            )

            self.volcano_plot.param.watch(
                self.update_selected_protein,
                'click_data'
            )

            self.layout[1] = pn.Card(
                pn.Pane(
                    aqplot.plot_withincond_fcs_plotly(
                        self.normalized_intensity_df,
                        'Log2 fold change distribution between all samples'
                    ),
                    align='center'

                ),
                pn.Pane(
                    aqplot.plot_withincond_fcs_plotly(
                        c1_normed,
                        f'Log2 fold change distribution between samples in condition {self.cond1}'
                    ),
                    align='center'
                ),
                pn.Pane(
                    aqplot.plot_withincond_fcs_plotly(
                        c2_normed,
                        f'Log2 fold change distribution between samples in condition {self.cond2}'
                    ),
                    align='center'
                ),
                pn.Pane(
                    aqplot.plot_betweencond_fcs_plotly(
                        c1_normed,
                        c2_normed,
                        'The distribution of direct peptide fold changes between conditions'
                    ),
                    align='center'
                ),
                header='Log2 fold change distribution plots',
                collapsed=True,
                align='center',
                sizing_mode='stretch_width'
            )
            self.layout[2] = self.volcano_plot

        else:
            self.layout[1] = None
            self.layout[2] = None

    def update_selected_protein(self, event):
        self.protein.value = str(event[2].click_data['points'][0]['text'])

    def visualize_after_protein_selection(self, *args):
        melted_df = self.protein_df_getter.get_melted_protein_ion_intensity_table(self.protein.value)
        self.protein_df = self.protein_df_getter.get_protein_diffresults(self.protein.value)

        #get fold change data table for protein
        self.fc_df = aqplot.get_betweencond_fcs_table(
            melted_df,
            self.cond1,
            self.cond2
        )

        if self.iontree_condpair != None:
            self.protnode = anytree.find(self.iontree_condpair, filter_= lambda x : x.name == self.protein.value, maxlevel=2)
        else:
            self.protnode = None
        
        self.layout[4][1] = pn.Pane(
            aqplot.beeswarm_ion_plot_plotly(melted_df, self.protein_df),
        )
        self.layout[4][2] = pn.Pane(
            aqplot.foldchange_ion_plot_plotly(self.fc_df, self.protein_df, protein_node=self.protnode),
        )



