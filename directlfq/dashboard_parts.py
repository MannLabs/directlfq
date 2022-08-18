import os
import re
import itertools

from click import option
import directlfq.gui_textfields as gui_textfields
import panel as pn
import directlfq.lfq_manager as lfq_manager


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
        description
    ):
        self.project_description = pn.pane.Markdown(
            description,
            margin=(10, 0, 10, 0),
            css_classes=['main-part'],
            align='start'
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

        self.dropdown_menu_for_input_type = pn.widgets.Select(name = "Type of input table", options = ['MaxQuant evidence.txt', 'MaxQuant peptides.txt',
        'Spectronaut fragment level', 'Spectronaut precursor level', 'DIANN fragment level', 'DIANN precursor level'])
        
        # self.path_output_folder = pn.widgets.TextInput(
        #     name='(optional) Specify a path to the output folder:',
        #     placeholder='(optional) Enter the whole path to the output folder',
        #     width=900,
        #     sizing_mode='stretch_width',
        #     margin=(15, 15, 0, 15)
        # )
        
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
        # self.visualize_data_button = pn.widgets.Button(
        #     name='Visualize data',
        #     button_type='primary',
        #     height=31,
        #     width=250,
        #     align='center',
        #     margin=(40, 0, 0, 0)
        # )
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
                    self.dropdown_menu_for_input_type,
                    self.run_pipeline_error,
                    #self.path_output_folder,

                    gui_textfields.Descriptions.project_instruction,
                    gui_textfields.Cards.spectronaut,
                    gui_textfields.Cards.diann,
                    gui_textfields.Cards.alphapept,
                    gui_textfields.Cards.maxquant,
                    self.run_pipeline_button,
                    self.run_pipeline_progress,
                   # self.visualize_data_button,
                    align='center',

                    #margin=(100, 30, 10, 10),
                ),
                pn.Spacer(sizing_mode='stretch_width'),
                # pn.Column(
                #     gui_textfields.Descriptions.project_instruction,
                #     gui_textfields.Cards.spectronaut,
                #     gui_textfields.Cards.diann,
                #     gui_textfields.Cards.alphapept,
                #     gui_textfields.Cards.maxquant,
                #     self.run_pipeline_button,
                #     self.run_pipeline_progress,
                #    # self.visualize_data_button,
                #     align='center',
                #     margin=(100, 40, 0, 0),
                # )
            ),
            title='Run directLFQ analysis',
            collapsed=False,
            header_background='#eaeaea',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            #height=1000,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )
        # self.path_analysis_file.param.watch(
        #     self.activate_after_analysis_file_upload,
        #     'value'
        # )
        self.run_pipeline_button.param.watch(
            self.run_pipeline,
            'clicks'
        )
        # self.visualize_data_button.param.watch(
        #     self.visualize_data,
        #     'clicks'
        # )
        return self.layout

    # def activate_after_analysis_file_upload(self, *args):
    #     self.run_pipeline()


    def run_pipeline(self, *args):
        self.run_pipeline_progress.active = True
        input_file = self.path_analysis_file.value
        input_type_to_use = self.get_intable_config_id_from_selected_input_type()

        lfq_manager.run_lfq(input_file = input_file, input_type_to_use = input_type_to_use, min_nonan = 1, maximum_number_of_quadratic_ions_to_use_per_protein = 10, number_of_quadratic_samples = 50)

        self.trigger_dependancy()
        self.run_pipeline_progress.active = False

    def get_intable_config_id_from_selected_input_type(self):
        input_list = ['MaxQuant evidence.txt', 'MaxQuant peptides.txt',
        'Spectronaut fragment level', 'Spectronaut precursor level', 'DIANN fragment level', 'DIANN precursor level']
        intable_ids = ['maxquant_evidence_proteins_column', 'maxquant_peptides', 'spectronaut_fragion_isotopes', 'spectronaut_precursor','diann_fragion_isotopes',
        'diann_precursor']
        mapping =  dict(zip(input_list, intable_ids))
        selected_value_by_user = self.dropdown_menu_for_input_type.value
        return mapping.get(selected_value_by_user)
 


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


