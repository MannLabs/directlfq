import os
import panel as pn
import directlfq.gui_textfields as gui_textfields
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
    """
    This class creates a layout for the header of the dashboard with the name of the tool and all links to the MPI website, the MPI Mann Lab page, and the GitHub repo.

    Parameters
    ----------
    title : str
        The name of the tool.

    Attributes
    ----------
    header_title : pn.pane.Markdown
        A Panel Markdown pane that returns the title of the tool.
    mpi_biochem_logo : pn.pane.PNG
        A Panel PNG pane that embeds a png image file of the MPI Biochemistry logo and makes the image clickable with the link to the official website.
    mpi_logo : pn.pane.JPG
        A Panel JPG pane that embeds a jpg image file of the MPI Biochemistry logo and makes the image clickable with the link to the official website.
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
            align='center',  # Center the title
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
            align='start',
            margin=(0, 10, 0, 0),
        )
        self.mpi_logo = pn.pane.JPG(
            self.mpi_logo_path,
            link_url='https://www.biochem.mpg.de/en',
            height=62,
            embed=True,
            width=62,
            margin=(5, 0, 0, 5),
            css_classes=['opt'],
            align='start',
        )
        self.github_logo = pn.pane.PNG(
            self.github_logo_path,
            link_url=github_url,
            height=70,
            align='end',
            margin=(0, 0, 0, 10),
        )

    def create(self):
        return pn.Row(
            self.mpi_biochem_logo,
            self.mpi_logo,
            self.header_title,
            self.github_logo,
            height=73,
            sizing_mode='stretch_width',
            align='center',
        )


class MainWidget(object):

    def __init__(
        self,
        description
    ):
        self.project_description = pn.pane.Markdown(
            description,
            margin=(10, 15, 10, 15),
            css_classes=['main-part'],
            sizing_mode='stretch_width',
        )

    def create(self):
        LAYOUT = pn.Row(
            self.project_description,
            align='center',
            sizing_mode='stretch_width',
            margin=(10, 15, 10, 15),
            #css_classes=['background']
        )
        return LAYOUT


class RunPipeline(BaseWidget):

    def __init__(self):
        super().__init__(name="Data")
        # DATA FILES
        self.path_analysis_file = pn.widgets.TextInput(
            name='Specify an analysis file (see detailed specifications below!):',
            placeholder='Enter the whole path to the AP | MQ | Spectronaut | DIA-NN | Fragpipe outputs according to the specifications below',
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )

        self.path_protein_groups_file = pn.widgets.TextInput(
            name='(optional) If you are using MaxQuant evidence.txt or peptides.txt files, you can add the link to the corresponding proteinGroups.txt file (will improve peptide-to-protein mapping)',
            placeholder='(optional) Enter the whole path to the MaxQuant proteinGroups.txt file',
            value='',
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )

        # Optional files
        self.additional_headers_title = pn.pane.Markdown(
            '* Add the names of columns that you want to keep in the directLFQ output file, separated by semicolons. Note that some basic additional columns such as gene names are always added to the output table by value.\nWARNING: Take care that columns you add are not ambiguous. For example, adding the peptide sequence column will not work, because there are multiple peptide sequences per protein.',
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.additional_headers = pn.widgets.TextInput(
            name='',
            placeholder='(optional) Enter the names of columns that you want to keep',
            sizing_mode='stretch_width',
            margin=(0, 15, 0, 15)
        )

        self.protein_subset_for_normalization_title = pn.pane.Markdown(
            '* Specify a list of proteins (no header, separated by linebreaks) that you want to use for normalization. This could, for example, be a list of housekeeping proteins:',
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.protein_subset_for_normalization_file = pn.widgets.TextInput(
            name='',
            value='',
            placeholder='(optional) Enter the whole path to the protein list file',
            sizing_mode='stretch_width',
            margin=(0, 15, 0, 15)
        )

        self.yaml_filt_dict_title = pn.pane.Markdown(
            '* In case you want to define specific filters in addition to the standard filters, you can add a yaml file where the filters are defined (see GitHub docs).',
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.yaml_filt_dict_path = pn.widgets.TextInput(
            name='',
            value='',
            placeholder='(optional) Enter the whole path to the yaml file with the filters',
            sizing_mode='stretch_width',
            margin=(0, 15, 0, 15)
        )

        self.dropdown_menu_for_input_type_title = pn.pane.Markdown(
            '* Specify the type of the input table you want to use from the dropdown menu. Applies only if you want to use non-default settings, for example if you want to use summarized precursor intensities instead of fragment ion intensities for DIA data:',
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.dropdown_menu_for_input_type = pn.widgets.Select(
            name='',
            options={
                'detect automatically': None,
                'Alphapept peptides.csv': 'alphapept_peptides',
                'MaxQuant evidence.txt': "maxquant_evidence",
                'MaxQuant peptides.txt': 'maxquant_peptides',
                'Spectronaut fragment level': 'spectronaut_fragion_isotopes',
                'Spectronaut precursor level': 'spectronaut_precursor',
                'DIANN fragment level': 'diann_fragion_isotopes',
                'DIANN fragment level raw': 'diann_fragion_isotopes_raw',
                'DIANN precursor level': 'diann_precursors',
                'DIANN precursor level MS1': 'diann_precursors_ms1',
                'DIANN precursor MS1 and MS2': 'diann_precursor_ms1_and_ms2'
            },
            sizing_mode='stretch_width',
            margin=(0, 15, 0, 15)
        )

        self.num_nonan_vals_title = pn.pane.Markdown(
            '* Specify the minimum number of non-NaN ion intensities required to derive a protein intensity. The higher this number, the more reliable the protein quantification at the cost of more missing values:',
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.num_nonan_vals = pn.widgets.IntInput(
            name='',
            value=1,
            step=1,
            start=0,
            end=1000,
            width=100,
            margin=(0, 15, 0, 15)
        )

        self.num_cores_title = pn.pane.Markdown(
            '* Specify the number of cores to use (default of 0 means multiprocessing):',
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.num_cores_vals = pn.widgets.IntInput(
            name='',
            value=0,
            step=1,
            start=0,
            end=1000,
            width=100,
            margin=(0, 15, 0, 15)
        )

        # RUN PIPELINE
        self.run_pipeline_button = pn.widgets.Button(
            name='Run pipeline',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            #margin=(20, 0, 0, 0)
        )
        self.run_pipeline_progress = pn.indicators.Progress(
            active=False,
            bar_color='light',
            width=250,
            align='center',
            #margin=(-10, 0, 20, 0)
        )
        self.run_pipeline_error = pn.pane.Alert(
            width=600,
            alert_type="danger",
            margin=(-20, 10, -5, 16),
        )
        self.run_pipeline_success = pn.pane.Alert(
            width=600,
            alert_type="success",
            margin=(-20, 10, -5, 16),
        )

    def create(self):
        optional_configs = pn.Card(
            self.additional_headers_title,
            self.additional_headers,
            self.dropdown_menu_for_input_type_title,
            self.dropdown_menu_for_input_type,
            self.protein_subset_for_normalization_title,
            self.protein_subset_for_normalization_file,
            self.num_nonan_vals_title,
            self.num_nonan_vals,
            self.num_cores_title,
            self.num_cores_vals,
            self.yaml_filt_dict_title,
            self.yaml_filt_dict_path,
            header='Optional Configurations',
            collapsed=True,
            header_background='#eaeaea',
            header_color='#333',
            sizing_mode='stretch_width',
            margin=(10, 15, 0, 15)
        )

        # Wrap the specifications in a column for consistent alignment
        specifications = pn.Column(
            gui_textfields.Descriptions.project_instruction,
            gui_textfields.Cards.spectronaut,
            gui_textfields.Cards.diann,
            gui_textfields.Cards.alphapept,
            gui_textfields.Cards.maxquant,
            gui_textfields.Cards.fragpipe,
            sizing_mode='stretch_width',
            margin=(10, 15, 0, 15),
        )

        self.layout = pn.Card(
            pn.Column(
                self.path_analysis_file,
                self.path_protein_groups_file,
                optional_configs,
                specifications,
                pn.Row(
                    pn.layout.HSpacer(),
                    self.run_pipeline_button,
                    pn.layout.HSpacer(),
                    sizing_mode='stretch_width',
                    align='center',
                ),
                self.run_pipeline_progress,
                self.run_pipeline_success,
                align='start',
                sizing_mode='stretch_width',
                margin=(0, 15, 15, 15)
            ),
            title='Run directLFQ analysis',
            collapsible=False,  # Make the card non-collapsible
            header_background='#eaeaea',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        self.run_pipeline_button.on_click(self.run_pipeline)
        return self.layout

    
    def run_pipeline(self, *args):
        self.run_pipeline_progress.active = True
        self.run_pipeline_success.object = ""  # Clear any previous success message
        input_file = self.path_analysis_file.value_input
        input_type_to_use = self.dropdown_menu_for_input_type.value
        mq_protein_groups_txt = None if self.path_protein_groups_file.value_input == '' else self.path_protein_groups_file.value_input
        additional_headers = None if self.additional_headers.value_input == '' else self.additional_headers.value_input
        min_nonan = self.num_nonan_vals.value
        file_of_proteins_for_normalization = None if self.protein_subset_for_normalization_file.value_input == '' else self.protein_subset_for_normalization_file.value_input
        num_cores = None if self.num_cores_vals.value == -1 else self.num_cores_vals.value
        yaml_filt_dict_path = None if self.yaml_filt_dict_path.value == '' else self.yaml_filt_dict_path.value
        if additional_headers is not None: #the user will enter a string with semicolon separated values
            additional_headers = additional_headers.split(';')

        try:
            lfq_manager.run_lfq(
                input_file=input_file,
                input_type_to_use=input_type_to_use,
                maximum_number_of_quadratic_ions_to_use_per_protein=10,
                number_of_quadratic_samples=50,
                mq_protein_groups_txt=mq_protein_groups_txt,
                columns_to_add=additional_headers,
                selected_proteins_file=file_of_proteins_for_normalization,
                min_nonan=min_nonan,
                num_cores=num_cores,
                filter_dict=yaml_filt_dict_path
            )
            # Show success message
            self.run_pipeline_success.object = "✅ **Pipeline completed successfully!** Analysis has finished and output files have been generated."
        except Exception as e:
            # Show error message if something goes wrong
            self.run_pipeline_error.object = f"❌ **Error:** {str(e)}"

        self.trigger_dependancy()
        self.run_pipeline_progress.active = False



class Tabs(object):

    def __init__(self, pipeline):
        self.layout = None
        self.pipeline = pipeline

    def create(
        self,
    ):
        return self.pipeline.depends(self.create_layout)

    def create_layout(self, *args):
        if self.pipeline.path_output_folder.value is not None and self.pipeline.visualize_data_button.clicks != 0:
            self.pipeline.layout.collapsed = True
            self.layout = pn.Tabs(
                tabs_location='above',
                margin=(30, 10, 5, 8),
                sizing_mode='stretch_width',
            )
            self.layout += [
                ('Single Comparison', SingleComparison(
                    self.pipeline.path_output_folder.value,
                    self.pipeline.samplemap_table.value,
                ).create()),
            ]
            self.active = 0
            return self.layout
