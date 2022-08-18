from turtle import width
import panel.pane
import panel as pn
import os
import pathlib


class Paths():
    CONFIGS_PATH = os.path.join(pathlib.Path(__file__).parent.absolute(), "configs")
    spectronaut_fragion_path = os.path.join(CONFIGS_PATH, "spectronaut_tableconfig_fragion.rs")
    
class ButtonConfiguration():
    width = 530


class Descriptions():


    project_instruction = panel.pane.Markdown("""
        #### How to use directLFQ:
        1. Provide the filepath to your proteomic datasets analyzed by
        Spectronaut, DIA-NN or MaxQuant (see instructions below!).
        2. Select the type of the input table you have specified from the dropdown menu
        3. Click on the _RUN PIPELINE_ button, you can follow the progress on the terminal window
        """,
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 10))

    single_comparison_instruction = panel.pane.Markdown("""
        Here you can visualize comparisons of two conditions as a volcano plot. You can click or search proteins of interest and detail plots of the quantification will be shown. 
        The displayed data is stored as text files in the output folder you specified.
        """,
        width=830,
        align='start',
        margin=(0, 80, 0, 10))

    alphapept = pn.pane.Markdown(
        """
        Provide the path to the AlphaPept results_peptides.csv output table.

        """,
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 20)
    )

    spectronaut = pn.pane.Markdown(
        """
        To get the most out of the Spectronaut data, AlphaQuant utilizes more than 30 different columns.
        These can be obtained by downloading the export scheme "spectronaut_tableconfig_fragion.rs", 
        which can then simply be loaded into Spectronaut as follows: 
        
        Go to the "Report" perspective in Spectronaut, click "Import Schema" and provide the file.

        The data needs to be exported in the **normal long** format as .tsv or .csv file. 

        """,
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 20)
    )

    diann = pn.pane.Markdown(
            """
            Provide the path to the DIANN report.tsv output table.
            """,
            width=ButtonConfiguration.width,
            align='start',
            margin=(0, 80, 0, 20)
        )

    maxquant = pn.pane.Markdown(
            """
            Provide the path to the MaxQuant peptides.txt output table.
            """,
            width=ButtonConfiguration.width,
            align='start',
            margin=(0, 80, 0, 20)
        )




class DownloadSchemes():

    spectronaut = pn.widgets.FileDownload(
    file=Paths.spectronaut_fragion_path,
    filename="spectronaut_tableconfig_fragion.rs",
    button_type='default',
    auto=True,
    css_classes=['button_options'],
)


class Cards():
    width = 530

    alphapept = pn.Card(
        Descriptions.alphapept,
        header='AlphaPept instructions',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )


    spectronaut = pn.Card(
        Descriptions.spectronaut,
        DownloadSchemes.spectronaut,
        header='Spectronaut instructions',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 5, 10),
        css_classes=['spectronaut_instr']
    )
    diann = pn.Card(
        Descriptions.diann,
        header='DIANN instructions',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )

    maxquant = pn.Card(
        Descriptions.maxquant,
        header='MaxQuant instructions',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )



