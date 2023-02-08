from turtle import width
import panel.pane
import panel as pn
import os
import pathlib


class Paths():
    CONFIGS_PATH = os.path.join(pathlib.Path(__file__).parent.absolute(), "configs")
    spectronaut_fragion_path = os.path.join(CONFIGS_PATH, "spectronaut_tableconfig_fragion.rs")
    spectronaut_precursor_path = os.path.join(CONFIGS_PATH, "spectronaut_tableconfig_precursor.rs")
    
class ButtonConfiguration():
    width = 530


class Descriptions():


    project_instruction = panel.pane.Markdown("""
        #### How to use directLFQ:
        1. Provide the filepath to your proteomic datasets analyzed by
        Alphapept, MaxQuant, Spectronaut, DIA-NN or FragPipe (see specifications below!).
        2. (optional) If you are using MaxQuant evidence.txt or peptides.txt files, you can add the link to the corresponding proteinGroups.txt file. Adding the proteinGroups.txt will improve peptide mapping.
        3. (optional) You can configure your run with some of the optional parameters accessible in the "optional configurations" section.
        4. Click on the _RUN PIPELINE_ button, **check the progress on the terminal window (this window opens together with this GUI) and ensure that 'Starting directLFQ analysis' message has appeared.**
        5. Once the analysis is finished, a results file with LFQ intensities ending with .protein_intensities.tsv is written in the same folder as the input file.
        """,
        width=ButtonConfiguration.width,
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
        directLFQ takes a Spectronaut .tsv table as input. When exporting from Spectronaut, the correct columns need to be selected. 
        These can be obtained by downloading one of the export schemes available below. We provide one export scheme for precursor quantification
        and one export scheme for fragment ion quantification. Fragment ion quantification shows slightly more accuracy, but the files are around 10 times larger.
        
        An export scheme can then simply be loaded into Spectronaut as follows:
        
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
            Provide the path to the MaxQuant peptides.txt output table or the MaxQuant evidence.txt output table. Additionally and if possible, provide the path to the corresponding proteinGroups.txt file.
            """,
            width=ButtonConfiguration.width,
            align='start',
            margin=(0, 80, 0, 20)
        )

    fragpipe = pn.pane.Markdown(
           """
        Provide the path to the "combined_ion.tsv" output table.

        """
        ,
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 20)
    )




class DownloadSchemes():

    spectronaut_fragion = pn.widgets.FileDownload(
    file=Paths.spectronaut_fragion_path,
    filename="spectronaut_tableconfig_fragion.rs",
    button_type='default',
    auto=True,
    css_classes=['button_options'],
)
    spectronaut_precursor = pn.widgets.FileDownload(
    file=Paths.spectronaut_precursor_path,
    filename="spectronaut_tableconfig_precursor.rs",
    button_type='default',
    auto=True,
    css_classes=['button_options'],
)


class Cards():
    width = 530

    alphapept = pn.Card(
        Descriptions.alphapept,
        header='AlphaPept input specifications',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )


    spectronaut = pn.Card(
        Descriptions.spectronaut,
        DownloadSchemes.spectronaut_fragion,
        DownloadSchemes.spectronaut_precursor,
        header='Spectronaut input specifications',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 5, 10),
        css_classes=['spectronaut_instr']
    )
    diann = pn.Card(
        Descriptions.diann,
        header='DIANN input specifications',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )

    maxquant = pn.Card(
        Descriptions.maxquant,
        header='MaxQuant input specifications',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )

    fragpipe = pn.Card(
        Descriptions.fragpipe,
        header='FragPipe input specifications',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )



