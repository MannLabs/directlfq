#!python


# external
import click
import os
import json
import copy
import contextlib

# local
import directlfq
import directlfq.utils


with open(f"{directlfq.utils.CONFIG_PATH}/interface_parameters.json") as in_file:
    INTERFACE_PARAMETERS = json.load(in_file)


@contextlib.contextmanager
def parse_cli_settings(command_name: str, **kwargs):
    """A context manager that parses and logs CLI settings.
    Parameters
    ----------
    command_name : str
        The name of the command that utilizes these CLI settings.
    **kwargs
        All values that need to be logged.
        Values (if included) that explicitly will be parsed are:
            output_folder
            threads
            disable_log_stream
            log_file
            parameter_file
            export_parameters
    Returns
    -------
    : dict
        A dictionary with parsed parameters.
    """
    try:
        kwargs = {key: arg for key, arg in kwargs.items() if arg is not None}
        yield kwargs
    except Exception:
        print("Something went wrong, execution incomplete!")





def cli_option(
    parameter_name: str,
    as_argument: bool = False,
    **kwargs
):
    """A wrapper for click.options and click.arguments using local defaults.
    Parameters
    ----------
    parameter_name : str
        The name of the parameter or argument.
        It's default values need to be present in
        lib/interface_parameters.json.
    as_argument : bool
        If True, a click.argument is returned.
        If False, a click.option is returned.
        Default is False.
    **kwargs
        Items that overwrite the default values of
        lib/interface_parameters.json.
        These need to be valid items for click.
        A special "type" dict can be used to pass a click.Path or click.Choice,
        that has the following format:
        type = {"name": "path" or "choice", **choice_or_path_kwargs}
    Returns
    -------
    : click.option, click.argument
        A click.option or click.argument decorator.
    """
    parameters = copy.deepcopy(
        INTERFACE_PARAMETERS[parameter_name]
    )
    parameters.update(kwargs)
    if "type" in parameters:
        if parameters["type"] == "int":
            parameters["type"] = int
        elif parameters["type"] == "float":
            parameters["type"] = float
        elif parameters["type"] == "str":
            parameters["type"] = str
        elif isinstance(parameters["type"], dict):
            parameter_type = parameters["type"].pop("name")
            if parameter_type == "path":
                parameters["type"] = click.Path(**parameters["type"])
                if ("default" in parameters) and (parameters["default"]):
                    parameters["default"] = os.path.join(
                        os.path.dirname(__file__),
                        parameters["default"]
                    )
            elif parameter_type == "choice":
                options = parameters["type"].pop("options")
                parameters["type"] = click.Choice(
                    options,
                    **parameters["type"]
                )
    if "default" in parameters:
        if "is_flag" in parameters:
            parameters["show_default"] = False
        else:
            parameters["show_default"] = True
    if not as_argument:
        if "short_name" in parameters:
            short_name = parameters.pop("short_name")
            return click.option(
                f"--{parameter_name}",
                f"-{short_name}",
                **parameters,
            )
        else:
            return click.option(
                f"--{parameter_name}",
                **parameters,
            )
    else:
        if "required" in parameters:
            required = parameters.pop("required")
        else:
            required = True
        if "nargs" in parameters:
            return click.argument(
                parameter_name,
                type=click.Path(exists=True),
                nargs=parameters["nargs"],
                required=required,
            )
        else:
            return click.argument(
                parameter_name,
                type=parameters["type"],
                required=required,
            )




@click.group(
    context_settings=dict(
        help_option_names=['-h', '--help'],
    ),
    invoke_without_command=True
)
@click.pass_context
@click.version_option(directlfq.__version__, "-v", "--version")
def run(ctx, **kwargs):


    name = f"directLFQ {directlfq.__version__}"
    ascii_art_name = """
     _ _               _   _     ______ _____
    | (_)             | | | |    |  ___|  _  |
  __| |_ _ __ ___  ___| |_| |    | |_  | | | |
 / _` | | '__/ _ \/ __| __| |    |  _| | | | |
| (_| | | | |  __/ (__| |_| |____| |   \ \/' /
 \__,_|_|_|  \___|\___|\__\_____/\_|    \_/\_\


"""
    width = 47
    centered_name = " "*13 +f"* {name} *"

    click.echo("\n")
    click.echo("*" * (width))
    click.echo(ascii_art_name)
    click.echo(centered_name)
    click.echo("*" * (width))
    click.echo("\n")
    if ctx.invoked_subcommand is None:
        click.echo(run.get_help(ctx))

@run.command("gui", help="Start graphical user interface.")
@click.option("--port", "-p", type=int, default=None, help="Port to run the GUI server on (default: 41215 or PORT environment variable)")
def gui(port):
    import directlfq.gui
    directlfq.gui.run(port=port)

list_of_format_names = ["alphapept_peptides","fragpipe_precursors","maxquant_evidence","maxquant_peptides","diann_fragion_isotopes","diann_precursors","spectronaut_fragion_isotopes","spectronaut_precursor"]

@run.command("lfq", help="Run directLFQ normalization on proteomics input table.", no_args_is_help=True)
@click.option("--input_file", "-i", type=click.Path(exists=True), required=True, help="The input file containing the ion intensities. Usually the output of a search engine.")
@click.option("--columns_to_add", "-ca", type=list, default=[], multiple=True, help="Additional columns to add to the LFQ intensity output table. They are extraced from the input file.")
@click.option("--selected_proteins_file", "-sp", type=click.Path(exists=True),
default=None, help="If you want to perform normalization only on a subset of proteins, you can pass a .txt file containing the protein IDs, separeted by line breaks. No header expected.")
@click.option("--mq_protein_groups_txt", "-mp", type=click.Path(exists=True), default=None,
help="In the case of using MaxQuant data, the proteinGroups.txt table is needed in order to map IDs analogous to MaxQuant. Adding this table improves protein mapping, but is not necessary.")
@click.option("--min_nonan", "-mn", type=int, default=1, help="Min number of ion intensities necessary in order to derive a protein intensity. Increasing the number results in more reliable protein quantification at the cost of losing IDs.")
@click.option("--input_type_to_use", "-it", type=click.Choice(list_of_format_names), default=None, help="The type of input file to use. This is used to determine the column names of the input file. Only change this if you want to use non-default settings.")
@click.option("--maximum_number_of_quadratic_ions_to_use_per_protein", "-mn", type= int, default = 10,  help="How many ions are used to create the anchor intensity trace (see paper). Increasing might marginally increase performance at the cost of runtime.")
@click.option("--number_of_quadratic_samples", "-nq", type = int, default = 50, help="How many samples are used to create the anchor intensity trace (see paper). Increasing might marginally increase performance at the cost of runtime.")
@click.option("--filename_suffix", "-fs", type=str, default="", help="A suffix to add to the output file name.")
@click.option("--num_cores",  "-nc", type = int, default = None, help="The number of cores to use (default is to use multiprocessing).")
@click.option("--deactivate_normalization",  "-dn", type = bool, default = False, help="If you want to deactivate the normalization step, you can set this flag to True.")
@click.option("--filter_dict",  "-fd", type = str, default = None, help="In case you want to define specific filters in addition to the standard filters, you can add a yaml file where the filters are defined (see GitHub docu for example).")

def run_directlfq(**kwargs):
    print("starting directLFQ")
    import directlfq.lfq_manager
    directlfq.lfq_manager.run_lfq(**kwargs)
