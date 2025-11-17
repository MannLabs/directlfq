
<!-- ![Pip installation](https://github.com/MannLabs/directlfq/workflows/Default%20installation%20and%20tests/badge.svg) -->
![Unit tests](https://github.com/MannLabs/directlfq/actions/workflows/nbdev_tests.yml/badge.svg)
![System tests](https://github.com/MannLabs/directlfq/actions/workflows/quick_tests.yml/badge.svg)

<!-- [![Downloads](https://pepy.tech/badge/directlfq)](https://pepy.tech/project/directlfq)
[![Downloads](https://pepy.tech/badge/directlfq/month)](https://pepy.tech/project/directlfq)
[![Downloads](https://pepy.tech/badge/directlfq/week)](https://pepy.tech/project/directlfq) -->



# directLFQ
directLFQ is an open-source Python package for quantifying protein intensities based on peptide intensities or fragment-ion intensities measured with from Mass Spectrometry-based proteomics. It preserves peptide ratios, shows very accurate quantification and has a robust normalization approach. Furthermore, it allows fast processing also of very large sample cohorts, as runtime increases linearly with sample number. It is part of the AlphaPept ecosystem from the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and the [University of Copenhagen](https://www.cpr.ku.dk/research/proteomics/mann/).

You can process DIA and DDA data analyzed by [AlphaPept](https://github.com/MannLabs/alphapept), MaxQuant, FragPipe, Spectronaut and DIANN as well as [generic formats](#generic-input-format), using a Graphical User Interface (GUI) or the python package.


- [directLFQ](#directlfq)
  - [About](#about)
  - [Installation](#installation)
    - [One-click GUI](#one-click-gui)
    - [Pip](#pip)
    - [Developer](#developer)
  - [Running directLFQ](#running-directlfq)
    - [GUI](#gui)
    - [CLI](#cli)
    - [Python and Jupyter notebooks](#python-and-jupyter-notebooks)
  - [Troubleshooting](#troubleshooting)
  - [Citations](#citations)
  - [How to contribute](#how-to-contribute)
  - [License](#license)
  - [directLFQ commands](#directlfq-commands)
  - [directLFQ output](#directlfq-output)
  - [preparing input files](#preparing-input-files)
    - [Spectronaut](#spectronaut)
    - [DIA-NN](#dia-nn)
    - [MaxQuant](#maxquant)
    - [FragPipe](#fragpipe)
    - [generic input format](#generic-input-format)
  - [reproducing data analyses from the paper](#reproducing-data-analyses-from-the-paper)

---
## About
Generating protein intensities from Mass Spectrometry proteomics data comes with a variety of challenges. Differing peptides that belong to the same protein can have strongly differing intensities, for example due to differing ionization efficiencies. Missing values (i.e. peptides that have been detected in one run but not in the other) make simple summarization of peptide intensities to protein intensities problematic. Differences in sample loading can introduce systematic biases into the analysis.
With directLFQ, we provide a novel algorithm for addressing these challenges in an efficient and accurate manner. directLFQ retains peptide ratios and uses them to infer protein ratios and uses the concept of *intensity traces* for it's main processing steps. For futher details on the algorithm, please refer to the [preprint](https://doi.org/10.1101/2023.02.17.528962).

![](./release/images/methods_figure.png)

---
## Installation

directLFQ can be installed and used on all major operating systems (Windows, macOS and Linux).
There are currently two different types of installation possible:

* [**One-click GUI installation:**](#one-click-gui-installation) Choose this installation if you only want the GUI and/or keep things as simple as possible.
* [**Pip installation:**](#pip-installation) Choose this installation if you want to use directlfq as a Python package in an existing Python 3.9 environment (e.g. a Jupyter notebook). If needed, the GUI and CLI can be installed with pip as well.
* [**Developer installation:**](#developer-installation) Choose this installation if you are familiar with CLI tools, [conda](https://docs.conda.io/en/latest/) and Python. This installation allows access to all available features of directlfq and even allows to modify its source code directly. Generally, the developer version of directlfq outperforms the precompiled versions which makes this the installation of choice for high-throughput experiments.
- [**Docker installation:**](#docker-installation) Choose this installation if you want to use directlfq without any changes to your system.

### One-click GUI installation

The GUI of directlfq is a completely stand-alone tool that requires no knowledge of Python or CLI tools. Click [**here**](https://github.com/MannLabs/directlfq/releases/latest/) to download the latest release for either Windows, macOS or Linux.

Older releases remain available on the [release page](https://github.com/MannLabs/directlfq/releases), but no backwards compatibility is guaranteed.

### Pip installation

directLFQ can be installed in an existing Python 3.9 environment with a single `bash` command.

```bash
pip install "directlfq[dask]"
```

This installs the core directLFQ without graphical user interface (GUI). If you want to install with additional dependencies for GUI support, you can do this with:

```bash
pip install "directlfq[gui,dask]"
```

For installation with stable dependencies, use:

```bash
pip install "directlfq[stable,gui-stable,dask-stable]"
```

NOTE: You might need to run `pip install pip==21.0` before installing directlfq like this. Also note the double quotes `"`.
Note: The 'dask' extra can be omitted for slimmer dependencies, but this comes with a performance penalty for large datasets.

For those who are really adventurous, it is also possible to directly install any branch (e.g. `@development`) with any extras (e.g. `#egg=directlfq[stable,development]`) from GitHub with e.g.

```bash
pip install "git+https://github.com/MannLabs/directlfq.git@development#egg=directlfq[stable,development]"
```

### Developer installation

directlfq can also be installed in editable (i.e. developer) mode with a few `bash` commands. This allows to fully customize the software and even modify the source code to your specific needs. When an editable Python package is installed, its source code is stored in a transparent location of your choice. While optional, it is advised to first (create and) navigate to e.g. a general software folder:

```bash
mkdir ~/folder/where/to/install/software
cd ~/folder/where/to/install/software
```

***The following commands assume you do not perform any additional `cd` commands anymore***.

Next, download the directlfq repository from GitHub either directly or with a `git` command. This creates a new directlfq subfolder in your current directory.

```bash
git clone https://github.com/MannLabs/directlfq.git
```

For any Python package, it is highly recommended to use a separate [conda virtual environment](https://docs.conda.io/en/latest/), as otherwise *dependancy conflicts can occur with already existing packages*.

```bash
conda create --name directlfq python=3.9 -y
conda activate directlfq
```

Finally, directlfq and all its [dependencies](requirements) need to be installed. To take advantage of all features and allow development (with the `-e` flag), this is best done by also installing the [development dependencies](requirements/requirements_development.txt) instead of only the [core dependencies](requirements/requirements.txt):

```bash
pip install -e "./directlfq[development,gui]"
```

By default this installs loose dependencies (no explicit versioning), although it is also possible to use stable dependencies (e.g. `pip install -e "./directlfq[stable,development]"`).

***By using the editable flag `-e`, all modifications to the [directlfq source code folder](directlfq) are directly reflected when running directlfq. Note that the directlfq folder cannot be moved and/or renamed if an editable version is installed. In case of confusion, you can always retrieve the location of any Python module with e.g. the command `import module` followed by `module.__file__`.***


### Docker installation
The containerized version can be used to run directlfq without any installation to your system.

#### 1. Setting up Docker
Install the latest version of docker (https://docs.docker.com/engine/install/).

#### 2. Prepare folder structure
Set up your data to match the expected folder structure:
create a folder and store its name in a variable, and specify a port
```
DATA_FOLDER=/home/username/data; mkdir -p $DATA_FOLDER
PORT=41215
```

#### 3. Start the container
```bash
docker run -v $DATA_FOLDER:/app/data -p $PORT:41215 mannlabs/directlfq:latest
```
After initial download of the container, directlfq will start running immediately,
and can be accessed under [localhost:$PORT](http://localhost:41215).

Note: in the app, the local `$DATA_FOLDER` needs to be referred to as "`/app/data`".

#### Alternatively: Build the image yourself
If you want to build the image yourself, you can do so by
```bash
docker build -t directlfq .
```
and run it with
```bash
docker run -p $PORT:8501 -v $DATA_FOLDER:/app/data -t directlfq
```

---

## Running directLFQ

There are three ways to use directlfq:

* [**GUI**](#gui)
* [**CLI**](#cli)
* [**Python**](#python-and-jupyter-notebooks)

NOTE: The first time you use a fresh installation of directlfq, it is often quite slow because some functions might still need compilation on your local operating system and architecture. Subsequent use should be a lot faster.

### GUI
![](./release/images/gui_screenshot.png)
If you have installed directlfq with the one-click GUI installer, you can run the GUI simply by clicking the directLFQ icon on your desktop/applications folder.

If the GUI was not installed through a one-click GUI installer, it can be activate with the following `bash` command:

```bash
directlfq gui
```

Note that this needs to be prepended with a `!` when you want to run this from within a Jupyter notebook. When the command is run directly from the command-line, make sure you use the right environment (activate it with e.g. `conda activate directlfq` or set an alias to the binary executable (can be obtained with `where directlfq` or `which directlfq`)).


### CLI
<img src="./release/images/cli_screenshot.png" width="600" />


The CLI can be run with the following command (after activating the `conda` environment with `conda activate directlfq` or if an alias was set to the directlfq executable):

```bash
directlfq -h
```

It is possible to get help about each function and their (required) parameters by using the `-h` flag.


### Python and Jupyter notebooks

directLFQ can be imported as a Python package into any Python script or notebook with the command `import directlfq`.
Running the standard analysis (with plots) can be done via the command:
```bash
import directlfq.lfq_manager as lfq_manager

example_input_file_diann = "/path/to/example_input_file_diann.tsv"

lfq_manager.run_lfq(example_input_file_diann)
```

Several use cases for applying directLFQ can be found as Jupyter Notebooks in the tests folder. See for example the [quicktests notebook](tests/quicktests/run_pipeline_w_different_input_formats.ipynb).

Note that the nbdev_nbs folder contains the source code as Jupyter notebooks. These notebooks are automatically converted to Python scripts using the nbdev package and stored in the directlfq folder. The notebooks contain additional documentation and comments as well as unit tests that can be executed directly from the notebooks themselves.

---
## Troubleshooting

In case of issues, check out the following:

* [Issues](https://github.com/MannLabs/directlfq/issues): Try a few different search terms to find out if a similar problem has been encountered before
* [Discussions](https://github.com/MannLabs/directlfq/discussions): Check if your problem or feature requests has been discussed before.

---
## Citations

In the case that directLFQ is useful to you, please consider supporting us by citing the [paper](https://doi.org/10.1016/j.mcpro.2023.100581)

Ammar, C., Schessner, J.P., Willems, S., Michaelis, A.C., and Mann, M. (2023). Accurate label-free quantification by directLFQ to compare unlimited numbers of proteomes.  Molecular & Cellular Proteomics, 100581.



---
## How to contribute

If you like this software, you can give us a [star](https://github.com/MannLabs/directlfq/stargazers) to boost our visibility! All direct contributions are also welcome. 
Feel free to post a new [issue](https://github.com/MannLabs/directlfq/issues) or clone the repository and create a [pull request](https://github.com/MannLabs/directlfq/pulls) with a new branch. 
For an even more interactive participation, check out the [discussions](https://github.com/MannLabs/directlfq/discussions) and the [the Contributors License Agreement](misc/CLA.md).

### Tagging of changes
In order to have release notes automatically generated, changes need to be tagged with labels.
The following labels are used (should be safe-explanatory):
`breaking-change`, `bug`, `enhancement`.

### Release a new version
This package uses a shared release process defined in the
[alphashared](https://github.com/MannLabs/alphashared) repository. Please see the instructions
[there](https://github.com/MannLabs/alphashared/blob/reusable-release-workflow/.github/workflows/README.md#release-a-new-version)

---

## License

directLFQ was developed by the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and the [University of Copenhagen](https://www.cpr.ku.dk/research/proteomics/mann/) and is freely available with an [Apache License](LICENSE.txt). External Python packages (available in the [requirements](requirements) folder) have their own licenses, which can be consulted on their respective websites.

---

## directLFQ commands
directLFQ is started internally via the directlfq.lfq_manager.run_lfq() command. In principle and for most use cases, you only need to provide the path to the AlphaPept/MaxQuant/DIA-NN etc. file of interest. However, there are several other options that can be used to customize the analysis. The following commands are available:

- **input_file**: The input file containing the ion intensities. Usually the output of a search engine.
- **columns_to_add**: Add the names of columns that are present in the output table and that you want to keep in the directLFQ output file. Separated by semicolons. Note that some basic additional columns such as gene names are always added to the output table by default. WARNING: Take care that columns you add are not ambigous. For example, adding the peptide sequence column will not work, because there are multiple peptide sequences per protein.
- **selected_proteins_file**: If you want to perform normalization only on a subset of proteins, you can pass a .txt file containing the protein IDs, separeted by line breaks. No header expected.
- **mq_protein_groups_txt**: In the case of using MaxQuant data, the proteinGroups.txt table is needed in order to map IDs analogous to MaxQuant. Adding this table improves protein mapping, but is not necessary.
- **min_nonan**: Min number of ion intensities necessary in order to derive a protein intensity. Increasing the number results in more reliable protein quantification at the cost of losing IDs.
- **input_type_to_use**: The type of input file to use. This is used to determine the column names of the input file. Only change this if you want to use non-default settings.
- **maximum_number_of_quadratic_ions_to_use_per_protein**: How many ions are used to create the anchor intensity trace (see paper). Increasing might marginally increase performance at the cost of runtime.
- **number_of_quadratic_samples**: How many samples are used to create the anchor intensity trace (see paper). Increasing might marginally increase performance at the cost of runtime
-  **num_cores**: The number of cores to use (default is to use multiprocessing).
-  **filename_suffix**: Suffix to append to the output files.
-  **deactivate_normalization**: Set to true, if no between-sample normalization should be performed before processing.
-  **filter_dict**: In case you want to define specific filters in addition to the standard filters, you can add a yaml file where the filters are defined (see example [here](misc/examples/filterdict.yaml)). In the Python API you can also directly put in the dictionary instead of the .yaml file.

---

## directLFQ output
directLFQ writes three output files into the directory of the input file:
1.  The main output file ends with .protein_intensities.tsv and contains the estimated LFQ protein intensities. 
2.  The second output file ends with .ion_intensities.tsv and contains the aligned intensity traces of all ions. This allows to compare profiles of different ions to each other. In particular, if you run directLFQ with peptide-level quantification, you can use this file to compare the intensity traces of different peptides of the same protein.
3.  The third output file ends with .aq_reformat.tsv and contains the reformatted input data in matrix format (ions are rows, samples are columns). The values are identical to the values of the original input file, just the format is different.

---

## preparing input files

### Spectronaut
directLFQ takes a Spectronaut .tsv table as input. When exporting from Spectronaut, the correct columns need to be selected. These can be obtained by downloading one of the export schemes available below. We provide one export scheme for sprecursor quantification and one export scheme for fragment ion quantification. Fragment ion quantification shows slightly more accuracy, but the files are around 10 times larger.

An export scheme can then simply be loaded into Spectronaut as follows:

Go to the "Report" perspective in Spectronaut, click "Import Schema" and provide the file.

The data needs to be exported in the normal long format as .tsv file


<a href="https://github.com/MannLabs/directlfq/raw/development/directlfq/configs/spectronaut_tableconfig_precursor.rs" download>Download Spectronaut export scheme for precursor quantification</a>

<a href="https://github.com/MannLabs/directlfq/raw/development/directlfq/configs/spectronaut_tableconfig_fragion.rs" download>Download Spectronaut export scheme for fragment ion quantification</a>



### DIA-NN
Provide the path to the DIANN "report.tsv" output table.

### MaxQuant
Provide the path to the MaxQuant "peptides.txt" output table or the MaxQuant evidence.txt output table. Additionally and if possible, provide the path to the corresponding "proteinGroups.txt" file.

### FragPipe
Provide the path to the "combined_ion.tsv" output table.


### generic input format
In the case that you working with a search engine that is not supported by directLFQ, you can use the generic input format. This format is a tab-separated quantity matrix file with the following columns: "protein", "ion", "run_id1", "run_id2", ..,"run_idN". 
Each row contains therefore all the ion intensities that were measured for an ion in each run (see examples below). The ion identifier only needs to be unique for each ion and can be on the level you want (peptide, charged peptide, or fragment ion). After reformatting your file into this format, save the file with the ending ".aq_reformat.tsv". Then you can simply give this file as input to directLFQ and it will automatically detect the generic input format.
<img src="./release/images/example_input_format_mq.png" width="700" />
<img src="./release/images/example_input_format_sn.png" width="700" />

---

## reproducing data analyses from the paper
If you want to reproduce data analyses presented in the [manuscript](#citations), you can first download the data by executing
```bash
python tests/download_testfiles.py all_tests
```
This will download the underlying datasets into the appropriate location. The notebooks carrying out the analyses themselves are located in the tests directory in the respective subfolders ratio_tests, normalization_tests, runtime_tests and organellar_maps.
