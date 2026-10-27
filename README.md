# Modeling Impact of Ser/Thr Phosphorylation on Transcription Factors in _Mtb_

<!--toc:start-->

- [Modeling Impact of Ser/Thr Phosphorylation on Transcription Factors in _Mtb_](#modeling-impact-of-serthr-phosphorylation-on-transcription-factors-in-mtb)
  - [Directory Structure](#directory-structure)
  - [Usage](#usage)
    - [Downloading External Data](#downloading-external-data)
    - [Installing Dependencies](#installing-dependencies)
      - [Conda Environment](#conda-environment)
      - [Pip](#pip)
      - [Pixi](#pixi)
      - [UV](#uv)
    - [Results](#results)
  - [Licensing](#licensing)
  - [References](#references)
    - [Software](#software)
    - [Data](#data)

<!--toc:end-->

Code associated with "Regulation of _Mycobacterium tuberculosis_ transcriptional
repressors by Ser/Thr phosphorylation"

## Directory Structure

<details>
<summary>Descriptions of the files in this repository</summary>
<details>
<summary>data/</summary>

- Mycobacterium_tuberculosis_H37Rv_txt_v5.txt: Gene information for _Mtb_ from
  [Mycobrowser](https://mycobrowser.epfl.ch/), see [reference](#references)
  below.
- bosch_vulnerability_index.xlsx: Gene essentiality and vulnerability
  information from (Bosch et al., 2021)
- tfoe_targets.xlsx: Targets of the _Mtb_ transcription factors from (Rustad et
  al., 2014)

</details>
<details>
<summary>results/</summary>

- phosphosite_linear_regression_results.csv: All transcription factor
  phosphosite predictions
- visualization: Directory of heatmap visualizations of the transcription factor
  phosphosite activity predictions, named according to the gene name of the
  transcription factor (or locus tag if no gene name has been assigned)

</details>
<details>
<summary>scripts/</summary>

- 01_kinase_tf_phosphosite_regression.py: Script for fitting the linear models
- 02_kinase_tf_phosphosite_visualization.py: Script for generating the heatmap
  visualizations for the transcription factor phosphosite activity predictions
- common_functions.py: Some utility functions used by the scripts

</details>

- environment.yml: Conda environment specification
- LICENSE: License for reusing the code in this repository (does NOT apply to
  external data in data directory, see [licensing](#licensing) for details)
- README.md: This file
- pixi.lock: Lock file for pixi project/package manager
- pixi.toml: Configuration for pixi project/package manager
- pylock.toml: Python lockfile
- requirements.txt: Python requirements file
- ruff.toml: Configuration for ruff linter
- uv.lock: UV lockfile

</details>

## Usage

### Downloading External Data

Some data couldn't be included in this repository, and so needs to be downloaded
in order for the analysis to be performed. Specifically, data from 'Frando, A.,
Boradia, V., Gritsenko, M., Beltejar, C., Day, L., Sherman, D. R., Ma, S.,
Jacobs, J. M., & Grundner, C. (2023). The Mycobacterium tuberculosis protein
O-phosphorylation landscape. Nature Microbiology, 8(3), 548–561.'.

These data need to be downloaded, and saved into the 'data' directory of this
repository **BEFORE** the analysis can be performed.

Download:

- 'Supplementary Table 2' from either
  [Nature Microbiology](https://www.nature.com/articles/s41564-022-01313-7#Sec26)
  or from [Pubmed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC11376436/).
- 'Supplementary Table 3' from either
  [Nature Microbiology](https://www.nature.com/articles/s41564-022-01313-7#Sec26)
  or from [Pubmed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC11376436/).
- 'GSE195959_Expression.MTB.GeneData.csv.gz' from
  [Gene Expression Omnibus GSE195959](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195959)

And place both of these files into the data directory (no need to unzip the
gzipped file, but you can if desired).

### Installing Dependencies

This project mainly uses [pixi](https://pixi.prefix.dev/latest/) to manage
dependencies, but includes several different alternatives for installing
dependencies, including [conda](#conda-environment), [pip](#pip), [uv](#uv), and
[pixi](#pixi).

#### Conda Environment

Conda can be installed using [Anaconda](https://www.anaconda.com/download),
[Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main), or
[Miniforge](https://github.com/conda-forge/miniforge). We recommend using
Miniforge as that is what this repository was tested with, but any method should
work.

After installing conda using one of the above methods, then the environment can
be created using the environment.yml file:

```{bash}
conda create --file environment.yml
```

and then activated with

```{bash}
conda activate stpk_tf_regression
```

The scripts can then be run by moving to the scripts directory, and running them
in the order of their prefixed number:

```{bash}
cd scripts 
python 01_kinase_tf_phosphosite_regression.py
python 02_kinase_tf_phosphosite_visualization.py
```

#### Pip

Make sure python is installed, you could use:

- [python download page](https://www.python.org/downloads/)
- [pyenv](https://github.com/pyenv/pyenv)
- [mise-en-place](https://mise.jdx.dev/lang/python.html)
- [brew](https://brew.sh/)
- [scoop](https://scoop.sh/)

or other package manager to install python.

Once python is installed, a virtual environment can be created by

```{bash}
python -m venv .venv
```

then activated with

```{bash}
# Windows
.venv\Scripts\activate
# Macos/Linux
source myfirstproject/bin/activate
```

and finally required dependencies can be installed with

```{bash}
pip install -r requirements.txt 
```

The scripts can then be run by moving to the scripts directory, and running them
in the order of their prefixed number:

```{bash}
cd scripts 
python 01_kinase_tf_phosphosite_regression.py
python 02_kinase_tf_phosphosite_visualization.py

```

#### Pixi

Pixi can be installed based on instructions from
[https://pixi.prefix.dev/latest/installation/](https://pixi.prefix.dev/latest/installation/).

The scripts can then be run through pixi tasks

```{bash}
pixi run regression # For running 01_kinase_tf_phosphosite_regression.py
pixi run visualization # For running 02_kinase_tf_phosphosite_visualization.py
pixi run all # For running the full analysis
```

#### UV

UV can be installed based on instructions from
[https://docs.astral.sh/uv/getting-started/installation/](https://docs.astral.sh/uv/getting-started/installation/).
Then the dependencies can be installed with

```{bash}
uv sync
```

The environment can then be activated with

```{bash}
# Windows
.venv\Scripts\activate
# Macos/Linux
source myfirstproject/bin/activate
```

The scripts can then be run by moving to the scripts directory, and running them
in the order of their prefixed number:

```{bash}
cd scripts 
python 01_kinase_tf_phosphosite_regression.py
python 02_kinase_tf_phosphosite_visualization.py
```

### Results

Results of the analysis can be found in the 'results/' directory

## Licensing

The code in this repository, specifically in the scripts directory, is licensed
under an MIT license, see the LICENSE file for more information.

In the data directory

- bosch_vi.xlsx: From
  [Genome-wide gene expression tuning reveals diverse vulnerabilities of M. tuberculosis](https://doi.org/10.1016/j.cell.2021.06.033),
  licensed under the CC-BY license
- tfoe_targets.xlsx: From
  [Mapping and manipulating the Mycobacterium tuberculosis transcriptome using a transcription factor overexpression-derived regulatory network](https://doi.org/10.1186/s13059-014-0502-3),
  licensed under the
  [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/)
- Mycobacterium_tuberculosis_H37Rv_txt_v5.txt: Release from
  [Mycobrowser](https://mycobrowser.epfl.ch/), reference in
  [References](#references) below

## References

### Software

- Matplotlib: Hunter, J. D. (2007). Matplotlib: A 2D graphics environment.
  Computing in Science & Engineering, 9(3), 90–95.
  [https://doi.org/10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55)
- NumPy: Harris, C. R., Millman, K. J., van der Walt, S. J., Gommers, R.,
  Virtanen, P., Cournapeau, D., Wieser, E., Taylor, J., Berg, S., Smith, N. J.,
  Kern, R., Picus, M., Hoyer, S., van Kerkwijk, M. H., Brett, M., Haldane, A.,
  del Río, J. F., Wiebe, M., Peterson, P., … Oliphant, T. E. (2020). Array
  programming with NumPy. Nature, 585(7825), 357–362.
  [https://doi.org/10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)
- pandas:
  - Software: team, T. pandas development. (2020). pandas-dev/pandas: Pandas
    (Version latest) [Computer software]. Zenodo.
    [https://doi.org/10.5281/zenodo.3509134](https://doi.org/10.5281/zenodo.3509134)
  - Paper: McKinney, W. (2010). Data Structures for Statistical Computing in
    Python. In S. van der Walt & J. Millman (Eds.), Proceedings of the 9th
    Python in Science Conference (pp. 56–61).
    [https://doi.org/10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)
- scikit-learn: Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion,
  B., Grisel, O., Blondel, M., Prettenhofer, P., Weiss, R., Dubourg, V.,
  Vanderplas, J., Passos, A., Cournapeau, D., Brucher, M., Perrot, M., &
  Duchesnay, E. (2011). Scikit-learn: Machine Learning in Python. Journal of
  Machine Learning Research, 12, 2825–2830.
- SciPy: Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler
  Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser,
  Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K.
  Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern,
  Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas,
  Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero,
  Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa,
  Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental
  Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.
  DOI: [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2).
- seaborn: Waskom, M. L. (2021). seaborn: Statistical data visualization.
  Journal of Open Source Software, 6(60), 3021.
  [https://doi.org/10.21105/joss.03021](https://doi.org/10.21105/joss.03021)
- statsmodels: Seabold, S., & Perktold, J. (2010). statsmodels: Econometric and
  statistical modeling with python. 9th Python in Science Conference.

### Data

- bosch_vi.xlsx: Data from 'Bosch B, DeJesus MA, Poulton NC, Zhang W, Engelhart
  CA, Zaveri A, Lavalette S, Ruecker N, Trujillo C, Wallach JB, Li S, Ehrt S,
  Chait BT, Schnappinger D, Rock JM. Genome-wide gene expression tuning reveals
  diverse vulnerabilities of M. tuberculosis. Cell. 2021 Aug
  19;184(17):4579-4592.e24. doi: 10.1016/j.cell.2021.06.033. Epub 2021 Jul 22.
  [PMID: 34297925](https://pubmed.ncbi.nlm.nih.gov/34297925/);
  [PMCID: PMC8382161](https://pmc.ncbi.nlm.nih.gov/articles/PMC8382161/).'
- Phosphorylation and Gene expression data from STPK mutant strains: Data from
  'Frando, A., Boradia, V., Gritsenko, M., Beltejar, C., Day, L., Sherman, D.
  R., Ma, S., Jacobs, J. M., & Grundner, C. (2023). The Mycobacterium
  tuberculosis protein O-phosphorylation landscape. Nature Microbiology, 8(3),
  548–561.
- tfoe_targets.xlsx: Data from 'Rustad, T.R., Minch, K.J., Ma, S. et al. Mapping
  and manipulating the Mycobacterium tuberculosis transcriptome using a
  transcription factor overexpression-derived regulatory network. Genome Biol
  15, 502 (2014).
  [https://doi.org/10.1186/s13059-014-0502-3](https://doi.org/10.1186/s13059-014-0502-3)'
  [PMID: 25380655](https://pmc.ncbi.nlm.nih.gov/articles/PMC4249609/)
- Mycobacterium_tuberculosis_H37Rv_txt_v5.txt: MycoBrowser release V5 from
  'Kapopoulou A, Lew JM, Cole ST. The MycoBrowser portal: a comprehensive and
  manually annotated resource for mycobacterial genomes. Tuberculosis (Edinb).
  Jan 91(1):8–13 (2011)' at
  [https://mycobrowser.epfl.ch/](https://mycobrowser.epfl.ch/)
