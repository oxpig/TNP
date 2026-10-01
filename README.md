-----------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------
                \/                       _____ _   _ ____                   \/
               ⊂'l                      |_   _| \ | |  _ \                 ⊂'l     
                ll                        | | |  \| | |_) |                 ll     
                llama~                    | | | |\  |  __/                  llama~ 
                || ||                     |_| |_| \_|_|                     || || 
                '' ''                                                       '' ''
-----------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------

<div align="center">    
 
# The Therapeutic Nanobody Profiler: characterising and predicting nanobody developability to improve therapeutic design

</div>

## About

This software was developed in the _Oxford Protein Informatics Group_ ([OPIG](http://opig.stats.ox.ac.uk/)), Department of Statistics, University of Oxford with the support of Twist Bioscience.

**Authors**

* Gemma L. Gordon (University of Oxford)
* João Gervasio (University of Oxford, Okinawa Institute of Science and Technology Graduate University)
* Colby Souders (Twist Bioscience)
* Charlotte M. Deane (University of Oxford)


## Abstract 

Developability optimisation is an important step for successful biotherapeutic design. For monoclonal antibodies, developability is relatively well characterised. However, progress for novel biotherapeutics such as nanobodies is more limited. Differences in structural features between antibodies and nanobodies render current antibody computational methods unsuitable for direct application to nanobodies. Following the principles of the Therapeutic Antibody Profiler (TAP), we have built the Therapeutic Nanobody Profiler (TNP), an open-source computational tool for predicting nanobody developability. Tailored specifically for nanobodies, it accounts for their unique properties compared to conventional antibodies for more efficient development of this novel therapeutic format. We calibrate TNP metrics using the 36 currently available clinical-stage nanobody sequences. We also collected experimental developability data for 108 nanobodies and examine how these results are related to the TNP guidelines. TNP is available as a web application which you can find <a href="https://opig.stats.ox.ac.uk/webapps/tnp">on the OPIG website.</a>

## Citing this work

The code and data in this package is based on the <a href="https://doi.org/10.1101/2025.08.11.669635">following paper.</a> If you use it, please cite:

```tex
@article{gordon2025,
  title={The Therapeutic Nanobody Profiler: characterizing and predicting nanobody developability to improve therapeutic design},
  author={},
  journal={bioRxiv},
  pages={},
  year={2025},
  publisher={Cold Spring Harbor Laboratory},
  doi={https://doi.org/10.1101/2025.08.11.669635}
}
```

## Web tool

This tool is also available as a live [web app](https://opig.stats.ox.ac.uk/webapps/tnp):

<https://opig.stats.ox.ac.uk/webapps/tnp>

If you use this code or the web tool in your work, please [cite the paper](#citing-this-work).

## Installation

TNP is distributed through [Bioconda](https://bioconda.github.io/):

```bash
conda install -c conda-forge -c bioconda tnp
```

TNP depends on [ANARCI](https://github.com/oxpig/ANARCI) and HMMER from
Bioconda, and on [NanoBodyBuilder2](https://github.com/oxpig/ImmuneBuilder),
DSSP and FreeSASA.  Do not install the `anarci` package from PyPI, which is not
an official release.

NanoBodyBuilder2 downloads its trained weights on first use.  Conda installs a
CUDA-enabled PyTorch where it detects a GPU, and a CPU-only one otherwise.

## Usage

For a single sequence:

```bash
TNP --name my_sequence --output my_sequence_output --seq [sequence]
```

For multiple sequences in a FASTA file:

```bash
TNP --output /path/to/output/directory --file /path/to/sequences.fasta
```

TNP models each sequence with NanoBodyBuilder2, then writes its results to the
output directory: a log and a JSON summary of the metrics and their flags, and
the models in `Final_Models/`.  With `--web`, it also writes the data for the
plots shown by the web application.

For more options and information, run `TNP --help`.

## Development

To set up a development environment with Conda (or Mamba/Micromamba), in the
repository:

```bash
conda env create -f environment.yml
conda activate tnp
pip install --no-deps --editable .
```

Run the tests with `pytest`.  The slow end-to-end tests, which model a nanobody
with NanoBodyBuilder2, are skipped unless you run `pytest --runslow`.  Lint with
`ruff check`, and install the pre-commit hooks with `uvx pre-commit install`.

`scripts/compare_sasa.py` compares the surface patch metrics on the clinical-
stage models in `paper/paper_data` with the published values, or with another
version of TNP.

### Releasing

Bump the version, which commits the change and makes a signed tag, then push
it:

```bash
uvx bump-my-version bump patch  # Or minor, or major.
git push --follow-tags
```

For each version tag, GitHub Actions runs the tests and makes a GitHub release
with a signed source distribution.  The Bioconda bot then updates the recipe.
