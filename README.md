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

The code and data in this package is based on the <a href="">following paper.</a> If you use it, please cite:

```tex
@article{gordon2025,
  title={The Therapeutic Nanobody Profiler: characterizing and predicting nanobody developability to improve therapeutic design},
  author={},
  journal={bioRxiv},
  pages={},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Installation

In the package directory:

`pip install .`

**Requirements**

- python 3.10
- [ANARCI](https://github.com/oxpig/ANARCI)
- [NanoBodyBuilder2](https://github.com/oxpig/ImmuneBuilder)
- [DSSP](https://anaconda.org/salilab/dssp)
- biopython 1.77


## Usage

* For a single sequence

`TNP --name my_sequence --output my_sequence_output --seq [sequence]`

* For multiple sequences in a FASTA file

`TNP --name my_fasta_file --output /path/to/output/directory --file /path/to/my/fasta/file/filename.fasta`

<!-- * For a folder with already available PDB models: 
  
  - These structures should NOT contain hydrogens
  - SAbDab needs to be installed

`TNP --name my_models --output /path/to/output/directory --models /path/to/my/models` -->


For more options and information, run `TNP --help`.


