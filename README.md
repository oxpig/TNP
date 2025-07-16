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
 
# The Therapeutic Nanobody Profiler: characterizing and predicting nanobody developability to improve therapeutic design

</div>

## About

This software was developed in the _Oxford Protein Informatics Group_ ([OPIG](http://opig.stats.ox.ac.uk/)), Department of Statistics, University of Oxford with the support of Twist Bioscience.

**Authors**

* Gemma L. Gordon (Oxford), 
* Charlotte M. Deane (Oxford)

## Citing this work

The code and data in this package is based on the <a href="">following paper.</a> If you use it, please cite:

[ADD IN PREPRINT]

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

- [ANARCI](https://github.com/oxpig/ANARCI)
- [SAbDab](https://github.com/oxpig/SAbDab)
- [NanoBodyBuilder2](https://github.com/oxpig/ImmuneBuilder)

Set PSA files as executable with chmod +x bin/TNP

## Usage

* For a single sequence

`TNP --name my_sequence --output my_sequence_output --seq [sequence]`

* For multiple sequences in a FASTA file

`TNP --name my_fasta_file --output /path/to/output/directory --file /path/to/my/fasta/file/filename.fasta`

* For a folder with already available PDB models

`TNP --name my_models --output /path/to/output/directory --models /path/to/my/models`

For more options and information, run `TNP --help`.


