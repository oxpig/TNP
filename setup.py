#!/usr/bin/env python
from distutils.core import setup

setup(
    name = "theraprofnano",
      version = "0.0.1",
      description = "TNP: Therapeutic Nanobody Profiler.",
      author = "Gemma Gordon",
      author_email = "gemma.gordon@wolfson.ox.ac.uk",
      license='BSD 3-clause license',               
      packages = [
          "theraprofnano",                                           # Entire module                 
          "theraprofnano.CDR_Profiler",                              # Analyses CDR properties                   
          "theraprofnano.Hydrophobicity_and_Charge_Profiler",        # Calculates hydrophobicity and charge metrics
          "theraprofnano.Hydrophobicity_and_Charge_Profiler.Common", # PDB atom/residue classes etc.
          "theraprofnano.Plotters", 
          "scripts"                                 # For plotting graphs for the web front.
      ],
      package_dir = {
          "theraprofnano": "lib/python/theraprofnano",
          "theraprofnano.CDR_Profiler": "lib/python/theraprofnano/CDR_Profiler",           
          "theraprofnano.Hydrophobicity_and_Charge_Profiler": "lib/python/theraprofnano/Hydrophobicity_and_Charge_Profiler",
          "theraprofnano.Hydrophobicity_and_Charge_Profiler.Common": "lib/python/theraprofnano/Hydrophobicity_and_Charge_Profiler/Common",
          "theraprofnano.Plotters": "lib/python/theraprofnano/Plotters"
      },
      package_data = {
          "theraprofnano.Hydrophobicity_and_Charge_Profiler": [ "dat/*" ],
      },
      scripts = ["bin/TNP"], #, "bin/psa", "bin/psa_mac"],
      data_files = [ ('bin', ['bin/psa', 'bin/psa_mac']) ] # set as executable with chmod +x bin/TNP
)
