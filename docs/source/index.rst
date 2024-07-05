.. |doi| image:: https://zenodo.org/badge/379064819.svg
   :target: https://zenodo.org/badge/latestdoi/379064819

Welcome to ASCENT's documentation!
==================================

This documentation is an adaptation and update of the supplements associated with the original ASCENT publication.

**Please check out the associated** `publication <https://doi.org/10.1371/journal.pcbi.1009285>`_ **in PLOS Computational Biology!**

**Cite both the paper and the DOI for the release of the repository used for your work. We encourage you to clone the most recent commit of the repository.**

* **Cite the paper:**

  .. details:: APA
     :open:

     **Musselman, E. D.**, **Cariello, J. E.**, Grill, W. M., & Pelot, N. A. (2021). ASCENT (Automated Simulations to Characterize Electrical Nerve Thresholds): A pipeline for sample-specific computational modeling of electrical stimulation of peripheral nerves. PLOS Computational Biology, 17(9), e1009285. https://doi.org/10.1371/journal.pcbi.1009285
     
     **Peña, E.**, Pelot, N. A., & Grill, W. M. (2024). Computational models of compound nerve action potentials: Efficient filter-based methods to quantify effects of tissue conductivities, conduction distance, and nerve fiber parameters. PLoS computational biology, 20(3), e1011833. https://doi.org/10.1371/journal.pcbi.1011833
  
  .. details:: MLA

      Musselman, Eric D., et al. "ASCENT (Automated Simulations to Characterize Electrical Nerve Thresholds): A Pipeline for Sample-Specific Computational Modeling of Electrical Stimulation of Peripheral Nerves." PLOS Computational Biology, vol. 17, no. 9, Sept. 2021, p. e1009285. PLoS Journals, https://doi.org/10.1371/journal.pcbi.1009285.
      
      Peña, Edgar et al. “Computational models of compound nerve action potentials: Efficient filter-based methods to quantify effects of tissue conductivities, conduction distance, and nerve fiber parameters.” PLoS computational biology vol. 20,3 e1011833. 1 Mar. 2024, doi:10.1371/journal.pcbi.1011833
  
  .. details:: BibTeX

    .. code-block:: BibTeX

         @article{Musselman2021,
          doi = {10.1371/journal.pcbi.1009285},
          url = {https://doi.org/10.1371/journal.pcbi.1009285},
          year = {2021},
          month = sep,
          publisher = {Public Library of Science ({PLoS})},
          volume = {17},
          number = {9},
          pages = {e1009285},
          author = {Eric D. Musselman and Jake E. Cariello and Warren M. Grill and Nicole A. Pelot},
          editor = {Dina Schneidman-Duhovny},
          title = {{ASCENT} (Automated Simulations to Characterize Electrical Nerve Thresholds): A pipeline for sample-specific computational modeling of electrical stimulation of peripheral nerves},
          journal = {{PLOS} Computational Biology}
        }

         @article{Pena2024,
          doi = {10.1371/journal.pcbi.1011833},
          author = {Peña, Edgar AND Pelot, Nicole A. AND Grill, Warren M.},
          journal = {PLOS Computational Biology},
          publisher = {Public Library of Science},
          title = {Computational models of compound nerve action potentials: Efficient filter-based methods to quantify effects of tissue conductivities, conduction distance, and nerve fiber parameters},
          year = {2024},
          month = {03},
          volume = {20},
          url = {https://doi.org/10.1371/journal.pcbi.1011833},
          pages = {1-35},
          abstract = {Background Peripheral nerve recordings can enhance the efficacy of neurostimulation therapies by providing a feedback signal to adjust stimulation settings for greater efficacy or reduced side effects. Computational models can accelerate the development of interfaces with high signal-to-noise ratio and selective recording. However, validation and tuning of model outputs against in vivo recordings remains computationally prohibitive due to the large number of fibers in a nerve.   Methods We designed and implemented highly efficient modeling methods for simulating electrically evoked compound nerve action potential (CNAP) signals. The method simulated a subset of fiber diameters present in the nerve using NEURON, interpolated action potential templates across fiber diameters, and filtered the templates with a weighting function derived from fiber-specific conduction velocity and electromagnetic reciprocity outputs of a volume conductor model. We applied the methods to simulate CNAPs from rat cervical vagus nerve.   Results Brute force simulation of a rat vagal CNAP with all 1,759 myelinated and 13,283 unmyelinated fibers in NEURON required 286 and 15,860 CPU hours, respectively, while filtering interpolated templates required 30 and 38 seconds on a desktop computer while maintaining accuracy. Modeled CNAP amplitude could vary by over two orders of magnitude depending on tissue conductivities and cuff opening within experimentally relevant ranges. Conduction distance and fiber diameter distribution also strongly influenced the modeled CNAP amplitude, shape, and latency. Modeled and in vivo signals had comparable shape, amplitude, and latency for myelinated fibers but not for unmyelinated fibers.   Conclusions Highly efficient methods of modeling neural recordings quantified the large impact that tissue properties, conduction distance, and nerve fiber parameters have on CNAPs. These methods expand the computational accessibility of neural recording models, enable efficient model tuning for validation, and facilitate the design of novel recording interfaces for neurostimulation feedback and understanding physiological systems.},
          number = {3},
         }

* **Cite the code (use the DOI for the version of code used):** |doi|

  .. details:: APA
     :open:

     **Musselman, E. D.**, **Cariello, J. E.**, Grill, W. M., & Pelot, N. A. (2023). wmglab-duke/ascent: ASCENT v1.2.2 (v1.2.2) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.8298703

  .. details:: MLA

      Musselman, Eric D., et al. Wmglab-Duke/Ascent: ASCENT v1.2.2. v1.2.2, Zenodo, 2023, doi:10.5281/ZENODO.8298703.


  .. details:: BibTeX

    .. code-block:: BibTeX

        @misc{https://doi.org/10.5281/zenodo.8298703,
          doi = {10.5281/ZENODO.8298703},
          url = {https://zenodo.org/record/8298703},
          author = {Musselman,  Eric D and Cariello,  Jake E and Grill,  Warren M and Pelot,  Nicole A},
          title = {wmglab-duke/ascent: ASCENT v1.2.2},
          publisher = {Zenodo},
          year = {2023},
          copyright = {MIT License}
        }

**ASCENT** is an open source platform for simulating peripheral nerve stimulation. To download the software, visit the `ASCENT GitHub repository <https://github.com/wmglab-duke/ascent>`_.

..  youtube:: rG-KU7wWcXY

.. toctree::
   :maxdepth: 2
   :caption: Basic ASCENT Usage

   Getting_Started
   Running_ASCENT/index
   JSON/index
   Publishing_with_ASCENT/index

.. toctree::
   :maxdepth: 2
   :caption: Advanced ASCENT Usage

   MockSample
   Primitives_and_Cuffs/index
   Modeling_Neural_Recording
   Convergence_Example
   Troubleshooting-Guide

.. toctree::
   :maxdepth: 2
   :caption: ASCENT Hierarchies

   Code_Hierarchy/index
   Data_Hierarchy

.. toctree::
   :maxdepth: 2
   :caption: Reference

   Publications_Using_ASCENT
   references
   Validation
   Changelog

.. toctree::
   :maxdepth: 2
   :caption: External Links

   ASCENT Publication <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009285>
   ASCENT on GitHub <https://github.com/wmglab-duke/ascent>
   The Grill Lab <https://grill-lab.pratt.duke.edu>
   NIH SPARC <https://commonfund.nih.gov/sparc>
