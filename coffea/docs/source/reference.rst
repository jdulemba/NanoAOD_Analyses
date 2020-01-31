API Reference Guide
*******************
Coffea: a column object framework for effective analysis.

When executing

    >>> import coffea

a subset of the full coffea package is imported into the python environment.
Some packages must be imported explicitly, so as to avoid importing unnecessary
and/or heavy dependencies.  Below lists the packages available in the ``coffea`` namespace.

.. autosummary::
    :toctree: modules
    :template: automodapi_templ.rst

    coffea.analysis_objects
    coffea.arrays
    coffea.nanoaod
    coffea.nanoaod.methods
    coffea.lookup_tools
    coffea.jetmet_tools
    coffea.lumi_tools
    coffea.hist
    coffea.util
    coffea.processor
