

:Python version: Python 3.6, 3.7.3; most modules are Python2.7 compatible.


User guide and reference
###########################

Installation
============

::

    git clone https://github.com/GillesFischerSorbonne/dubii_project/
    cd dubii_project/
    python setupy.py install

Usage
=====

From a python shell:

.. plot::
    :include-source:

    from telofinder import get_data
    from telofinder import analyze_telom_length  as alt
    filename = get_data("test.telo.blocks.fasta")

    df = alt.run_telofinder(filename)
    df.loc["test.telo.blocks", "utg84"][["entropy", "pattern"]].loc[0:4000].plot()



Reference
=========

.. automodule:: telofinder.analyze_telom_length
    :members:
    :undoc-members:





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

