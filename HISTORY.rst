History
-------

1.0.0
~~~~~

Packaging
^^^^^^^^^
* First release of the ``streptocad`` library on PyPI (``pip install streptocad``).
* Packaged the ``streptocad`` toolbox with a hatchling build backend; the Dash web
  app, notebooks and test fixtures stay in the repository and are not distributed.
* Declared ``nbformat`` and ``nbconvert`` as runtime dependencies, since
  ``streptocad.utils`` imports them at module level.
* Relaxed dependency pins to compatible ranges for library use, so security
  updates can be picked up with ``uv lock --upgrade``. ``uv.lock`` still records
  exact versions, and ``requirements.txt`` is generated from it for the Docker image.

Corrected primer melting temperatures (changes designed primers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* **Designed primers may differ from those produced by earlier versions, usually by
  one base.** This corrects a real error and is not a change of preference.
* Biopython 1.82 corrected two entropy values in the SantaLucia & Hicks (2004)
  nearest-neighbour table (``DNA_NN4``): ``TA/AT`` from -20.4 to -21.3 and ``GG/CC``
  from -19.0 to -19.9. The corrected values reproduce the published dG37 for those
  pairs to within 0.02 kcal/mol; the earlier ones are off by about 0.28.
* StreptoCAD previously required ``biopython==1.80`` (via ``teemi`` 0.3.4), so it
  used the erroneous values. They overestimate primer melting temperatures by
  roughly 2.3 C, which made ``pydna``'s ``primer_design`` stop one base short of the
  requested ``target_tm``. Primers designed before this release therefore tend to be
  one base too short and anneal about 2 C below the temperature that was asked for.
* Raised the floors to ``biopython>=1.82``, ``teemi>=1.0.5`` and ``pydna>=5.5``.
* Updated the primer fixtures in ``tests/`` accordingly. Each expected primer was
  checked to be the candidate whose Tm lies closest to the requested ``target_tm``
  among its neighbouring lengths, rather than simply re-recording the new output.
* Users reproducing designs from before this release, including those in the
  StreptoCAD paper, should expect slightly different primer sequences. The earlier
  sequences remain valid primers; they were simply designed against a Tm that was
  overestimated.
