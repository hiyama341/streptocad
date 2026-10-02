History
-------

1.0.1
~~~~~
Documentation only; no change to the library.

* Corrected the local setup instructions in ``README.md``, which predated the move to
  uv. They told readers to install ``requirements.txt`` into a conda environment, when
  that file is generated from ``uv.lock`` for the Docker image, and both the local and
  Docker sections ended with ``python3 application.py``, which does not exist at the
  repository root. Replaced with ``uv sync --group app`` and ``uv run --group app
  python web_app/application.py``; conda and pip remain as a documented fallback.
* Replaced the library usage example with one that finds and filters Cas9 guide RNAs,
  including the returned columns and real output.
* Released so that the project page on PyPI, which renders the README baked into the
  uploaded artifact, no longer shows the superseded instructions.

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
  one or two bases.** This corrects real errors and is not a change of preference.
  There are two independent causes, described below.

Primers are now designed with the melting-temperature model they are reported with
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
* Reaction conditions are passed to ``pydna.design.primer_design`` as keyword
  arguments. pydna accepts ``**kwargs`` but never forwards them to the melting
  temperature function, so the selected polymerase and primer concentration were
  silently discarded and design fell back to pydna's generic Taq-buffer defaults,
  while the reported melting temperatures came from NEB. Design and reporting were
  therefore using two different thermodynamic models.
* ``streptocad.cloning.pcr_simulation.make_amplicons`` additionally passed the
  argument as ``tm_function`` rather than ``tm_func``, so it was absorbed by
  ``**kwargs`` and the NEB model was never used there at all.
* All design call sites now bind the conditions into the Tm function through
  :func:`streptocad.primers.tm.neb_tm_function`, and pass ``estimate_function`` so a
  fast local estimate seeds the search. Results are identical to querying NEB at
  every step while making far fewer API calls, and responses are cached.
* Effect on the workflow-1 regulator set, scoring each designed primer with the NEB
  model it is reported with: mean deviation from the requested melting temperature
  falls from 1.69 C to 1.00 C, primers within 1 C of target rise from 14/36 to 26/36,
  and the most common value moves from 58 C to the requested 60 C. The previous
  design sat systematically below target; that bias is gone.
* The Gibson repair-template primers were affected more severely. There the discarded
  conditions left design running against NEB's ``q5-0`` default, which reads several
  degrees higher than Phusion GC, so three of the four primers never grew past the
  13 nt ``min_primer_length`` floor. For the SCO5892 deletion they came out near 50 C
  against a requested 60 C; they now read 60-61 C. Anyone who ordered Gibson repair
  primers from an earlier version may want to re-check them, since a primer designed
  for a 60 C anneal that actually melts at 50 C can underperform at the bench.
* ``streptocad.primers.tm.neb_tm_function`` raises a clear error when the NEB API
  returns no result, instead of yielding ``None`` and failing later as a ``TypeError``
  deep inside pydna's search. Note the API rejects primers shorter than 8 nt and any
  product code outside those in ``streptocad.utils.polymerase_dict``.

Corrected nearest-neighbour thermodynamic parameters
""""""""""""""""""""""""""""""""""""""""""""""""""""
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
