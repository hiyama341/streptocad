History
-------

1.0.0
~~~~~
* First release of the ``streptocad`` library on PyPI (``pip install streptocad``).
* Packaged the ``streptocad`` toolbox with a hatchling build backend; the Dash web
  app, notebooks and test fixtures stay in the repository and are not distributed.
* Declared ``nbformat`` and ``nbconvert`` as runtime dependencies, since
  ``streptocad.utils`` imports them at module level.
* Relaxed dependency pins to compatible ranges for library use. The deployed web
  app keeps its exact pins in ``requirements.txt``.
