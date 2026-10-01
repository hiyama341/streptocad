"""Melting-temperature functions used for primer design.

pydna's ``primer_design`` takes the Tm function as ``tm_func`` and calls it with a
single argument. It accepts ``**kwargs`` but never forwards them to ``tm_func``, so
reaction conditions passed that way are silently discarded. Binding them here with a
closure is what makes the chosen polymerase and primer concentration actually reach
the calculation.

Use :func:`neb_tm_function` for both design and reporting so the two cannot drift
apart: designing against one model while reporting another produces primers whose
reported Tm does not match the requested ``target_tm``.
"""

from functools import lru_cache
from typing import Callable

from pydna.tm import tm_default
from teemi.build.PCR import primer_tm_neb

__all__ = ["neb_tm_function", "tm_default", "clear_tm_cache", "tm_cache_info"]


@lru_cache(maxsize=8192)
def _cached_neb_tm(primer: str, conc: float, prodcode: str) -> float:
    """One NEB API call per distinct (primer, conc, prodcode).

    Primer design walks outwards one base at a time and re-evaluates overlapping
    candidates, so caching cuts the number of requests substantially.

    Raises
    ------
    RuntimeError
        If NEB returns no melting temperature. ``primer_tm_neb`` returns ``None``
        in that case, which would otherwise surface far away as a ``TypeError``
        inside pydna's search loop when the result is compared with a number.
    """
    tm = primer_tm_neb(primer, conc=conc, prodcode=prodcode)
    if tm is None:
        raise RuntimeError(
            f"The NEB melting temperature API returned no result for {primer!r} "
            f"(length {len(primer)}, conc={conc}, prodcode={prodcode!r}). "
            "Check that the product code is one of the values in "
            "streptocad.utils.polymerase_dict (for example 'phusion-1', not "
            "'phusion'), that the primer is at least 8 nt, and that the API is "
            "reachable. See https://tmapi.neb.com/docs/productcodes"
        )
    return tm


def neb_tm_function(conc: float = 0.4, prodcode: str = "q5-0") -> Callable[[str], float]:
    """Return a single-argument Tm function bound to these reaction conditions.

    Parameters
    ----------
    conc : float
        Primer concentration in micromolar.
    prodcode : str
        NEB product code for the polymerase and buffer, for example ``"phusion-1"``.
        Codes are listed at https://tmapi.neb.com/docs/productcodes

    Returns
    -------
    callable
        A function of one argument suitable for pydna's ``tm_func``.

    Examples
    --------
    >>> from pydna.design import primer_design                    # doctest: +SKIP
    >>> tm = neb_tm_function(conc=0.4, prodcode="phusion-1")      # doctest: +SKIP
    >>> amplicon = primer_design(                                 # doctest: +SKIP
    ...     template, target_tm=65, limit=10,
    ...     tm_func=tm, estimate_function=tm_default,
    ... )
    """

    def tm(primer) -> float:
        return _cached_neb_tm(str(primer), conc, prodcode)

    tm.__name__ = f"neb_tm({prodcode}, {conc}uM)"
    return tm


def clear_tm_cache() -> None:
    """Drop cached NEB melting temperatures."""
    _cached_neb_tm.cache_clear()


def tm_cache_info():
    """Return cache statistics, useful for checking how many API calls were saved."""
    return _cached_neb_tm.cache_info()
