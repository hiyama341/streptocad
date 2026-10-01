"""Tests for the melting-temperature helper used by primer design.

The point of ``neb_tm_function`` is that reaction conditions survive the trip into
pydna's ``primer_design``, which accepts ``**kwargs`` but never forwards them to the
Tm function. These tests pin that behaviour down so the conditions cannot silently
start being ignored again.
"""

import pytest

from streptocad.primers.tm import (
    clear_tm_cache,
    neb_tm_function,
    tm_cache_info,
    tm_default,
)

PRIMER = "GACGATTCGGCCCGTGC"


def test_returns_a_single_argument_callable():
    """pydna calls tm_func with one argument, so the conditions must be bound."""
    tm = neb_tm_function(conc=0.4, prodcode="phusion-1")
    assert callable(tm)
    # One positional argument only; anything else means the conditions leaked back
    # into the call signature and would be dropped by pydna.
    import inspect

    assert len(inspect.signature(tm).parameters) == 1


def test_distinct_conditions_give_distinct_functions():
    a = neb_tm_function(conc=0.4, prodcode="phusion-1")
    b = neb_tm_function(conc=0.4, prodcode="q5-0")
    assert a is not b
    assert a.__name__ != b.__name__


def test_tm_default_is_importable_for_use_as_estimate_function():
    """primer_design takes this as estimate_function to seed the search."""
    assert callable(tm_default)
    assert tm_default(PRIMER) > 0


@pytest.mark.integration
def test_conditions_actually_reach_the_api():
    """Different polymerases must give different melting temperatures.

    If the product code were being discarded, these would be equal. That is exactly
    the failure this module exists to prevent.
    """
    phusion = neb_tm_function(conc=0.4, prodcode="phusion-1")(PRIMER)
    q5 = neb_tm_function(conc=0.4, prodcode="q5-0")(PRIMER)
    assert phusion != q5, (
        "Phusion and Q5 returned the same Tm, which means the product code is not "
        "reaching the NEB API."
    )


@pytest.mark.integration
def test_repeated_lookups_are_cached():
    clear_tm_cache()
    tm = neb_tm_function(conc=0.4, prodcode="phusion-1")
    first = tm(PRIMER)
    before = tm_cache_info()
    second = tm(PRIMER)
    after = tm_cache_info()
    assert first == second
    assert after.hits == before.hits + 1
    assert after.misses == before.misses


@pytest.mark.integration
def test_unusable_request_raises_instead_of_returning_none():
    """A bad product code must fail loudly.

    primer_tm_neb returns None when NEB rejects a request; left alone that surfaces
    far away as a TypeError inside pydna when the result is compared with a number.
    """
    tm = neb_tm_function(conc=0.4, prodcode="phusion")  # not a valid code
    with pytest.raises(RuntimeError, match="NEB melting temperature API"):
        tm(PRIMER)
