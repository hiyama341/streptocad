from Bio.Seq import Seq
from pydna.dseqrecord import Dseqrecord


def extract_feature_sequence(feature, sequence):
    """Extract a feature from pydna/Biopython sequence objects across Biopython versions."""
    return feature.extract(Seq(str(sequence)))


def patch_dseqrecord_linear_argument():
    """Allow legacy Dseqrecord(..., linear=True) calls with newer pydna."""
    if getattr(Dseqrecord, "_streptocad_linear_patch", False):
        return

    original_init = Dseqrecord.__init__

    def patched_init(self, record, *args, linear=None, circular=None, **kwargs):
        if linear is not None and circular is None:
            circular = not linear
        original_init(self, record, *args, circular=circular, **kwargs)

    Dseqrecord.__init__ = patched_init
    Dseqrecord._streptocad_linear_patch = True
