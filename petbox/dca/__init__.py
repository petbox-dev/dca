from importlib.metadata import version as _get_version

__version__ = _get_version("petbox-dca")

from .base import (
    get_time,
    get_time_monthly_vol,
    DeclineCurve,
    PrimaryPhase,
    AssociatedPhase,
    BothAssociatedPhase,
    SecondaryPhase,
    WaterPhase,
    DAYS_PER_MONTH,
    DAYS_PER_YEAR,
)
from .primary import (
    NullPrimaryPhase,
    MultisegmentHyperbolic,
    Hyperbolic,
    MH,
    THM,
    HyperbolicSegment,
    GeneralizedHyperbolic,
    IncliningHyperbolic,
    PLE,
    SE,
    Duong,
)
from .associated import (
    NullAssociatedPhase,
    MultisegmentPLYield,
    PLYield,
    PLYieldSegment,
    GeneralizedPLYield,
)
from .bourdet import bourdet

# The public API, stated explicitly. This package ships `py.typed`, so a consumer's own
# `mypy --strict` run type-checks against it -- and strict enables `--no-implicit-reexport`,
# under which a name merely imported into this module is NOT re-exported. Without `__all__`,
# every `dca.MH(...)` in downstream code raised
# `Module "petbox.dca" does not explicitly export attribute "MH"`.
#
# Anything listed here is API. `ParamDesc` is deliberately absent: it is documented as
# `petbox.dca.base.ParamDesc` and reached through the submodule.
__all__ = [
    "__version__",
    # helpers
    "get_time",
    "get_time_monthly_vol",
    "bourdet",
    # constants
    "DAYS_PER_MONTH",
    "DAYS_PER_YEAR",
    # abstract interfaces
    "DeclineCurve",
    "PrimaryPhase",
    "AssociatedPhase",
    "BothAssociatedPhase",
    "SecondaryPhase",
    "WaterPhase",
    # primary phase models
    "NullPrimaryPhase",
    "MultisegmentHyperbolic",
    "Hyperbolic",
    "MH",
    "THM",
    "GeneralizedHyperbolic",
    "IncliningHyperbolic",
    "PLE",
    "SE",
    "Duong",
    # associated phase models
    "NullAssociatedPhase",
    "MultisegmentPLYield",
    "PLYield",
    "GeneralizedPLYield",
    # segment types for the generalized models
    "HyperbolicSegment",
    "PLYieldSegment",
]
