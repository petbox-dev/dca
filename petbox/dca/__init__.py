from importlib.metadata import version as _get_version

__version__ = _get_version("petbox-dca")

from .associated import (
    GeneralizedPLYield,
    MultisegmentPLYield,
    NullAssociatedPhase,
    PLYield,
    PLYieldSegment,
)
from .base import (
    DAYS_PER_MONTH,
    DAYS_PER_YEAR,
    AssociatedPhase,
    BothAssociatedPhase,
    DeclineCurve,
    PrimaryPhase,
    SecondaryPhase,
    WaterPhase,
    get_time,
    get_time_monthly_vol,
)
from .bourdet import bourdet
from .primary import (
    MH,
    PLE,
    SE,
    THM,
    Duong,
    GeneralizedHyperbolic,
    Hyperbolic,
    HyperbolicSegment,
    IncliningHyperbolic,
    MultisegmentHyperbolic,
    NullPrimaryPhase,
)

# The public API, stated explicitly. This package ships `py.typed`, so a consumer's own
# `mypy --strict` run type-checks against it -- and strict enables `--no-implicit-reexport`,
# under which a name merely imported into this module is NOT re-exported. Without `__all__`,
# every `dca.MH(...)` in downstream code raised
# `Module "petbox.dca" does not explicitly export attribute "MH"`.
#
# Anything listed here is API. `ParamDesc` is deliberately absent: it is documented as
# `petbox.dca.base.ParamDesc` and reached through the submodule.
# Grouped by role rather than sorted (RUF022): the groups are what make a 28-name list
# readable, and a sort would interleave models, interfaces and constants.
__all__ = [  # noqa: RUF022
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
