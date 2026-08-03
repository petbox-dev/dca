from importlib.metadata import version as _get_version

__version__ = _get_version('petbox-dca')

from .base import (get_time, get_time_monthly_vol,
                   DeclineCurve, PrimaryPhase,
                   AssociatedPhase, BothAssociatedPhase,
                   SecondaryPhase, WaterPhase,
                   DAYS_PER_MONTH, DAYS_PER_YEAR)
from .primary import (NullPrimaryPhase, MultisegmentHyperbolic, MH, THM,
                      HyperbolicSegment, GeneralizedHyperbolic, IncliningHyperbolic,
                      PLE, SE, Duong)
from .associated import (NullAssociatedPhase, MultisegmentPLYield,
                         PLYield, PLYieldSegment, GeneralizedPLYield)
from .bourdet import bourdet
