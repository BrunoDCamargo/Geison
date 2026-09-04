"""Primer-design API backed by the shared region-selection module."""

from __future__ import annotations

from . import legacy as _legacy
from qpcr_pipeline.region_selection import (
    CandidateRegion,
    RegionSelectionError,
    candidate_region_from_window,
    is_target_eligible,
    overlap_fraction,
    select_conservation_candidate_regions,
)


for _name in dir(_legacy):
    if not _name.startswith("_") and _name != "CandidateRegion":
        globals()[_name] = getattr(_legacy, _name)

PrimerDesignError = _legacy.PrimerDesignError


def _select_candidate_regions(conservation, config):
    try:
        return select_conservation_candidate_regions(conservation, config)
    except RegionSelectionError as error:
        raise PrimerDesignError(str(error)) from error


# The legacy implementation resolves these globals at call time. Point it at
# the shared primitives so public design behavior stays unchanged while the
# geometry has one source of truth.
_legacy.CandidateRegion = CandidateRegion
_legacy._select_candidate_regions = _select_candidate_regions

__all__ = [
    name
    for name in globals()
    if not name.startswith("_")
]
