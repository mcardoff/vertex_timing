# Vertex Timing Analysis

HGTD vertex timing reconstruction and clustering analysis for ATLAS.

See `CLAUDE.md` for full project documentation.

## Disabled Event Selection Cuts

The following cuts are commented out in `src/event_processing.h` and can be re-enabled for specific studies:

- **VBF signal region** (`branch->passVBFSignalRegion()`): Applies VBF H→invisible selection cuts. Disabled to run over the full inclusive sample; re-enable for VBF-specific efficiency studies.
- **Forward HS track requirement** (`branch->pass_forward_hs_tracks(nForwardTrack_HS)`): Requires a minimum number of hard-scatter tracks in the forward region. Useful for studies of track multiplicity dependence; disabled by default to retain low-track events in the denominator.
