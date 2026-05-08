"""Shared matplotlib helpers for the quality-check plots.

Lazy-imports matplotlib / seaborn so unit tests on the rest of the package
don't pay the import cost.
"""
from __future__ import annotations

import math
import os
from typing import List, Sequence


def _epoch_to_datetime_strings(epochs: Sequence[str]) -> List[str]:
    """Convert epoch labels ``"yy mm dd hh mi ss.s"`` to a string accepted
    by :func:`dateutil.parser.parse`."""
    return [
        epoch.split()[2] + " " + ":".join(epoch.split()[3:6]) for epoch in epochs
    ]


def scatter_per_satellite(
    filename: str,
    title: str,
    prns: Sequence[str],
    epochs: Sequence[str],
    dataset,
    band_index: int,
    n_legend_columns: int,
    limit: float = 100.0,
    y_label: str | None = None,
) -> None:
    """Generic scatter plot of a per-satellite, per-epoch quantity.

    ``dataset`` has shape ``(n_epochs, n_sats, n_bands)``. Points whose
    absolute value exceeds ``limit`` (or whose value is exactly ``0``) are
    skipped — RINEX files frequently contain corrupted samples that swamp
    the y-axis if rendered.
    """
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    import seaborn as sns
    from dateutil import parser

    times = [parser.parse(s) for s in _epoch_to_datetime_strings(epochs)]

    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_title(title)

    palette = sns.color_palette("viridis", len(prns))
    ax.set_prop_cycle("color", palette)

    for sat_idx, prn in enumerate(prns):
        xs, ys = [], []
        for epoch_idx, t in enumerate(times):
            value = dataset[epoch_idx][sat_idx][band_index]
            if value != 0 and abs(value) < limit:
                xs.append(t)
                ys.append(value)
        ax.scatter(xs, ys, s=1, label=prn)

    plt.legend(
        bbox_to_anchor=(1.05, 0),
        loc=3,
        borderaxespad=0,
        ncol=n_legend_columns,
    )
    if y_label:
        ax.set_ylabel(y_label)
    plt.savefig(filename, dpi=360, bbox_inches="tight")
    plt.close(fig)


def signal_plot(opath: str, epochs: Sequence[str], satellites: Sequence[str]) -> None:
    """Per-PRN signal-presence raster (Figure 12 in the manual)."""
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from dateutil import parser

    times = [parser.parse(s) for s in _epoch_to_datetime_strings(epochs)]

    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %H:%M"))
    ax.set_title(os.path.basename(opath)[:-4] + " Signal Plot")
    ax.scatter(times, satellites, s=2, marker=".", alpha=0.75, c="g")
    ax.grid(linestyle="--", linewidth=0.5)
    ax.set_ylabel("SATELLITE NO")
    fig.set_size_inches(15, 7)
    ax.tick_params(axis="y", labelsize=5)
    plt.savefig(opath[:-4] + "SignalPlot.jpg", dpi=360)
    plt.close(fig)


def sky_plot(
    opath: str,
    azi_ele_records: Sequence[dict],
    prns: Sequence[str],
) -> None:
    """Polar plot of azimuth vs elevation for each GPS PRN."""
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(subplot_kw={"polar": True})
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location("N")
    ax.set_rlim(90, 0)
    ax.set_title(os.path.basename(opath)[:-4] + "SkyPlot")

    gps_prns = [p for p in prns if p.startswith("G")]
    if not gps_prns:
        plt.close(fig)
        return

    palette = plt.cm.viridis(np.linspace(0, 1, len(gps_prns)))
    ax.set_prop_cycle("color", palette)

    series: dict[str, list[tuple[float, float]]] = {p: [] for p in gps_prns}
    for record in azi_ele_records:
        for prn in gps_prns:
            entry = record.get(prn)
            if entry is not None:
                series[prn].append(entry)

    for prn in gps_prns:
        if not series[prn]:
            continue
        rho = [math.degrees(point[1]) for point in series[prn]]
        theta = [point[0] for point in series[prn]]
        ax.scatter(theta, rho, s=10, label=prn)

    plt.legend(bbox_to_anchor=(0, 1.05), borderaxespad=3)
    plt.savefig(opath[:-4] + "skyplot.png", dpi=360, bbox_inches="tight")
    plt.close(fig)


def cycle_slip_plot(
    opath: str,
    epochs: Sequence[str],
    cycleslips,
    prns: Sequence[str],
    n_legend_columns: int,
) -> None:
    """Plot the carrier-phase second-difference cycle-slip indicator."""
    name = os.path.basename(opath)
    scatter_per_satellite(
        filename=opath[:-4] + "CycleSlipCarrier.png",
        title=name[:-4] + " Cycle Slip of Carrier-Phase observations plot",
        prns=prns,
        epochs=epochs,
        dataset=cycleslips,
        band_index=3,
        n_legend_columns=n_legend_columns,
        y_label="cycle slip (m)",
    )
