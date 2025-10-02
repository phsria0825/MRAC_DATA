"""
ALL-μ KDE-HDR (50/80/95%) + Utilization(λ) Scatter
- Inputs:  LQR_TV_SWA_Tire.csv, MRAC_TV_SWA_Tire.csv
- Outputs: Individual figures (PNG + PDF + SVG) for:
    LQR_FRONT, LQR_REAR, MRAC_FRONT, MRAC_REAR
  And a 1x4 panel (PNG + PDF + SVG) composed as a single vector figure.

Notes
- Normalized coordinates: (Fx/(μFz), Fy/(μFz))
- KDE on standardized space (StandardScaler), Scott bandwidth (n^(-1/6))
- HDR thresholds: 50%, 80%, 95%
- Scatter color = λ = sqrt(Fx^2 + Fy^2) / (μ|Fz|)
"""

from __future__ import annotations
import os
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn.preprocessing import StandardScaler

# =========================
# Config (adjust as needed)
# =========================
LQR_CSV  = "LQR_TV_SWA_Tire.csv"
MRAC_CSV = "MRAC_TV_SWA_Tire.csv"
OUT_DIR  = "fig_out"                     # output folder

LIM: float = 2.0                         # axis limit (±LIM)
GRID_SIZE: int = 200                     # KDE grid resolution
HDR_LEVELS = (0.5, 0.8, 0.95)            # HDR mass levels
SCATTER_MAX_N: int = 3000                # scatter subsample size
SCATTER_SIZE: float = 8.0                # scatter dot size
SCATTER_ALPHA: float = 0.35              # scatter alpha
PNG_DPI: int = 180                       # for PNG only

TITLE_FONT = 14
LABEL_FONT = 12
TICK_FONT  = 10
GRID_ALPHA = 0.6

# ===================================
# Utilities: I/O, math, KDE-HDR, etc.
# ===================================

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def save_all_formats(fig: plt.Figure, base_path_wo_ext: str) -> None:
    """Save figure as PNG, PDF, and SVG (PDF/SVG are true vector)."""
    fig.savefig(f"{base_path_wo_ext}.png", dpi=PNG_DPI, bbox_inches="tight")
    fig.savefig(f"{base_path_wo_ext}.pdf", bbox_inches="tight")
    fig.savefig(f"{base_path_wo_ext}.svg", bbox_inches="tight")

def load_xy_lambda(csv_path: str, wheels: Tuple[str, str]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load and return:
      X: (N,2) of [Fx/(μFz), Fy/(μFz)] concatenated for given wheels
      lam: (N,) utilization lambda for each point
    wheels example: ("FL","FR") or ("RL","RR")
    """
    df = pd.read_csv(csv_path)
    eps = 1e-9
    xs, ys, lams = [], [], []
    for w in wheels:
        Fx = df[f"{w}_x"].astype(float).values
        Fy = df[f"{w}_y"].astype(float).values
        Fz = np.abs(df[f"{w}_z"].astype(float).values)
        mu = df[f"Road_{w}_mu"].astype(float).values

        denom = np.maximum(mu * Fz, eps)
        xs.append(Fx / denom)
        ys.append(Fy / denom)
        lams.append(np.sqrt(Fx*Fx + Fy*Fy) / denom)

    X = np.column_stack([np.concatenate(xs), np.concatenate(ys)])
    lam = np.concatenate(lams)

    good = np.isfinite(X).all(axis=1) & np.isfinite(lam)
    return X[good], lam[good]

@dataclass
class KDEGrid:
    XX: np.ndarray
    YY: np.ndarray
    Z: np.ndarray
    thr: Dict[float, float]  # q -> density threshold
    dx: float
    dy: float

def kde_hdr_grid(
    X: np.ndarray,
    lim: float = LIM,
    gridsize: int = GRID_SIZE,
    levels = HDR_LEVELS
) -> KDEGrid:
    """
    Fit KDE on standardized space, evaluate density on uniform grid (original axes),
    normalize density to integrate to 1, and compute HDR thresholds.
    """
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)

    # Scott bandwidth (d=2 => n^(-1/6))
    n = Xs.shape[0]
    bw = n ** (-1.0 / 6.0)
    kde = KernelDensity(kernel="gaussian", bandwidth=bw).fit(Xs)

    xs = np.linspace(-lim, lim, gridsize)
    ys = np.linspace(-lim, lim, gridsize)
    XX, YY = np.meshgrid(xs, ys)
    grid_xy  = np.column_stack([XX.ravel(), YY.ravel()])
    grid_xys = scaler.transform(grid_xy)

    # raw KDE and normalization to probability density
    p = np.exp(kde.score_samples(grid_xys)).reshape(gridsize, gridsize)
    dx = (2 * lim) / (gridsize - 1)
    dy = (2 * lim) / (gridsize - 1)
    Z  = p / (p.sum() * dx * dy)

    # HDR thresholds (mass above threshold equals q)
    flat  = Z.ravel()
    order = np.argsort(flat)[::-1]
    csum  = np.cumsum(flat[order] * dx * dy)

    thr: Dict[float, float] = {}
    for q in levels:
        idx = np.searchsorted(csum, q)
        thr[q] = flat[order[idx]] if idx < len(order) else flat[order[-1]]

    return KDEGrid(XX=XX, YY=YY, Z=Z, thr=thr, dx=dx, dy=dy)

def sorted_levels(thr: Dict[float, float], levels = HDR_LEVELS) -> np.ndarray:
    """Return strictly increasing contour levels from HDR thresholds."""
    vals = sorted({thr[q] for q in levels})
    if len(vals) == 1:  # pathological guard
        vals = [vals[0]*0.999, vals[0]*1.001]
    return np.asarray(vals, dtype=float)

def subsample(X: np.ndarray, y: np.ndarray, k: int = SCATTER_MAX_N, seed: int = 42) -> Tuple[np.ndarray, np.ndarray]:
    """Uniform random subsample for scatter readability."""
    n = len(X)
    if n <= k: 
        return X, y
    rng = np.random.default_rng(seed)
    idx = rng.choice(n, size=k, replace=False)
    return X[idx], y[idx]

# ==========================
# Drawing (matplotlib, 100% vector)
# ==========================

def style_axes(ax: plt.Axes, title: str, lim: float = LIM) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.grid(True, linestyle=":", alpha=GRID_ALPHA)
    ax.set_xlabel("Fx / (μ·Fz)", fontsize=LABEL_FONT)
    ax.set_ylabel("Fy / (μ·Fz)", fontsize=LABEL_FONT)
    ax.set_title(title, fontsize=TITLE_FONT)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontsize(TICK_FONT)
    # unit circle
    t = np.linspace(0, 2*np.pi, 360)
    ax.plot(np.cos(t), np.sin(t), linewidth=1.3)

def draw_hdr_scatter(ax: plt.Axes, X: np.ndarray, lam: np.ndarray) -> None:
    """Draw HDR (50/80/95) with light fill + contour, and utilization scatter overlay."""
    grid = kde_hdr_grid(X)
    levels = sorted_levels(grid.thr)
    fill_levels = list(levels) + [grid.Z.max() * 1.001]

    # light fill and contour lines (default colormap)
    ax.contourf(grid.XX, grid.YY, grid.Z, levels=fill_levels, alpha=0.28)
    cs = ax.contour(grid.XX, grid.YY, grid.Z, levels=levels, linewidths=1.8)
    # label contours with HDR percentages
    fmt = {grid.thr[q]: f"{int(q*100)}%" for q in HDR_LEVELS if grid.thr[q] in cs.levels}
    try:
        ax.clabel(cs, inline=1, fmt=fmt, fontsize=9)
    except Exception:
        pass

    # utilization scatter (slightly darker as requested)
    Xs, ls = subsample(X, lam, k=SCATTER_MAX_N)
    ax.scatter(Xs[:, 0], Xs[:, 1], s=SCATTER_SIZE, c=ls, alpha=SCATTER_ALPHA, zorder=3)

def render_single(X: np.ndarray, lam: np.ndarray, title: str, out_base: str) -> None:
    """Render one figure and save as PNG/PDF/SVG."""
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    draw_hdr_scatter(ax, X, lam)
    style_axes(ax, title)
    fig.tight_layout()
    save_all_formats(fig, out_base)
    plt.close(fig)

def render_panel_horizontal(data_map: Dict[str, Tuple[np.ndarray, np.ndarray]], out_base: str) -> None:
    """
    Render a 1x4 panel (true vector): keys order defines position.
    data_map: Ordered dict-like {label: (X, lam), ...} expected length 4.
    """
    labels = list(data_map.keys())
    fig, axs = plt.subplots(nrows=1, ncols=len(labels), figsize=(6.2*len(labels), 6.0), squeeze=True)
    if len(labels) == 1:
        axs = [axs]  # safety for single axis
    for ax, lab in zip(axs, labels):
        X, lam = data_map[lab]
        draw_hdr_scatter(ax, X, lam)
        # Put concise titles (panel already labeled in caption)
        style_axes(ax, lab)
    fig.tight_layout()
    save_all_formats(fig, out_base)
    plt.close(fig)

# ===========
# Main script
# ===========
def main() -> None:
    ensure_dir(OUT_DIR)

    # Load data
    X_lqr_front, lam_lqr_front = load_xy_lambda(LQR_CSV,  ("FL", "FR"))
    X_lqr_rear,  lam_lqr_rear  = load_xy_lambda(LQR_CSV,  ("RL", "RR"))
    X_mrac_front,lam_mrac_front= load_xy_lambda(MRAC_CSV, ("FL", "FR"))
    X_mrac_rear, lam_mrac_rear = load_xy_lambda(MRAC_CSV, ("RL", "RR"))

    # Individual figures (PNG + PDF + SVG)
    render_single(X_lqr_front, lam_lqr_front, "LQR_TV — FRONT (HDR + utilization scatter)",
                  os.path.join(OUT_DIR, "LQR_ALL_front_cmap_scatter"))
    render_single(X_lqr_rear,  lam_lqr_rear,  "LQR_TV — REAR (HDR + utilization scatter)",
                  os.path.join(OUT_DIR, "LQR_ALL_rear_cmap_scatter"))
    render_single(X_mrac_front,lam_mrac_front,"MRAC_TV — FRONT (HDR + utilization scatter)",
                  os.path.join(OUT_DIR, "MRAC_ALL_front_cmap_scatter"))
    render_single(X_mrac_rear, lam_mrac_rear, "MRAC_TV — REAR (HDR + utilization scatter)",
                  os.path.join(OUT_DIR, "MRAC_ALL_rear_cmap_scatter"))

    # 1x4 vector panel (LQR Front | LQR Rear | MRAC Front | MRAC Rear)
    data_map = {
        "LQR — Front":  (X_lqr_front,  lam_lqr_front),
        "LQR — Rear":   (X_lqr_rear,   lam_lqr_rear),
        "MRAC — Front": (X_mrac_front, lam_mrac_front),
        "MRAC — Rear":  (X_mrac_rear,  lam_mrac_rear),
    }
    render_panel_horizontal(data_map, os.path.join(OUT_DIR, "ALLmu_panel_1x4_scatter"))

if __name__ == "__main__":
    main()
