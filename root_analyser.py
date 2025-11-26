#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Geant4 ROOT zırhlama analiz GUI'si (PyQt5 + uproot + matplotlib)
Gelişmiş Sürüm (TRACKS + GEOMETRY + PDG + Çoklu Grafikler + Gelişmiş 3D Hit Haritası + Event Filtresi)
"""
import sys
import os
import re
from typing import List, Tuple, Optional, Dict

import numpy as np
import uproot

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QFileDialog,
    QListWidget,
    QListWidgetItem,
    QLineEdit,
    QLabel,
    QFormLayout,
    QSpinBox,
    QDoubleSpinBox,
    QTextEdit,
    QTabWidget,
    QMessageBox,
    QSizePolicy,
    QCheckBox,
    QTableWidget,
    QTableWidgetItem,
    QComboBox,
)
from PyQt5.QtCore import Qt


PDG_NAME_MAP = {
    11: "e⁻",
    -11: "e⁺",
    13: "μ⁻",
    -13: "μ⁺",
    22: "γ",
    2112: "n",
    2212: "p",
    0: "optik foton",
}


def pdg_name(pdg: int) -> str:
    return PDG_NAME_MAP.get(pdg, f"PDG {pdg}")


def pdg_color(pdg: int) -> str:
    if pdg == 11:
        return "tab:blue"
    if pdg == -11:
        return "tab:cyan"
    if pdg == 22:
        return "tab:orange"
    if pdg == 2212:
        return "tab:red"
    if pdg == 2112:
        return "tab:green"
    if pdg in (13, -13):
        return "tab:purple"
    if pdg == 0:
        return "tab:pink"
    return "tab:gray"


def infer_thickness_from_filename(basename: str) -> Optional[float]:
    m = re.search(r"(\d+(?:\.\d*)?)mm", basename, flags=re.IGNORECASE)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def read_tree(filename: str, tree_name: str):
    f = uproot.open(filename)
    keys = list(f.keys())
    clean_keys = [k.split(";")[0] for k in keys]

    ttree_names = []
    for k, obj in f.items():
        if isinstance(obj, uproot.behaviors.TTree.TTree):
            ttree_names.append(k.split(";")[0])

    if tree_name not in clean_keys:
        msg_lines = [
            f"'{filename}' içinde '{tree_name}' isimli TTree bulunamadı.",
            "Mevcut TTree'ler:",
        ]
        if ttree_names:
            for n in ttree_names:
                msg_lines.append(f"  - {n}")
        else:
            msg_lines.append("  (Bu dosyada TTree bulunamadı)")
        raise RuntimeError("\n".join(msg_lines))

    for k in keys:
        if k.startswith(tree_name + ";"):
            tree = f[k]
            if not isinstance(tree, uproot.behaviors.TTree.TTree):
                raise RuntimeError(
                    f"'{tree_name}' bulundu ama TTree değil (tip: {type(tree)})."
                )
            return tree

    raise RuntimeError(f"'{tree_name}' anahtarı çözümlenemedi (ROOT key).")


def list_ttrees_in_file(filename: str) -> List[str]:
    f = uproot.open(filename)
    names = []
    for k, obj in f.items():
        if isinstance(obj, uproot.behaviors.TTree.TTree):
            names.append(k.split(";")[0])
    return names


def list_branches_in_tree(filename: str, tree_name: str) -> List[str]:
    tree = read_tree(filename, tree_name)
    return list(tree.keys())


def get_branch(tree, branch_name: str) -> np.ndarray:
    try:
        arr = tree[branch_name].array(library="np")
    except Exception as e:
        raise RuntimeError(
            f"TTree içinde '{branch_name}' branch'i okunamadı: {e}\n"
            f"Mevcut branch'leri görmek için GUI'deki 'Branch'leri listele' butonunu kullan."
        )
    return arr


def make_energy_spectrum(edep: np.ndarray, nbins: int = 200,
                         emin: Optional[float] = None,
                         emax: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
    edep = np.asarray(edep, dtype=float)
    mask = np.isfinite(edep) & (edep >= 0.0)
    edep = edep[mask]

    if edep.size == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 0.0])

    if emin is None:
        emin = float(np.min(edep))
    if emax is None:
        emax = float(np.max(edep))

    if np.isclose(emin, emax):
        emax = emin + 1.0

    hist, edges = np.histogram(edep, bins=nbins, range=(emin, emax))
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, hist


def make_depth_dose(z: np.ndarray, edep: np.ndarray,
                    nbins: int = 100,
                    zmin: Optional[float] = None,
                    zmax: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
    z = np.asarray(z, dtype=float)
    edep = np.asarray(edep, dtype=float)

    mask = np.isfinite(z) & np.isfinite(edep)
    z = z[mask]
    edep = edep[mask]

    if z.size == 0:
        return np.array([]), np.array([])

    if zmin is None:
        zmin = float(np.min(z))
    if zmax is None:
        zmax = float(np.max(z))

    if np.isclose(zmin, zmax):
        zmax = zmin + 1.0

    bins = np.linspace(zmin, zmax, nbins + 1)
    dose, edges = np.histogram(z, bins=bins, weights=edep)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, dose


def make_2d_dose_map(x: np.ndarray, z: np.ndarray, edep: np.ndarray,
                     nx: int = 80, nz: int = 80) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    z = np.asarray(z, dtype=float)
    w = np.asarray(edep, dtype=float)

    mask = np.isfinite(x) & np.isfinite(z) & np.isfinite(w)
    x = x[mask]
    z = z[mask]
    w = w[mask]

    if x.size == 0 or z.size == 0:
        return np.zeros((2, 2)), np.array([0, 1]), np.array([0, 1])

    H, xedges, zedges = np.histogram2d(x, z, bins=[nx, nz], weights=w)
    return H.T, xedges, zedges


def compute_total_output_from_window(edep: np.ndarray,
                                     z: Optional[np.ndarray] = None,
                                     zmin: Optional[float] = None,
                                     zmax: Optional[float] = None) -> float:
    edep = np.asarray(edep, dtype=float)
    mask = np.isfinite(edep)

    if z is not None and zmin is not None and zmax is not None and zmax > zmin:
        z = np.asarray(z, dtype=float)
        mask = mask & np.isfinite(z) & (z >= zmin) & (z <= zmax)

    return float(np.sum(edep[mask]))


def fit_attenuation(thickness_mm: List[float],
                    outputs: List[float]) -> Tuple[float, float, float, np.ndarray, float, np.ndarray]:
    t_full = np.asarray(thickness_mm, dtype=float)
    I_full = np.asarray(outputs, dtype=float)

    mask = np.isfinite(t_full) & np.isfinite(I_full) & (I_full > 0.0)
    t = t_full[mask]
    I = I_full[mask]

    if t.size < 2:
        raise RuntimeError("Geçerli noktası (I>0) olan en az 2 dosya yok.")
    if np.unique(t).size < 2:
        raise RuntimeError("Tüm kalınlıklar aynı, HVL/TVL fit'i anlamsız.")

    idx = np.argsort(t)
    t = t[idx]
    I = I[idx]

    I0 = I[0]
    if I0 <= 0:
        positives = np.where(I > 0)[0]
        if positives.size == 0:
            raise RuntimeError("Pozitif I yok, HVL/TVL fit'i yapılamaz.")
        I0 = I[positives[0]]

    lnI_ratio = np.log(I / I0)
    coeffs = np.polyfit(t, lnI_ratio, 1)
    b = coeffs[0]
    a = coeffs[1]

    mu = -b
    HVL = np.log(2.0) / mu
    TVL = np.log(10.0) / mu

    t_fit = np.linspace(float(np.min(t)), float(np.max(t)), 200)
    lnI_fit = a + b * t_fit
    I_fit = I0 * np.exp(lnI_fit)

    fit_curve = np.vstack([t_fit, I_fit])
    return mu, HVL, TVL, fit_curve, I0, mask


def save_csv(filename: str, header: str, data: np.ndarray):
    np.savetxt(filename, data, delimiter=",", header=header, comments="")


def parse_int_list(text: str) -> List[int]:
    text = text.strip()
    if not text:
        return []
    parts = [p.strip() for p in text.replace(";", ",").split(",") if p.strip()]
    out = []
    for p in parts:
        try:
            out.append(int(p))
        except ValueError:
            continue
    return out


def autodetect_xyz_branches(branches: List[str]) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    lower_map = {b.lower(): b for b in branches}
    candidates_x = ["posx", "x", "prex", "postx"]
    candidates_y = ["posy", "y", "prey", "posty"]
    candidates_z = ["posz", "z", "prez", "postz"]

    def pick(cands):
        for c in cands:
            if c in lower_map:
                return lower_map[c]
        return None

    return pick(candidates_x), pick(candidates_y), pick(candidates_z)


class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure(figsize=(5, 4), dpi=120)
        super().__init__(fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.updateGeometry()

    def plot_energy(self, centers, hist_list, labels, title, xlabel, ylog: bool = False):
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        if not isinstance(hist_list, (list, tuple)):
            hist_list = [hist_list]
        if not isinstance(centers[0], np.ndarray):
            centers_list = [centers]
        else:
            centers_list = centers

        for i, h in enumerate(hist_list):
            if isinstance(centers_list, list) and len(centers_list) == len(hist_list):
                c = centers_list[i]
            else:
                c = centers
            lbl = labels[i] if labels and i < len(labels) else None
            ax.plot(c, h, marker="o", markersize=3, linewidth=1.5, label=lbl)

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Counts")
        ax.set_title(title)
        if ylog and any(np.any(h > 0) for h in hist_list):
            ax.set_yscale("log")
        else:
            ax.set_yscale("linear")
        ax.grid(True, linestyle="--", alpha=0.5)
        if labels and any(labels):
            ax.legend(fontsize=8)
        self.figure.tight_layout()
        self.draw()

    def plot_depth_dose(self, z_list, dose_list, labels, title, xlabel):
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        if not isinstance(z_list, (list, tuple)):
            z_list = [z_list]
        if not isinstance(dose_list, (list, tuple)):
            dose_list = [dose_list]

        for i, zc in enumerate(z_list):
            if zc is None or len(zc) == 0:
                continue
            d = dose_list[i]
            lbl = labels[i] if labels and i < len(labels) else None
            ax.plot(zc, d, marker="o", markersize=3, linewidth=1.5, label=lbl)

        ax.set_xlabel(xlabel)
        ax.set_ylabel("ΣEdep [a.u.]")
        ax.set_title(title)
        ax.grid(True, linestyle="--", alpha=0.5)
        if labels and any(labels):
            ax.legend(fontsize=8)
        self.figure.tight_layout()
        self.draw()

    def plot_attenuation(self, t, I_norm, t_fit, I_fit_norm, mu, HVL, TVL, title):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.scatter(t, I_norm, label="Simülasyon Noktaları", s=20)
        ax.plot(t_fit, I_fit_norm, linestyle="--", linewidth=1.6, label="Fit")
        ax.set_xlabel("Kalınlık [mm]")
        ax.set_ylabel("I / I0")
        ax.set_title(title)
        ax.grid(True, linestyle="--", alpha=0.5)

        text = (
            f"µ = {mu:.4f} 1/mm\n"
            f"HVL = {HVL:.3f} mm\n"
            f"TVL = {TVL:.3f} mm"
        )
        ax.annotate(
            text,
            xy=(0.05, 0.95),
            xycoords="axes fraction",
            va="top",
            ha="left",
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )
        ax.legend(fontsize=8)
        self.figure.tight_layout()
        self.draw()

    def plot_2d_dose(self, H, xedges, zedges, title, xlabel="X", ylabel="Z"):
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        if H.size > 0:
            extent = [xedges[0], xedges[-1], zedges[0], zedges[-1]]
            im = ax.imshow(
                H,
                origin="lower",
                extent=extent,
                aspect="auto",
            )
            self.figure.colorbar(im, ax=ax, label="ΣEdep [a.u.]")

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        self.figure.tight_layout()
        self.draw()


class PlotCanvas3D(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure(figsize=(5, 4), dpi=120)
        super().__init__(fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.updateGeometry()

    def _set_equal_aspect(self, ax, x, y, z):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        if x.size == 0 or y.size == 0 or z.size == 0:
            return
        x_min, x_max = np.min(x), np.max(x)
        y_min, y_max = np.min(y), np.max(y)
        z_min, z_max = np.min(z), np.max(z)
        max_range = max(x_max - x_min, y_max - y_min, z_max - z_min)
        if max_range <= 0:
            max_range = 1.0
        cx = 0.5 * (x_max + x_min)
        cy = 0.5 * (y_max + y_min)
        cz = 0.5 * (z_max + z_min)
        ax.set_xlim(cx - max_range / 2, cx + max_range / 2)
        ax.set_ylim(cy - max_range / 2, cy + max_range / 2)
        ax.set_zlim(cz - max_range / 2, cz + max_range / 2)

    def plot_hits_3d(
        self,
        x,
        y,
        z,
        color=None,
        edep=None,
        title="",
        xlabel="X",
        ylabel="Y",
        zlabel="Z",
        cbar_label="",
        marker_size: int = 5,
    ):
        self.figure.clear()
        ax = self.figure.add_subplot(111, projection="3d")

        if x is not None and y is not None and z is not None and np.size(x) > 0:
            x = np.asarray(x, dtype=float)
            y = np.asarray(y, dtype=float)
            z = np.asarray(z, dtype=float)

            if color is not None:
                c = np.asarray(color, dtype=float)
                finite = np.isfinite(c)
                if np.any(finite):
                    c_min = float(np.min(c[finite]))
                    c_max = float(np.max(c[finite]))
                    if c_max > c_min:
                        c_norm = (c - c_min) / (c_max - c_min)
                    else:
                        c_norm = np.zeros_like(c)
                else:
                    c_norm = np.zeros_like(c)
                sc = ax.scatter(x, y, z, s=marker_size, c=c_norm, depthshade=True)
                label = cbar_label or "Değer"
                self.figure.colorbar(sc, ax=ax, label=label)
            elif edep is not None and np.size(edep) == np.size(x):
                e = np.asarray(edep, dtype=float)
                finite = np.isfinite(e)
                if np.any(finite):
                    e_min = float(np.min(e[finite]))
                    e_max = float(np.max(e[finite]))
                    if e_max > e_min:
                        e_norm = (e - e_min) / (e_max - e_min)
                    else:
                        e_norm = np.zeros_like(e)
                else:
                    e_norm = np.zeros_like(e)
                sc = ax.scatter(x, y, z, s=marker_size, c=e_norm, depthshade=True)
                self.figure.colorbar(sc, ax=ax, label="Normalized Edep")
            else:
                ax.scatter(x, y, z, s=marker_size, depthshade=True)

            self._set_equal_aspect(ax, x, y, z)
        else:
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_zlim(0, 1)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        ax.set_title(title or "3D Hit Haritası")
        self.draw()


class TrackCanvas3D(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure(figsize=(5, 4), dpi=120)
        super().__init__(fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.updateGeometry()

    def _set_equal_aspect(self, ax, x, y, z):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        if x.size == 0 or y.size == 0 or z.size == 0:
            return
        x_min, x_max = np.min(x), np.max(x)
        y_min, y_max = np.min(y), np.max(y)
        z_min, z_max = np.min(z), np.max(z)
        max_range = max(x_max - x_min, y_max - y_min, z_max - z_min)
        if max_range <= 0:
            max_range = 1.0
        cx = 0.5 * (x_max + x_min)
        cy = 0.5 * (y_max + y_min)
        cz = 0.5 * (z_max + z_min)
        ax.set_xlim(cx - max_range / 2, cx + max_range / 2)
        ax.set_ylim(cy - max_range / 2, cy + max_range / 2)
        ax.set_zlim(cz - max_range / 2, cz + max_range / 2)

    def _draw_geometry_boxes(self, ax, volume_bounds: Dict[int, Tuple[float, float, float, float, float, float]]):
        for vid, (xmin, xmax, ymin, ymax, zmin, zmax) in volume_bounds.items():
            corners = [
                (xmin, ymin, zmin),
                (xmax, ymin, zmin),
                (xmax, ymax, zmin),
                (xmin, ymax, zmin),
                (xmin, ymin, zmax),
                (xmax, ymin, zmax),
                (xmax, ymax, zmax),
                (xmin, ymax, zmax),
            ]
            faces = [
                [corners[i] for i in [0, 1, 2, 3]],
                [corners[i] for i in [4, 5, 6, 7]],
                [corners[i] for i in [0, 1, 5, 4]],
                [corners[i] for i in [1, 2, 6, 5]],
                [corners[i] for i in [2, 3, 7, 6]],
                [corners[i] for i in [3, 0, 4, 7]],
            ]
            poly = Poly3DCollection(faces, alpha=0.15)
            ax.add_collection3d(poly)
            cx = 0.5 * (xmin + xmax)
            cy = 0.5 * (ymin + ymax)
            cz = 0.5 * (zmin + zmax)
            ax.text(cx, cy, cz, f"V{vid}", fontsize=8, ha="center", va="center")

    def plot_tracks(self,
                    x: np.ndarray,
                    y: np.ndarray,
                    z: np.ndarray,
                    event_id: np.ndarray,
                    track_id: np.ndarray,
                    parent_id: Optional[np.ndarray],
                    volume_id: Optional[np.ndarray],
                    pdg_id: Optional[np.ndarray],
                    selected_event: int,
                    max_tracks_to_draw: int = 50,
                    only_primary: bool = False,
                    only_secondary: bool = False,
                    pdg_filter: Optional[List[int]] = None,
                    title: str = "Track Viewer"):
        self.figure.clear()
        ax = self.figure.add_subplot(111, projection="3d")

        if x is None or y is None or z is None:
            ax.set_title("Track verisi yok (eventID/trackID/posX/Y/Z bulunamadı).")
            self.draw()
            return

        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        event_id = np.asarray(event_id, dtype=int)
        track_id = np.asarray(track_id, dtype=int)

        if parent_id is not None:
            parent_id = np.asarray(parent_id, dtype=int)
        if volume_id is not None:
            volume_id = np.asarray(volume_id, dtype=int)
        if pdg_id is not None:
            pdg_id = np.asarray(pdg_id, dtype=int)

        mask_ev = (event_id == selected_event)
        if not np.any(mask_ev):
            ax.set_title(f"EventID={selected_event} için adım bulunamadı.")
            self.draw()
            return

        x_ev = x[mask_ev]
        y_ev = y[mask_ev]
        z_ev = z[mask_ev]
        tid_ev = track_id[mask_ev]
        pid_ev = pdg_id[mask_ev] if pdg_id is not None else None
        par_ev = parent_id[mask_ev] if parent_id is not None else None
        vol_ev = volume_id[mask_ev] if volume_id is not None else None

        if pdg_filter and (pid_ev is not None):
            mask_pdg = np.isin(pid_ev, np.array(pdg_filter, dtype=int))
            x_ev, y_ev, z_ev = x_ev[mask_pdg], y_ev[mask_pdg], z_ev[mask_pdg]
            tid_ev = tid_ev[mask_pdg]
            if par_ev is not None:
                par_ev = par_ev[mask_pdg]
            if vol_ev is not None:
                vol_ev = vol_ev[mask_pdg]
            if pid_ev is not None:
                pid_ev = pid_ev[mask_pdg]

        if par_ev is not None:
            if only_primary and not only_secondary:
                mask_prim = (par_ev == 0)
                x_ev, y_ev, z_ev = x_ev[mask_prim], y_ev[mask_prim], z_ev[mask_prim]
                tid_ev = tid_ev[mask_prim]
                par_ev = par_ev[mask_prim]
                if vol_ev is not None:
                    vol_ev = vol_ev[mask_prim]
                if pid_ev is not None:
                    pid_ev = pid_ev[mask_prim]
            elif only_secondary and not only_primary:
                mask_sec = (par_ev != 0)
                x_ev, y_ev, z_ev = x_ev[mask_sec], y_ev[mask_sec], z_ev[mask_sec]
                tid_ev = tid_ev[mask_sec]
                par_ev = par_ev[mask_sec]
                if vol_ev is not None:
                    vol_ev = vol_ev[mask_sec]
                if pid_ev is not None:
                    pid_ev = pid_ev[mask_sec]

        if x_ev.size == 0:
            ax.set_title(f"EventID={selected_event} için filtreden geçen step yok.")
            self.draw()
            return

        unique_tracks = np.unique(tid_ev)
        if unique_tracks.size > max_tracks_to_draw:
            unique_tracks = unique_tracks[:max_tracks_to_draw]

        if vol_ev is not None:
            volume_bounds: Dict[int, Tuple[float, float, float, float, float, float]] = {}
            for vid in np.unique(vol_ev):
                mask_v = (vol_ev == vid)
                xv = x_ev[mask_v]
                yv = y_ev[mask_v]
                zv = z_ev[mask_v]
                if xv.size == 0:
                    continue
                xmin, xmax = float(np.min(xv)), float(np.max(xv))
                ymin, ymax = float(np.min(yv)), float(np.max(yv))
                zmin, zmax = float(np.min(zv)), float(np.max(zv))
                if np.isclose(xmin, xmax):
                    delta = 0.5
                    xmin -= delta
                    xmax += delta
                if np.isclose(ymin, ymax):
                    delta = 0.5
                    ymin -= delta
                    ymax += delta
                if np.isclose(zmin, zmax):
                    delta = 0.5
                    zmin -= delta
                    zmax += delta
                volume_bounds[int(vid)] = (xmin, xmax, ymin, ymax, zmin, zmax)
            self._draw_geometry_boxes(ax, volume_bounds)

        legend_done: Dict[int, bool] = {}
        for tid in unique_tracks:
            m_t = (tid_ev == tid)
            x_t = x_ev[m_t]
            y_t = y_ev[m_t]
            z_t = z_ev[m_t]

            if x_t.size < 2:
                continue

            order = np.argsort(z_t)
            x_t = x_t[order]
            y_t = y_t[order]
            z_t = z_t[order]

            if pid_ev is not None:
                pid_t = pid_ev[m_t][order][0]
                color = pdg_color(pid_t)
                label = None
                if pid_t not in legend_done:
                    label = pdg_name(pid_t)
                    legend_done[pid_t] = True
            else:
                color = None
                label = None
                pid_t = None

            if par_ev is not None:
                par_for_track = par_ev[m_t][order][0]
                if par_for_track == 0:
                    lw = 2.3
                    ls = "-"
                else:
                    lw = 1.2
                    ls = "--"
            else:
                par_for_track = 0
                lw = 1.5
                ls = "-"

            ax.plot(x_t, y_t, z_t, linewidth=lw, linestyle=ls, color=color, label=label)

            if par_for_track != 0:
                x0, y0, z0 = x_t[0], y_t[0], z_t[0]
                ax.scatter(
                    [x0],
                    [y0],
                    [z0],
                    s=40,
                    color=color,
                    edgecolors="k",
                    alpha=0.9,
                )

        if legend_done:
            ax.legend(title="Parçacık türü", fontsize=8)

        self._set_equal_aspect(ax, x_ev, y_ev, z_ev)

        ax.set_xlabel("X [mm]")
        ax.set_ylabel("Y [mm]")
        ax.set_zlabel("Z [mm]")
        ax.set_title(f"{title} (EventID={selected_event}, Tracks={unique_tracks.size})")
        self.draw()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Geant4 Shielding Analyzer GUI (Multi-Plot + Tracks + PDG + Advanced 3D + Event Filter)")
        self.resize(1900, 960)

        self.total_outputs: List[float] = []
        self.thickness_list: List[float] = []
        self.file_names: List[str] = []

        self.particle_dose_summary: Dict[int, float] = {}
        self.volume_dose_summary: Dict[int, Tuple[float, int]] = {}

        self.track_last_x = None
        self.track_last_y = None
        self.track_last_z = None
        self.track_last_event = None
        self.track_last_track = None
        self.track_last_parent = None
        self.track_last_volume = None
        self.track_last_pdg = None

        self.hit3d_last_x = None
        self.hit3d_last_y = None
        self.hit3d_last_z = None
        self.hit3d_last_edep = None
        self.hit3d_last_pdg = None
        self.hit3d_last_volume = None
        self.hit3d_last_event = None
        self.hit3d_last_xname = "X"
        self.hit3d_last_yname = "Y"
        self.hit3d_last_zname = "Z"

        self.energy_results = []
        self.depth_results = []
        self.dose2d_results = []

        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)

        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)

        self.file_list = QListWidget()
        self.file_list.currentItemChanged.connect(self.on_file_selection_changed)

        btn_add = QPushButton("ROOT Dosyası Ekle")
        btn_add.clicked.connect(self.add_files)

        btn_remove = QPushButton("Seçileni Sil")
        btn_remove.clicked.connect(self.remove_selected)

        btn_clear = QPushButton("Listeyi Temizle")
        btn_clear.clicked.connect(self.clear_list)

        thickness_label = QLabel("Seçili dosya için kalınlık [mm]:")
        self.thickness_spin = QDoubleSpinBox()
        self.thickness_spin.setRange(0.0, 1e6)
        self.thickness_spin.setDecimals(3)
        self.thickness_spin.setSingleStep(1.0)
        self.thickness_spin.valueChanged.connect(self.on_thickness_changed)

        self.btn_list_ttrees = QPushButton("Seçili dosyadaki TTree'leri listele")
        self.btn_list_ttrees.clicked.connect(self.list_ttrees)

        self.btn_show_fileinfo = QPushButton("Seçili dosya içeriğini göster")
        self.btn_show_fileinfo.clicked.connect(self.show_file_info)

        left_layout.addWidget(self.file_list)
        left_layout.addWidget(btn_add)

        row_btns = QHBoxLayout()
        row_btns.addWidget(btn_remove)
        row_btns.addWidget(btn_clear)
        left_layout.addLayout(row_btns)

        left_layout.addWidget(self.btn_list_ttrees)
        left_layout.addWidget(self.btn_show_fileinfo)
        left_layout.addWidget(thickness_label)
        left_layout.addWidget(self.thickness_spin)

        middle_widget = QWidget()
        middle_layout = QVBoxLayout(middle_widget)

        form = QFormLayout()
        self.tree_edit = QLineEdit("Hits")
        self.edep_edit = QLineEdit("Edep")
        self.x_edit = QLineEdit("")
        self.y_edit = QLineEdit("")
        self.z_edit = QLineEdit("")
        self.pid_branch_edit = QLineEdit("particleID")
        self.volume_branch_edit = QLineEdit("volumeID")
        self.event_branch_edit = QLineEdit("eventID")
        self.track_branch_edit = QLineEdit("trackID")
        self.parent_branch_edit = QLineEdit("parentID")

        self.nbins_edep_spin = QSpinBox()
        self.nbins_edep_spin.setRange(10, 10000)
        self.nbins_edep_spin.setValue(200)

        self.nbins_depth_spin = QSpinBox()
        self.nbins_depth_spin.setRange(10, 10000)
        self.nbins_depth_spin.setValue(100)

        self.nx_2d_spin = QSpinBox()
        self.nx_2d_spin.setRange(10, 400)
        self.nx_2d_spin.setValue(80)
        self.nz_2d_spin = QSpinBox()
        self.nz_2d_spin.setRange(10, 400)
        self.nz_2d_spin.setValue(80)

        self.energy_title_edit = QLineEdit("Enerji Spektrumu")
        self.depth_title_edit = QLineEdit("Derinlik-Doz Profili")
        self.atten_title_edit = QLineEdit("Attenuation Curve (Total ΣEdep as I)")
        self.dose2d_title_edit = QLineEdit("2D Doz Haritası (X-Z)")
        self.hit3d_title_edit = QLineEdit("3D Hit Haritası")

        form.addRow("TTree adı:", self.tree_edit)
        form.addRow("Edep branch adı:", self.edep_edit)
        form.addRow("X branch adı (opsiyonel):", self.x_edit)
        form.addRow("Y branch adı (opsiyonel):", self.y_edit)
        form.addRow("Z branch adı (opsiyonel):", self.z_edit)
        form.addRow("particleID branch adı:", self.pid_branch_edit)
        form.addRow("volumeID branch adı:", self.volume_branch_edit)
        form.addRow("eventID branch adı:", self.event_branch_edit)
        form.addRow("trackID branch adı:", self.track_branch_edit)
        form.addRow("parentID branch adı:", self.parent_branch_edit)
        form.addRow("Edep bin sayısı:", self.nbins_edep_spin)
        form.addRow("Derinlik-doz bin sayısı:", self.nbins_depth_spin)
        form.addRow("2D doz haritası X bin sayısı:", self.nx_2d_spin)
        form.addRow("2D doz haritası Z bin sayısı:", self.nz_2d_spin)
        form.addRow(QLabel("---- Grafik Başlıkları ----"))
        form.addRow("Enerji grafiği başlığı:", self.energy_title_edit)
        form.addRow("Derinlik-doz başlığı:", self.depth_title_edit)
        form.addRow("Attenuation başlığı:", self.atten_title_edit)
        form.addRow("2D doz başlığı:", self.dose2d_title_edit)
        form.addRow("3D hit başlığı:", self.hit3d_title_edit)

        self.use_zrange_cb = QCheckBox("I hesabı için çıkış yüzeyi Z aralığı kullan")
        self.zmin_out_spin = QDoubleSpinBox()
        self.zmin_out_spin.setRange(-1e6, 1e6)
        self.zmin_out_spin.setDecimals(3)
        self.zmin_out_spin.setSingleStep(0.1)
        self.zmax_out_spin = QDoubleSpinBox()
        self.zmax_out_spin.setRange(-1e6, 1e6)
        self.zmax_out_spin.setDecimals(3)
        self.zmax_out_spin.setSingleStep(0.1)

        form.addRow(self.use_zrange_cb)
        form.addRow("Z_min (çıkış yüzeyi):", self.zmin_out_spin)
        form.addRow("Z_max (çıkış yüzeyi):", self.zmax_out_spin)

        self.use_pid_filter_cb = QCheckBox("Sadece seçilen PDG ID'leri analiz et")
        self.pid_list_edit = QLineEdit("")
        form.addRow(self.use_pid_filter_cb)
        form.addRow("PDG ID listesi (virgülle):", self.pid_list_edit)

        self.use_volume_filter_cb = QCheckBox("Sadece seçilen volumeID'leri analiz et")
        self.volume_list_edit = QLineEdit("")
        form.addRow(self.use_volume_filter_cb)
        form.addRow("Volume ID listesi (virgülle):", self.volume_list_edit)

        self.track_event_spin = QSpinBox()
        self.track_event_spin.setRange(0, 1_000_000)
        self.track_event_spin.setValue(0)

        self.track_only_primary_cb = QCheckBox("Yalnız primary track'ler")
        self.track_only_secondary_cb = QCheckBox("Yalnız secondary track'ler")

        self.track_max_spin = QSpinBox()
        self.track_max_spin.setRange(1, 1000)
        self.track_max_spin.setValue(50)

        self.track_pdg_filter_edit = QLineEdit("")

        form.addRow(QLabel("--- Track Viewer Ayarları ---"))
        form.addRow("Görüntülenecek EventID:", self.track_event_spin)
        form.addRow("Max track sayısı:", self.track_max_spin)
        form.addRow("Track PDG filtresi (virgülle):", self.track_pdg_filter_edit)
        form.addRow(self.track_only_primary_cb)
        form.addRow(self.track_only_secondary_cb)

        self.hit3d_max_points_spin = QSpinBox()
        self.hit3d_max_points_spin.setRange(100, 200000)
        self.hit3d_max_points_spin.setValue(10000)

        self.hit3d_pdg_filter_edit = QLineEdit("")
        self.hit3d_volume_filter_edit = QLineEdit("")

        self.hit3d_color_mode_combo = QComboBox()
        self.hit3d_color_mode_combo.addItems(["Edep", "Z", "PDG", "Volume"])

        self.hit3d_use_event_cb = QCheckBox("3D'de sadece seçilen EventID'deki noktaları göster")
        self.hit3d_point_size_spin = QSpinBox()
        self.hit3d_point_size_spin.setRange(1, 50)
        self.hit3d_point_size_spin.setValue(5)

        form.addRow(QLabel("--- 3D Hit Haritası Ayarları ---"))
        form.addRow("Max nokta sayısı:", self.hit3d_max_points_spin)
        form.addRow("3D PDG filtresi (virgülle):", self.hit3d_pdg_filter_edit)
        form.addRow("3D Volume filtresi (virgülle):", self.hit3d_volume_filter_edit)
        form.addRow("Renk modu:", self.hit3d_color_mode_combo)
        form.addRow(self.hit3d_use_event_cb)
        form.addRow("Nokta boyutu (pixel):", self.hit3d_point_size_spin)

        self.btn_list_branches = QPushButton("Seçili dosyada branch'leri listele")
        self.btn_list_branches.clicked.connect(self.list_branches)

        self.run_button = QPushButton("Analizi Çalıştır")
        self.run_button.clicked.connect(self.run_analysis)

        self.update_tracks_button = QPushButton("Track Viewer'ı EventID'ye göre güncelle")
        self.update_tracks_button.clicked.connect(self.update_track_viewer_from_ui)

        self.update_hit3d_button = QPushButton("3D Hit Haritasını Güncelle")
        self.update_hit3d_button.clicked.connect(self.update_hit3d_from_ui)

        self.save_log_button = QPushButton("Log'u kaydet")
        self.save_log_button.clicked.connect(self.save_log_to_file)

        self.log_edit = QTextEdit()
        self.log_edit.setReadOnly(True)

        middle_layout.addLayout(form)
        middle_layout.addWidget(self.btn_list_branches)
        middle_layout.addWidget(self.run_button)
        middle_layout.addWidget(self.update_tracks_button)
        middle_layout.addWidget(self.update_hit3d_button)
        middle_layout.addWidget(self.save_log_button)
        middle_layout.addWidget(QLabel("Log:"))
        middle_layout.addWidget(self.log_edit)

        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)

        self.tabs = QTabWidget()

        self.energy_tab = QTabWidget()
        self.energy_canvas_single = PlotCanvas()
        self.energy_canvas_multi = PlotCanvas()
        energy_single_widget = QWidget()
        evs_layout = QVBoxLayout(energy_single_widget)
        evs_layout.addWidget(self.energy_canvas_single)
        energy_multi_widget = QWidget()
        evm_layout = QVBoxLayout(energy_multi_widget)
        evm_layout.addWidget(self.energy_canvas_multi)
        self.energy_tab.addTab(energy_single_widget, "Tekil grafik")
        self.energy_tab.addTab(energy_multi_widget, "Kıyaslama (tüm dosyalar)")
        energy_widget = QWidget()
        ev_layout = QVBoxLayout(energy_widget)
        ev_layout.addWidget(self.energy_tab)
        self.tabs.addTab(energy_widget, "Enerji Spektrumu")

        self.depth_tab = QTabWidget()
        self.depth_canvas_single = PlotCanvas()
        self.depth_canvas_multi = PlotCanvas()
        depth_single_widget = QWidget()
        dps_layout = QVBoxLayout(depth_single_widget)
        dps_layout.addWidget(self.depth_canvas_single)
        depth_multi_widget = QWidget()
        dpm_layout = QVBoxLayout(depth_multi_widget)
        dpm_layout.addWidget(self.depth_canvas_multi)
        self.depth_tab.addTab(depth_single_widget, "Tekil grafik")
        self.depth_tab.addTab(depth_multi_widget, "Kıyaslama (tüm dosyalar)")
        depth_widget = QWidget()
        dd_layout = QVBoxLayout(depth_widget)
        dd_layout.addWidget(self.depth_tab)
        self.tabs.addTab(depth_widget, "Derinlik-Doz")

        self.atten_canvas = PlotCanvas()
        atten_widget = QWidget()
        at_layout = QVBoxLayout(atten_widget)
        at_layout.addWidget(self.atten_canvas)
        self.tabs.addTab(atten_widget, "Zayıflama (HVL/TVL)")

        self.hit3d_canvas = PlotCanvas3D()
        hit3d_widget = QWidget()
        h3_layout = QVBoxLayout(hit3d_widget)
        h3_layout.addWidget(self.hit3d_canvas)
        self.tabs.addTab(hit3d_widget, "3D Hit Haritası")

        self.dose2d_tab = QTabWidget()
        self.dose2d_canvas_single = PlotCanvas()
        self.dose2d_canvas_multi = PlotCanvas()
        dose2d_single_widget = QWidget()
        d2s_layout = QVBoxLayout(dose2d_single_widget)
        d2s_layout.addWidget(self.dose2d_canvas_single)
        dose2d_multi_widget = QWidget()
        d2m_layout = QVBoxLayout(dose2d_multi_widget)
        d2m_layout.addWidget(self.dose2d_canvas_multi)
        self.dose2d_tab.addTab(dose2d_single_widget, "Tekil harita (son dosya)")
        self.dose2d_tab.addTab(dose2d_multi_widget, "Birleşik harita (sum)")
        dose2d_widget = QWidget()
        d2_layout = QVBoxLayout(dose2d_widget)
        d2_layout.addWidget(self.dose2d_tab)
        self.tabs.addTab(dose2d_widget, "2D Doz Haritası (X-Z)")

        self.summary_table = QTableWidget(0, 4)
        self.summary_table.setHorizontalHeaderLabels(["Dosya", "Kalınlık [mm]", "ΣEdep", "I/I0"])
        summary_widget = QWidget()
        s_layout = QVBoxLayout(summary_widget)
        s_layout.addWidget(self.summary_table)
        self.tabs.addTab(summary_widget, "Özet")

        self.particle_table = QTableWidget(0, 4)
        self.particle_table.setHorizontalHeaderLabels(["PDG ID", "İsim", "ΣEdep", "Toplam İçindeki Payı"])
        particle_widget = QWidget()
        p_layout = QVBoxLayout(particle_widget)
        p_layout.addWidget(self.particle_table)
        self.tabs.addTab(particle_widget, "Parçacık Özeti")

        self.volume_table = QTableWidget(0, 4)
        self.volume_table.setHorizontalHeaderLabels(
            ["Volume ID", "ΣEdep", "Hit Sayısı", "Toplam ΣEdep İçindeki Payı"]
        )
        volume_widget = QWidget()
        v_layout = QVBoxLayout(volume_widget)
        v_layout.addWidget(self.volume_table)
        self.tabs.addTab(volume_widget, "Volume Özeti")

        self.track_canvas = TrackCanvas3D()
        track_widget = QWidget()
        tv_layout = QVBoxLayout(track_widget)
        tv_layout.addWidget(self.track_canvas)
        self.tabs.addTab(track_widget, "Track Viewer")

        self.energy_log_cb = QCheckBox("Enerji spektrumu için log Y ekseni")
        self.energy_log_cb.stateChanged.connect(self.redraw_energy_plots)

        right_layout.addWidget(self.energy_log_cb)
        right_layout.addWidget(self.tabs)

        main_layout.addWidget(left_widget, 2)
        main_layout.addWidget(middle_widget, 3)
        main_layout.addWidget(right_widget, 7)

        self._last_energy_single = None
        self._last_depth_single = None
        self._last_2d_single = None

    def log(self, msg: str):
        self.log_edit.append(msg)
        self.log_edit.moveCursor(self.log_edit.textCursor().End)
        QApplication.processEvents()

    def current_item(self):
        return self.file_list.currentItem()

    def save_log_to_file(self):
        text = self.log_edit.toPlainText()
        if not text.strip():
            QMessageBox.information(self, "Bilgi", "Log boş, kaydedilecek içerik yok.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Log'u kaydet",
            "shielding_log.txt",
            "Text files (*.txt);;All files (*)",
        )
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(text)
            QMessageBox.information(self, "Bilgi", f"Log kaydedildi:\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Hata", f"Log kaydedilemedi:\n{e}")

    def add_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "ROOT Dosyası Seç",
            "",
            "ROOT files (*.root);;All files (*)",
        )
        if not files:
            return
        for f in files:
            basename = os.path.basename(f)
            inferred_thickness = infer_thickness_from_filename(basename)
            item = QListWidgetItem(basename)
            item.setData(Qt.UserRole, {"path": f, "thickness_mm": inferred_thickness or 0.0})
            self.file_list.addItem(item)
            if inferred_thickness is not None:
                self.log(f"{basename} için dosya adından kalınlık tahmin edildi: {inferred_thickness} mm")
        self.log(f"{len(files)} dosya eklendi.")

    def remove_selected(self):
        item = self.current_item()
        if item is None:
            return
        row = self.file_list.row(item)
        self.file_list.takeItem(row)
        self.log("Seçili dosya listeden silindi.")

    def clear_list(self):
        self.file_list.clear()
        self.log("Dosya listesi temizlendi.")

    def on_file_selection_changed(self, current, previous):
        if current is None:
            self.thickness_spin.blockSignals(True)
            self.thickness_spin.setValue(0.0)
            self.thickness_spin.blockSignals(False)
            return
        data = current.data(Qt.UserRole) or {}
        thickness = data.get("thickness_mm", 0.0) or 0.0
        self.thickness_spin.blockSignals(True)
        self.thickness_spin.setValue(thickness)
        self.thickness_spin.blockSignals(False)

    def on_thickness_changed(self, value: float):
        item = self.current_item()
        if item is None:
            return
        data = item.data(Qt.UserRole) or {}
        data["thickness_mm"] = float(value)
        item.setData(Qt.UserRole, data)

    def list_ttrees(self):
        item = self.current_item()
        if item is None:
            QMessageBox.warning(self, "Uyarı", "Önce sol listeden bir dosya seç.")
            return
        data = item.data(Qt.UserRole) or {}
        path = data.get("path")
        try:
            names = list_ttrees_in_file(path)
        except Exception as e:
            self.log(f"[HATA] TTree listelenemedi: {e}")
            QMessageBox.critical(self, "Hata", str(e))
            return

        if not names:
            text = "Bu dosyada TTree bulunamadı."
        else:
            text = "Bu dosyadaki TTree'ler:\n" + "\n".join(f"  - {n}" for n in names)
        self.log(text)
        QMessageBox.information(self, "TTrees", text)

    def list_branches(self):
        item = self.current_item()
        if item is None:
            QMessageBox.warning(self, "Uyarı", "Önce sol listeden bir dosya seç.")
            return
        tree_name = self.tree_edit.text().strip()
        if not tree_name:
            QMessageBox.warning(self, "Uyarı", "Önce TTree adını gir.")
            return
        data = item.data(Qt.UserRole) or {}
        path = data.get("path")
        try:
            branches = list_branches_in_tree(path, tree_name)
        except Exception as e:
            self.log(f"[HATA] Branch listelenemedi: {e}")
            QMessageBox.critical(self, "Hata", str(e))
            return

        if not branches:
            text = f"'{tree_name}' içinde branch bulunamadı."
        else:
            text = f"'{tree_name}' içindeki branch'ler:\n" + "\n".join(f"  - {b}" for b in branches)

        self.log(text)
        QMessageBox.information(self, "Branch'ler", text)

    def show_file_info(self):
        item = self.current_item()
        if item is None:
            QMessageBox.warning(self, "Uyarı", "Önce sol listeden bir dosya seç.")
            return
        data = item.data(Qt.UserRole) or {}
        path = data.get("path")

        try:
            f = uproot.open(path)
        except Exception as e:
            self.log(f"[HATA] Dosya açılamadı: {e}")
            QMessageBox.critical(self, "Hata", str(e))
            return

        lines = [f"Dosya: {path}", "İçerik:"]
        for key, obj in f.items():
            name = key.split(";")[0]
            tname = type(obj).__name__
            extra = ""
            if isinstance(obj, uproot.behaviors.TTree.TTree):
                extra = f" (TTree, entries={obj.num_entries})"
            lines.append(f"  - {name}: {tname}{extra}")

        text = "\n".join(lines)
        self.log(text)
        QMessageBox.information(self, "Dosya İçeriği", text)

    def run_analysis(self):
        n_items = self.file_list.count()
        if n_items == 0:
            QMessageBox.warning(self, "Uyarı", "Önce en az bir ROOT dosyası eklemelisin.")
            return

        tree_name = self.tree_edit.text().strip()
        edep_branch = self.edep_edit.text().strip()
        x_branch_user = self.x_edit.text().strip() or None
        y_branch_user = self.y_edit.text().strip() or None
        z_branch_user = self.z_edit.text().strip() or None
        pid_branch = self.pid_branch_edit.text().strip() or None
        volume_branch = self.volume_branch_edit.text().strip() or None
        event_branch = self.event_branch_edit.text().strip() or None
        track_branch = self.track_branch_edit.text().strip() or None
        parent_branch = self.parent_branch_edit.text().strip() or None

        nbins_edep = self.nbins_edep_spin.value()
        nbins_depth = self.nbins_depth_spin.value()
        nx_2d = self.nx_2d_spin.value()
        nz_2d = self.nz_2d_spin.value()

        if not tree_name:
            QMessageBox.warning(self, "Uyarı", "TTree adı boş olamaz.")
            return
        if not edep_branch:
            QMessageBox.warning(self, "Uyarı", "Edep branch adı boş olamaz.")
            return

        use_zwindow = self.use_zrange_cb.isChecked()
        zmin_out = self.zmin_out_spin.value()
        zmax_out = self.zmax_out_spin.value()

        use_pid_filter = self.use_pid_filter_cb.isChecked()
        pdg_list = parse_int_list(self.pid_list_edit.text())
        if use_pid_filter and (pid_branch is None or not pdg_list):
            self.log("[UYARI] Parçacık filtresi aktif ama particleID branch veya PDG listesi tanımlı değil. Filtre uygulanmayacak.")
            use_pid_filter = False

        use_volume_filter = self.use_volume_filter_cb.isChecked()
        volume_list = parse_int_list(self.volume_list_edit.text())
        if use_volume_filter and (volume_branch is None or not volume_list):
            self.log("[UYARI] Volume filtresi aktif ama volumeID branch veya volume listesi tanımlı değil. Filtre uygulanmayacak.")
            use_volume_filter = False

        self.log("=== Analiz başlıyor ===")
        self.log(
            f"TTree: {tree_name}, Edep: {edep_branch}, "
            f"X: {x_branch_user}, Y: {y_branch_user}, Z: {z_branch_user}, "
            f"PID: {pid_branch}, Volume: {volume_branch}, "
            f"eventID: {event_branch}, trackID: {track_branch}, parentID: {parent_branch}"
        )

        if use_zwindow and z_branch_user is not None and zmax_out > zmin_out:
            self.log(f"I hesabı için Z aralığı kullanılacak: [{zmin_out},{zmax_out}]")
        else:
            if use_zwindow and z_branch_user is None:
                self.log("[UYARI] Z branch tanımlı değil, Z aralığı kullanılamaz. Tüm hacim kullanılacak.")
            elif use_zwindow and not (zmax_out > zmin_out):
                self.log("[UYARI] Z_min < Z_max değil, Z aralığı geçersiz. Tüm hacim kullanılacak.")
            self.log("I hesabı tüm Edep üzerinden yapılacak (Z bağımsız).")

        if use_pid_filter:
            self.log(f"Parçacık filtresi aktif. PDG ID'ler: {pdg_list}")
        if use_volume_filter:
            self.log(f"Volume filtresi aktif. Volume ID'ler: {volume_list}")

        self.total_outputs = []
        self.thickness_list = []
        self.file_names = []
        self.particle_dose_summary = {}
        self.volume_dose_summary = {}
        self.summary_table.setRowCount(0)
        self.particle_table.setRowCount(0)
        self.volume_table.setRowCount(0)

        self.track_last_x = None
        self.track_last_y = None
        self.track_last_z = None
        self.track_last_event = None
        self.track_last_track = None
        self.track_last_parent = None
        self.track_last_volume = None
        self.track_last_pdg = None

        self.hit3d_last_x = None
        self.hit3d_last_y = None
        self.hit3d_last_z = None
        self.hit3d_last_edep = None
        self.hit3d_last_pdg = None
        self.hit3d_last_volume = None
        self.hit3d_last_event = None
        self.hit3d_last_xname = x_branch_user or "X"
        self.hit3d_last_yname = y_branch_user or "Y"
        self.hit3d_last_zname = z_branch_user or "Z"

        self.energy_results.clear()
        self.depth_results.clear()
        self.dose2d_results.clear()

        out_dir = os.getcwd()

        self._last_energy_single = None
        self._last_depth_single = None
        self._last_2d_single = None

        for idx in range(n_items):
            item = self.file_list.item(idx)
            data = item.data(Qt.UserRole) or {}
            path = data.get("path")
            thickness_mm = float(data.get("thickness_mm", 0.0) or 0.0)
            base = os.path.splitext(os.path.basename(path))[0]

            if thickness_mm == 0.0:
                inferred = infer_thickness_from_filename(base)
                if inferred is not None:
                    thickness_mm = inferred
                    data["thickness_mm"] = inferred
                    item.setData(Qt.UserRole, data)
                    self.log(f"{base} için kalınlık isimden otomatik alındı: {inferred} mm")

            self.log(f"[{idx+1}/{n_items}] Dosya işleniyor: {path}")

            try:
                tree = read_tree(path, tree_name)
            except Exception as e:
                self.log(f"[HATA] TTree okunamadı: {e}")
                QMessageBox.critical(self, "Hata", str(e))
                continue

            try:
                all_branches = list(tree.keys())
            except Exception:
                all_branches = []

            try:
                edep_raw = get_branch(tree, edep_branch)
            except Exception as e:
                self.log(f"[HATA] {e}")
                QMessageBox.critical(self, "Hata", str(e))
                continue

            x_branch = x_branch_user
            y_branch = y_branch_user
            z_branch = z_branch_user
            if x_branch is None or y_branch is None or z_branch is None:
                if all_branches:
                    ax_b, ay_b, az_b = autodetect_xyz_branches(all_branches)
                    if x_branch is None and ax_b is not None:
                        x_branch = ax_b
                        self.log(f"  -> X branch auto-detect: {x_branch}")
                    if y_branch is None and ay_b is not None:
                        y_branch = ay_b
                        self.log(f"  -> Y branch auto-detect: {y_branch}")
                    if z_branch is None and az_b is not None:
                        z_branch = az_b
                        self.log(f"  -> Z branch auto-detect: {z_branch}")

            x_arr_raw = y_arr_raw = z_arr_raw = None
            if x_branch is not None:
                try:
                    x_arr_raw = get_branch(tree, x_branch)
                except Exception as e:
                    self.log(f"[UYARI] X branch okunamadı ({x_branch}): {e}")
            if y_branch is not None:
                try:
                    y_arr_raw = get_branch(tree, y_branch)
                except Exception as e:
                    self.log(f"[UYARI] Y branch okunamadı ({y_branch}): {e}")
            if z_branch is not None:
                try:
                    z_arr_raw = get_branch(tree, z_branch)
                except Exception as e:
                    self.log(f"[UYARI] Z branch okunamadı ({z_branch}): {e}")

            pid_arr = None
            if pid_branch is not None:
                try:
                    pid_arr = get_branch(tree, pid_branch)
                except Exception as e:
                    self.log(f"[UYARI] particleID branch okunamadı ({pid_branch}): {e}")
                    pid_arr = None

            volume_arr = None
            if volume_branch is not None:
                try:
                    volume_arr = get_branch(tree, volume_branch)
                except Exception as e:
                    self.log(f"[UYARI] volumeID branch okunamadı ({volume_branch}): {e}")
                    volume_arr = None

            event_arr = None
            if event_branch is not None:
                try:
                    event_arr = get_branch(tree, event_branch)
                except Exception as e:
                    self.log(f"[UYARI] eventID branch okunamadı ({event_branch}): {e}")
                    event_arr = None

            track_arr = None
            if track_branch is not None:
                try:
                    track_arr = get_branch(tree, track_branch)
                except Exception as e:
                    self.log(f"[UYARI] trackID branch okunamadı ({track_branch}): {e}")
                    track_arr = None

            parent_arr = None
            if parent_branch is not None:
                try:
                    parent_arr = get_branch(tree, parent_branch)
                except Exception as e:
                    self.log(f"[UYARI] parentID branch okunamadı ({parent_branch}): {e}")
                    parent_arr = None

            edep = np.asarray(edep_raw, dtype=float)
            x_arr = np.asarray(x_arr_raw, dtype=float) if x_arr_raw is not None else None
            y_arr = np.asarray(y_arr_raw, dtype=float) if y_arr_raw is not None else None
            z_arr = np.asarray(z_arr_raw, dtype=float) if z_arr_raw is not None else None

            if pid_arr is not None:
                pid_arr = np.asarray(pid_arr, dtype=int)
            if volume_arr is not None:
                volume_arr = np.asarray(volume_arr, dtype=int)
            if event_arr is not None:
                event_arr = np.asarray(event_arr, dtype=int)
            if track_arr is not None:
                track_arr = np.asarray(track_arr, dtype=int)
            if parent_arr is not None:
                parent_arr = np.asarray(parent_arr, dtype=int)

            if use_pid_filter and (pid_arr is not None):
                mask = np.isin(pid_arr, np.array(pdg_list, dtype=int))
                if x_arr is not None:
                    x_arr = x_arr[mask]
                if y_arr is not None:
                    y_arr = y_arr[mask]
                if z_arr is not None:
                    z_arr = z_arr[mask]
                if volume_arr is not None:
                    volume_arr = volume_arr[mask]
                if event_arr is not None:
                    event_arr = event_arr[mask]
                if track_arr is not None:
                    track_arr = track_arr[mask]
                if parent_arr is not None:
                    parent_arr = parent_arr[mask]
                edep = edep[mask]
                pid_arr = pid_arr[mask]
                self.log(f"  -> Parçacık filtresi sonrası event sayısı: {np.count_nonzero(mask)}")

            if use_volume_filter and (volume_arr is not None):
                mask_v = np.isin(volume_arr, np.array(volume_list, dtype=int))
                if x_arr is not None:
                    x_arr = x_arr[mask_v]
                if y_arr is not None:
                    y_arr = y_arr[mask_v]
                if z_arr is not None:
                    z_arr = z_arr[mask_v]
                if pid_arr is not None:
                    pid_arr = pid_arr[mask_v]
                if event_arr is not None:
                    event_arr = event_arr[mask_v]
                if track_arr is not None:
                    track_arr = track_arr[mask_v]
                if parent_arr is not None:
                    parent_arr = parent_arr[mask_v]
                edep = edep[mask_v]
                volume_arr = volume_arr[mask_v]
                self.log(f"  -> Volume filtresi sonrası event sayısı: {np.count_nonzero(mask_v)}")

            centers, hist = make_energy_spectrum(edep, nbins=nbins_edep)
            self.energy_results.append(
                {"file": base, "thickness": thickness_mm, "centers": centers, "hist": hist}
            )
            self._last_energy_single = (centers, hist, base)

            energy_title = self.energy_title_edit.text().strip() or "Enerji Spektrumu"
            xlabel_e = f"{edep_branch} [a.u.]"
            self.energy_canvas_single.plot_energy(
                centers, [hist], [f"{base} ({thickness_mm:.1f} mm)"],
                title=f"{energy_title} - {base}",
                xlabel=xlabel_e,
                ylog=self.energy_log_cb.isChecked(),
            )

            png_edep = os.path.join(out_dir, f"{base}_EdepSpectrum.png")
            csv_edep = os.path.join(out_dir, f"{base}_EdepSpectrum.csv")
            try:
                self.energy_canvas_single.figure.savefig(png_edep, dpi=300)
                data_edep = np.vstack([centers, hist]).T
                save_csv(csv_edep, "E_center,Counts", data_edep)
                self.log(f"  -> Enerji spektrumu kaydedildi: {png_edep}, {csv_edep}")
            except Exception as e:
                self.log(f"[UYARI] Enerji spektrumu kaydedilemedi: {e}")

            z_centers = dose = None
            if z_arr is not None:
                z_centers, dose = make_depth_dose(z_arr, edep, nbins=nbins_depth)
                self.depth_results.append(
                    {"file": base, "thickness": thickness_mm, "z_centers": z_centers, "dose": dose}
                )
                self._last_depth_single = (z_centers, dose, base)
                depth_title = self.depth_title_edit.text().strip() or "Derinlik-Doz Profili"
                self.depth_canvas_single.plot_depth_dose(
                    [z_centers], [dose],
                    [f"{base} ({thickness_mm:.1f} mm)"],
                    title=f"{depth_title} - {base}",
                    xlabel=z_branch or "Z",
                )
                png_dd = os.path.join(out_dir, f"{base}_DepthDose.png")
                csv_dd = os.path.join(out_dir, f"{base}_DepthDose.csv")
                try:
                    self.depth_canvas_single.figure.savefig(png_dd, dpi=300)
                    data_dd = np.vstack([z_centers, dose]).T
                    save_csv(csv_dd, "Z_center,Dose_sumEdep", data_dd)
                    self.log(f"  -> Derinlik-doz kaydedildi: {png_dd}, {csv_dd}")
                except Exception as e:
                    self.log(f"[UYARI] Derinlik-doz kaydedilemedi: {e}")
            else:
                self.log("  -> Z branch yok, derinlik-doz grafiği çizilemedi.")

            if (x_arr is not None) and (z_arr is not None):
                H2d, xedges, zedges = make_2d_dose_map(x_arr, z_arr, edep, nx=nx_2d, nz=nz_2d)
                self.dose2d_results.append(
                    {"file": base, "thickness": thickness_mm, "H": H2d, "xedges": xedges, "zedges": zedges}
                )
                self._last_2d_single = (H2d, xedges, zedges, base)
                dose2d_title = self.dose2d_title_edit.text().strip() or "2D Doz Haritası (X-Z)"
                self.dose2d_canvas_single.plot_2d_dose(
                    H2d,
                    xedges,
                    zedges,
                    title=f"{dose2d_title} - {base}",
                    xlabel=x_branch or "X",
                    ylabel=z_branch or "Z",
                )
                png_2d = os.path.join(out_dir, f"{base}_DoseMap_XZ.png")
                csv_2d = os.path.join(out_dir, f"{base}_DoseMap_XZ.csv")
                try:
                    self.dose2d_canvas_single.figure.savefig(png_2d, dpi=300)
                    with open(csv_2d, "w", encoding="utf-8") as fcsv:
                        fcsv.write("# X edges:\n")
                        fcsv.write(",".join(str(v) for v in xedges) + "\n")
                        fcsv.write("# Z edges:\n")
                        fcsv.write(",".join(str(v) for v in zedges) + "\n")
                        fcsv.write("# Dose matrix (rows: Z bins, cols: X bins):\n")
                        np.savetxt(fcsv, H2d, delimiter=",")
                    self.log(f"  -> 2D doz haritası kaydedildi: {png_2d}, {csv_2d}")
                except Exception as e:
                    self.log(f"[UYARI] 2D doz haritası kaydedilemedi: {e}")
            else:
                self.log("  -> 2D doz haritası için X veya Z eksik, 2D grafik çizilmedi.")

            if use_zwindow and (z_arr is not None) and (zmax_out > zmin_out):
                total_I = compute_total_output_from_window(edep, z_arr, zmin_out, zmax_out)
                self.log(f"  -> ΣEdep (Z penceresi [{zmin_out},{zmax_out}]) = {total_I:.6e}")
            else:
                total_I = compute_total_output_from_window(edep, None, None, None)
                self.log(f"  -> ΣEdep (tüm hacim) = {total_I:.6e}")

            self.total_outputs.append(total_I)
            self.thickness_list.append(thickness_mm)
            self.file_names.append(base)
            self.log(f"  -> Kalınlık = {thickness_mm:.3f} mm")

            if pid_arr is not None:
                for pdg in np.unique(pid_arr):
                    mask_p = (pid_arr == pdg)
                    dose_p = float(np.sum(edep[mask_p]))
                    self.particle_dose_summary[pdg] = self.particle_dose_summary.get(pdg, 0.0) + dose_p

            if volume_arr is not None:
                for vid in np.unique(volume_arr):
                    mask_v2 = (volume_arr == vid)
                    dose_v = float(np.sum(edep[mask_v2]))
                    count_v = int(np.count_nonzero(mask_v2))
                    total_dose, total_count = self.volume_dose_summary.get(vid, (0.0, 0))
                    self.volume_dose_summary[vid] = (total_dose + dose_v, total_count + count_v)

            if (event_arr is not None) and (track_arr is not None) and (x_arr is not None) and (y_arr is not None) and (z_arr is not None):
                self.track_last_x = x_arr
                self.track_last_y = y_arr
                self.track_last_z = z_arr
                self.track_last_event = event_arr
                self.track_last_track = track_arr
                self.track_last_parent = parent_arr
                self.track_last_volume = volume_arr
                self.track_last_pdg = pid_arr
                try:
                    ev_min = int(np.min(event_arr))
                    self.track_event_spin.setValue(ev_min)
                    self.log(f"  -> Track viewer için eventID aralığı: [{int(np.min(event_arr))}, {int(np.max(event_arr))}]")
                except Exception:
                    pass

            if (x_arr is not None) and (y_arr is not None) and (z_arr is not None):
                self.hit3d_last_x = x_arr
                self.hit3d_last_y = y_arr
                self.hit3d_last_z = z_arr
                self.hit3d_last_edep = edep
                self.hit3d_last_pdg = pid_arr
                self.hit3d_last_volume = volume_arr
                self.hit3d_last_event = event_arr
                self.hit3d_last_xname = x_branch or "X"
                self.hit3d_last_yname = y_branch or "Y"
                self.hit3d_last_zname = z_branch or "Z"
            else:
                self.log("  -> 3D hit haritası için X/Y/Z eksik, 3D grafik çizimi atlandı.")

            QApplication.processEvents()

        self.update_energy_multi_plot(edep_branch)
        self.update_depth_multi_plot(z_branch_user or z_branch)
        self.update_dose2d_multi_plot(x_branch_user or x_branch, z_branch_user or z_branch)

        self.update_summary_table(None)
        self.update_particle_summary_table()
        self.update_volume_summary_table()

        self.run_attenuation_analysis(out_dir)

        self.update_track_viewer_from_ui()
        self.update_hit3d_from_ui()

        self.log("=== Analiz tamamlandı ===")
        self.log("Not: 3D veya Track grafiğini göremiyorsan log'daki 'X/Y/Z eksik' mesajına bakıp branch isimlerini düzelt.")

    def redraw_energy_plots(self):
        if self._last_energy_single is not None:
            centers, hist, base = self._last_energy_single
            title = self.energy_title_edit.text().strip() or "Enerji Spektrumu"
            xlabel_e = f"{self.edep_edit.text().strip() or 'Edep'} [a.u.]"
            self.energy_canvas_single.plot_energy(
                centers, [hist], [f"{base}"],
                title=f"{title} - {base}",
                xlabel=xlabel_e,
                ylog=self.energy_log_cb.isChecked(),
            )
        self.update_energy_multi_plot(self.edep_edit.text().strip() or "Edep")

    def update_energy_multi_plot(self, edep_branch: str):
        if not self.energy_results:
            self.energy_canvas_multi.figure.clear()
            self.energy_canvas_multi.draw()
            return

        x_min = min(res["centers"][0] for res in self.energy_results if res["centers"].size > 0)
        x_max = max(res["centers"][-1] for res in self.energy_results if res["centers"].size > 0)
        if not np.isfinite(x_min) or not np.isfinite(x_max) or x_max <= x_min:
            x_min, x_max = 0.0, 1.0
        x_common = np.linspace(x_min, x_max, 400)

        hist_list = []
        labels = []
        centers_list = []
        for res in self.energy_results:
            c = res["centers"]
            h = res["hist"]
            if c.size == 0:
                continue
            y_interp = np.interp(x_common, c, h)
            hist_list.append(y_interp)
            centers_list.append(x_common)
            labels.append(f"{res['file']} ({res['thickness']:.1f} mm)")

        if not hist_list:
            self.energy_canvas_multi.figure.clear()
            self.energy_canvas_multi.draw()
            return

        title = self.energy_title_edit.text().strip() or "Enerji Spektrumu"
        xlabel_e = f"{edep_branch} [a.u.]"
        self.energy_canvas_multi.plot_energy(
            centers_list,
            hist_list,
            labels,
            title=f"{title} - Kıyaslama",
            xlabel=xlabel_e,
            ylog=self.energy_log_cb.isChecked(),
        )

    def update_depth_multi_plot(self, z_branch: Optional[str]):
        if not self.depth_results:
            self.depth_canvas_multi.figure.clear()
            self.depth_canvas_multi.draw()
            return

        z_list = []
        dose_list = []
        labels = []
        for res in self.depth_results:
            zc = res["z_centers"]
            d = res["dose"]
            if zc is None or len(zc) == 0:
                continue
            z_list.append(zc)
            dose_list.append(d)
            labels.append(f"{res['file']} ({res['thickness']:.1f} mm)")

        if not z_list:
            self.depth_canvas_multi.figure.clear()
            self.depth_canvas_multi.draw()
            return

        title = self.depth_title_edit.text().strip() or "Derinlik-Doz Profili"
        xlabel = z_branch or "Z"
        self.depth_canvas_multi.plot_depth_dose(
            z_list,
            dose_list,
            labels,
            title=f"{title} - Kıyaslama",
            xlabel=xlabel,
        )

    def update_dose2d_multi_plot(self, x_branch: Optional[str], z_branch: Optional[str]):
        if not self.dose2d_results:
            self.dose2d_canvas_multi.figure.clear()
            self.dose2d_canvas_multi.draw()
            return

        H_acc = None
        xedges = None
        zedges = None
        for res in self.dose2d_results:
            H = res["H"]
            xe = res["xedges"]
            ze = res["zedges"]
            if H_acc is None:
                H_acc = H.copy()
                xedges = xe
                zedges = ze
            else:
                if H.shape == H_acc.shape:
                    H_acc += H
                else:
                    continue

        if H_acc is None or xedges is None or zedges is None:
            self.dose2d_canvas_multi.figure.clear()
            self.dose2d_canvas_multi.draw()
            return

        dose2d_title = self.dose2d_title_edit.text().strip() or "2D Doz Haritası (X-Z)"
        self.dose2d_canvas_multi.plot_2d_dose(
            H_acc,
            xedges,
            zedges,
            title=f"{dose2d_title} - Birleşik (sum)",
            xlabel=x_branch or "X",
            ylabel=z_branch or "Z",
        )

    def update_summary_table(self, I_norm_full: Optional[np.ndarray]):
        n = len(self.file_names)
        self.summary_table.setRowCount(n)
        for i in range(n):
            name = self.file_names[i] if i < len(self.file_names) else ""
            thickness = self.thickness_list[i] if i < len(self.thickness_list) else 0.0
            totalI = self.total_outputs[i] if i < len(self.total_outputs) else 0.0

            self.summary_table.setItem(i, 0, QTableWidgetItem(str(name)))
            self.summary_table.setItem(i, 1, QTableWidgetItem(f"{thickness:.3f}"))
            self.summary_table.setItem(i, 2, QTableWidgetItem(f"{totalI:.6e}"))

            if I_norm_full is not None and i < len(I_norm_full) and np.isfinite(I_norm_full[i]):
                self.summary_table.setItem(i, 3, QTableWidgetItem(f"{I_norm_full[i]:.6f}"))
            else:
                self.summary_table.setItem(i, 3, QTableWidgetItem("N/A"))

        self.summary_table.resizeColumnsToContents()

    def update_particle_summary_table(self):
        summary = self.particle_dose_summary
        if not summary:
            self.particle_table.setRowCount(0)
            return

        pdgs = sorted(summary.keys())
        total_all = float(sum(summary.values())) if summary else 0.0

        self.particle_table.setRowCount(len(pdgs))
        for i, pdg in enumerate(pdgs):
            dose = summary[pdg]
            frac = dose / total_all if total_all > 0.0 else 0.0
            self.particle_table.setItem(i, 0, QTableWidgetItem(str(pdg)))
            self.particle_table.setItem(i, 1, QTableWidgetItem(pdg_name(pdg)))
            self.particle_table.setItem(i, 2, QTableWidgetItem(f"{dose:.6e}"))
            self.particle_table.setItem(i, 3, QTableWidgetItem(f"{frac:.4f}"))
        self.particle_table.resizeColumnsToContents()

    def update_volume_summary_table(self):
        summary = self.volume_dose_summary
        if not summary:
            self.volume_table.setRowCount(0)
            return

        vids = sorted(summary.keys())
        total_all = float(sum(v[0] for v in summary.values())) if summary else 0.0

        self.volume_table.setRowCount(len(vids))
        for i, vid in enumerate(vids):
            dose, count = summary[vid]
            frac = dose / total_all if total_all > 0.0 else 0.0
            self.volume_table.setItem(i, 0, QTableWidgetItem(str(vid)))
            self.volume_table.setItem(i, 1, QTableWidgetItem(f"{dose:.6e}"))
            self.volume_table.setItem(i, 2, QTableWidgetItem(str(count)))
            self.volume_table.setItem(i, 3, QTableWidgetItem(f"{frac:.4f}"))
        self.volume_table.resizeColumnsToContents()

    def run_attenuation_analysis(self, out_dir: str):
        if not self.thickness_list or not self.total_outputs:
            self.log("[BİLGİ] HVL/TVL için veri yok.")
            self.update_summary_table(None)
            return

        t_full = np.asarray(self.thickness_list, dtype=float)
        I_full = np.asarray(self.total_outputs, dtype=float)

        try:
            mu, HVL, TVL, fit_curve, I0, mask = fit_attenuation(t_full.tolist(), I_full.tolist())
        except Exception as e:
            self.log(f"[BİLGİ] HVL/TVL hesabı atlandı: {e}")
            self.update_summary_table(None)
            return

        t_valid = t_full[mask]
        I_valid = I_full[mask]
        I_norm_valid = I_valid / I0

        I_norm_full = np.full_like(I_full, np.nan, dtype=float)
        I_norm_full[mask] = I_norm_valid

        t_fit = fit_curve[0]
        I_fit = fit_curve[1]
        I_fit_norm = I_fit / I0

        atten_title = self.atten_title_edit.text().strip() or "Attenuation Curve (Total ΣEdep as I)"
        self.atten_canvas.plot_attenuation(
            t_valid,
            I_norm_valid,
            t_fit,
            I_fit_norm,
            mu,
            HVL,
            TVL,
            title=atten_title,
        )

        png_att = os.path.join(out_dir, "AttenuationCurve.png")
        csv_att = os.path.join(out_dir, "AttenuationData.csv")
        try:
            self.atten_canvas.figure.savefig(png_att, dpi=300)
            data_att = np.vstack([t_valid, I_valid, I_norm_valid]).T
            header = (
                "Thickness_mm,TotalEdep,I_over_I0\n"
                f"# mu[1/mm]={mu:.6f}, HVL[mm]={HVL:.6f}, TVL[mm]={TVL:.6f}"
            )
            save_csv(csv_att, header, data_att)
            self.log("=== ZAYIFLAMA SONUÇLARI ===")
            self.log(f"µ   = {mu:.6f} 1/mm")
            self.log(f"HVL = {HVL:.6f} mm")
            self.log(f"TVL = {TVL:.6f} mm")
            self.log(f"AttenuationCurve.png ve AttenuationData.csv kaydedildi.")
        except Exception as e:
            self.log(f"[UYARI] Zayıflama verileri kaydedilemedi: {e}")

        self.update_summary_table(I_norm_full)

    def update_track_viewer_from_ui(self):
        if self.track_last_x is None or self.track_last_event is None or self.track_last_track is None:
            self.track_canvas.plot_tracks(
                None, None, None,
                None, None, None, None, None,
                selected_event=0,
                title="Track Viewer (track verisi yok)",
            )
            return

        selected_event = int(self.track_event_spin.value())
        max_tracks = int(self.track_max_spin.value())
        only_prim = self.track_only_primary_cb.isChecked()
        only_sec = self.track_only_secondary_cb.isChecked()
        pdg_filter = parse_int_list(self.track_pdg_filter_edit.text())

        self.track_canvas.plot_tracks(
            self.track_last_x,
            self.track_last_y,
            self.track_last_z,
            self.track_last_event,
            self.track_last_track,
            self.track_last_parent,
            self.track_last_volume,
            self.track_last_pdg,
            selected_event=selected_event,
            max_tracks_to_draw=max_tracks,
            only_primary=only_prim,
            only_secondary=only_sec,
            pdg_filter=pdg_filter if pdg_filter else None,
            title="Track Viewer",
        )

    def update_hit3d_from_ui(self):
        if self.hit3d_last_x is None or self.hit3d_last_y is None or self.hit3d_last_z is None:
            self.hit3d_canvas.plot_hits_3d(
                None, None, None,
                color=None,
                edep=None,
                title="3D Hit Haritası (veri yok)",
            )
            return

        x = np.asarray(self.hit3d_last_x, dtype=float)
        y = np.asarray(self.hit3d_last_y, dtype=float)
        z = np.asarray(self.hit3d_last_z, dtype=float)
        e = np.asarray(self.hit3d_last_edep, dtype=float) if self.hit3d_last_edep is not None else None
        pid = np.asarray(self.hit3d_last_pdg, dtype=int) if self.hit3d_last_pdg is not None else None
        vol = np.asarray(self.hit3d_last_volume, dtype=int) if self.hit3d_last_volume is not None else None
        ev  = np.asarray(self.hit3d_last_event, dtype=int) if self.hit3d_last_event is not None else None

        mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        if e is not None:
            mask = mask & np.isfinite(e)
        if pid is not None:
            mask = mask & np.isfinite(pid)
        if vol is not None:
            mask = mask & np.isfinite(vol)
        if ev is not None:
            mask = mask & np.isfinite(ev)

        pdg_filter = parse_int_list(self.hit3d_pdg_filter_edit.text())
        if pdg_filter and pid is not None:
            mask = mask & np.isin(pid, np.array(pdg_filter, dtype=int))

        vol_filter = parse_int_list(self.hit3d_volume_filter_edit.text())
        if vol_filter and vol is not None:
            mask = mask & np.isin(vol, np.array(vol_filter, dtype=int))

        if self.hit3d_use_event_cb.isChecked() and ev is not None:
            selected_event = int(self.track_event_spin.value())
            mask = mask & (ev == selected_event)

        x = x[mask]
        y = y[mask]
        z = z[mask]
        if e is not None:
            e = e[mask]
        if pid is not None:
            pid = pid[mask]
        if vol is not None:
            vol = vol[mask]

        n_points = x.size
        if n_points == 0:
            self.hit3d_canvas.plot_hits_3d(
                None, None, None,
                color=None,
                edep=None,
                title="3D Hit Haritası (filtre sonrası nokta yok)",
            )
            self.log("[BİLGİ] 3D hit haritası filtresi sonrası nokta kalmadı.")
            return

        max_points = int(self.hit3d_max_points_spin.value())
        if n_points > max_points:
            idx_sample = np.random.choice(n_points, size=max_points, replace=False)
            x = x[idx_sample]
            y = y[idx_sample]
            z = z[idx_sample]
            if e is not None:
                e = e[idx_sample]
            if pid is not None:
                pid = pid[idx_sample]
            if vol is not None:
                vol = vol[idx_sample]
            self.log(f"[BİLGİ] 3D hit haritası için {n_points} noktadan {max_points} tanesi örneklendi.")

        color_mode = self.hit3d_color_mode_combo.currentText()
        color_array = None
        cbar_label = ""

        if color_mode == "Edep" and e is not None:
            color_array = e
            cbar_label = "Edep"
        elif color_mode == "Z":
            color_array = z
            cbar_label = "Z"
        elif color_mode == "PDG" and pid is not None:
            color_array = pid.astype(float)
            cbar_label = "PDG ID"
        elif color_mode == "Volume" and vol is not None:
            color_array = vol.astype(float)
            cbar_label = "Volume ID"
        else:
            if e is not None:
                color_array = e
                cbar_label = "Edep"

        hit_title = self.hit3d_title_edit.text().strip() or "3D Hit Haritası"
        marker_size = int(self.hit3d_point_size_spin.value())

        self.hit3d_canvas.plot_hits_3d(
            x,
            y,
            z,
            color=color_array,
            edep=None,
            title=hit_title,
            xlabel=self.hit3d_last_xname,
            ylabel=self.hit3d_last_yname,
            zlabel=self.hit3d_last_zname,
            cbar_label=cbar_label,
            marker_size=marker_size,
        )


def main():
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
