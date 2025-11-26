#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Geant4-benzeri test verisi üreten ROOT dosyası generatorü (v3, track destekli).

Önceki v2'ye ek olarak:
- Her dosya için eventID + trackID + parentID + stepIndex üretilir.
- Her event içinde birden fazla track vardır (primary + secondary gibi).
- Her track birden fazla adımdan oluşur (step), böylece "yol" çizilebilir.
- volumeID, posX/Y/Z, Edep, particleID da her adım için verilir.

TTree:
    Adı: "Hits"  (steps/track hits gibi düşünülebilir)

Branch'ler:
    - eventID   (int32)
    - trackID   (int32)
    - parentID  (int32)  : 0 => primary, >0 => secondary trackID
    - stepIndex (int32)  : 0..Nstep-1
    - particleID(int32)  : PDG (11, 22, 13, 211 vs.)
    - volumeID  (int32)  : 0/1/2 region
    - posX,Y,Z  (float64, mm)
    - Edep      (float64)

Bu veriyle:
- 3D/2D haritalar,
- Zayıflama,
- Volume ve parçacık özetleri,
- Event bazlı track görüntüleme (GUI tarafında)
yapılabilir.

ROOT gerektirmez, sadece 'uproot' + 'numpy' yeterli.
"""

import os
import argparse
import numpy as np
import uproot


def generate_tracks_for_thickness(
    thickness_mm: float,
    n_events: int,
    avg_tracks_per_event: float = 2.0,
    mu_per_mm: float = 0.1,
    beam_radius_mm: float = 10.0,
    z_sigma_mm: float = 2.0,
    rng: np.random.Generator = None,
):
    """
    Belirli bir kalınlık için event + track + step yapısında sentetik veri üretir.

    Parametreler
    -----------
    thickness_mm : float
        Malzeme kalınlığı [mm].
    n_events : int
        Event sayısı.
    avg_tracks_per_event : float
        Event başına ortalama track sayısı (primary + secondary).
    mu_per_mm : float
        Zayıflama katsayısı [1/mm].
    beam_radius_mm : float
        X/Y beam yarıçapı (gaussian sigma ~ beam_radius_mm/2).
    z_sigma_mm : float
        Z etrafında saçılma genişliği.
    rng : np.random.Generator
        NumPy RNG.

    Dönüş
    -----
    dict
        Yukarıda açıklanan branch'lerin hepsini içerir.
    """
    if rng is None:
        rng = np.random.default_rng()

    event_ids = []
    track_ids = []
    parent_ids = []
    step_indices = []
    pdgs = []
    volume_ids = []
    posX = []
    posY = []
    posZ = []
    edeps = []

    track_counter = 0

    # Region-parametreli zayıflama
    def region_mu(volume_id):
        if volume_id == 0:
            return 0.3 * mu_per_mm
        elif volume_id == 1:
            return 1.0 * mu_per_mm
        else:
            return 1.5 * mu_per_mm

    for ev in range(n_events):
        # event başına track sayısı (en az 1)
        n_tracks = max(1, int(rng.poisson(lam=avg_tracks_per_event)))
        primary_track_ids = []

        for t in range(n_tracks):
            track_id = track_counter
            track_counter += 1

            # Primary / secondary ayrımı
            if t == 0:
                parent_id = 0  # primary
                primary_track_ids.append(track_id)
            else:
                if primary_track_ids and rng.random() < 0.7:
                    parent_id = int(rng.choice(primary_track_ids))
                else:
                    parent_id = 0

            # PDG seçimi
            pdg_options = np.array([11, 22, 13, 211], dtype=np.int32)
            pdg_probs = np.array([0.5, 0.3, 0.15, 0.05])
            pdg = int(rng.choice(pdg_options, p=pdg_probs))

            # Step sayısı
            n_steps = rng.integers(5, 30)

            # Başlangıç X/Y
            sigma_xy = beam_radius_mm / 2.0
            x0 = rng.normal(loc=0.0, scale=sigma_xy)
            y0 = rng.normal(loc=0.0, scale=sigma_xy)

            # Z boyunca yol:
            if thickness_mm <= 1e-6:
                # neredeyse boşluk
                z_steps = rng.normal(loc=0.0, scale=z_sigma_mm, size=n_steps)
                z_steps = np.clip(z_steps, 0.0, 1e-3)
            else:
                z0 = 0.0
                z1 = thickness_mm
                z_steps = np.linspace(z0, z1, n_steps)
                z_steps += rng.normal(loc=0.0, scale=z_sigma_mm * 0.3, size=n_steps)
                z_steps = np.clip(z_steps, 0.0, thickness_mm)

            # X/Y boyunca saçılma
            x_steps = x0 + rng.normal(loc=0.0, scale=sigma_xy * 0.2, size=n_steps)
            y_steps = y0 + rng.normal(loc=0.0, scale=sigma_xy * 0.2, size=n_steps)

            # volumeID ve step bazlı µ
            volume_step = np.zeros(n_steps, dtype=np.int32)
            if thickness_mm > 1e-6:
                z_norm = z_steps / thickness_mm
                volume_step[(z_norm >= 0.3) & (z_norm < 0.7)] = 1
                volume_step[z_norm >= 0.7] = 2

            mu_step = np.array([region_mu(int(v)) for v in volume_step], dtype=float)

            # Edep: track boyunca exponential step enerji depoziti
            base_edep = rng.exponential(scale=0.02, size=n_steps)
            edep_steps = base_edep * np.exp(-mu_step * z_steps)
            # kalınlığa bağlı global zayıflama
            edep_steps *= np.exp(-0.5 * mu_per_mm * thickness_mm)

            # Biraz threshold
            mask_e = edep_steps > 1e-5
            z_steps = z_steps[mask_e]
            x_steps = x_steps[mask_e]
            y_steps = y_steps[mask_e]
            edep_steps = edep_steps[mask_e]
            volume_step = volume_step[mask_e]

            n_keep = edep_steps.size
            if n_keep == 0:
                continue

            # Append
            for s in range(n_keep):
                event_ids.append(ev)
                track_ids.append(track_id)
                parent_ids.append(parent_id)
                step_indices.append(s)
                pdgs.append(pdg)
                volume_ids.append(int(volume_step[s]))
                posX.append(float(x_steps[s]))
                posY.append(float(y_steps[s]))
                posZ.append(float(z_steps[s]))
                edeps.append(float(edep_steps[s]))

    # Array'lere çevir
    event_ids = np.asarray(event_ids, dtype=np.int32)
    track_ids = np.asarray(track_ids, dtype=np.int32)
    parent_ids = np.asarray(parent_ids, dtype=np.int32)
    step_indices = np.asarray(step_indices, dtype=np.int32)
    pdgs = np.asarray(pdgs, dtype=np.int32)
    volume_ids = np.asarray(volume_ids, dtype=np.int32)
    posX = np.asarray(posX, dtype=float)
    posY = np.asarray(posY, dtype=float)
    posZ = np.asarray(posZ, dtype=float)
    edeps = np.asarray(edeps, dtype=float)

    return {
        "eventID": event_ids,
        "trackID": track_ids,
        "parentID": parent_ids,
        "stepIndex": step_indices,
        "particleID": pdgs,
        "volumeID": volume_ids,
        "posX": posX,
        "posY": posY,
        "posZ": posZ,
        "Edep": edeps,
    }


def write_root_file(filename: str, hits_dict: dict):
    """Verilen sözlüğü ROOT dosyasına yazar (TTree adı: 'Hits')."""
    entries = len(hits_dict["Edep"])
    with uproot.recreate(filename) as f:
        f["Hits"] = hits_dict
    total_edep = float(np.sum(hits_dict["Edep"]))
    print(f"[OK] Yazıldı: {filename} (TTree: Hits, entries={entries}, ΣEdep={total_edep:.3e})")


def parse_thickness_list(text_list):
    """Komut satırı string listesini float kalınlık listesine çevir."""
    out = []
    for t in text_list:
        try:
            out.append(float(t))
        except ValueError:
            raise ValueError(f"Kalınlık sayıya çevrilemedi: '{t}'")
    return out


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Geant4-benzeri shielding test ROOT verisi üreteci (v3, track destekli)."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./root_test_data_v3",
        help="Çıktı klasörü (varsayılan: ./root_test_data_v3)",
    )
    parser.add_argument(
        "--thicknesses",
        nargs="+",
        default=["0", "10", "20", "30", "50"],
        help="Kalınlık listesi [mm] (ör: --thicknesses 0 10 20 30 50)",
    )
    parser.add_argument(
        "--mu",
        type=float,
        default=0.1,
        help="Zayıflama katsayısı µ [1/mm] (varsayılan: 0.1)",
    )
    parser.add_argument(
        "--events",
        type=int,
        default=200,
        help="Event sayısı (varsayılan: 200)",
    )
    parser.add_argument(
        "--avg-tracks",
        type=float,
        default=2.0,
        help="Event başına ortalama track sayısı (varsayılan: 2.0)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Rastgele sayı üreteci seed (varsayılan: 12345)",
    )

    args = parser.parse_args()

    try:
        thickness_list_mm = parse_thickness_list(args.thicknesses)
    except ValueError as e:
        print(f"[HATA] {e}")
        return

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    n_events = max(1, args.events)
    mu_per_mm = args.mu
    avg_tracks = max(0.1, args.avg_tracks)

    rng = np.random.default_rng(args.seed)

    print(f"Çıktı klasörü: {output_dir}")
    print(f"Kalınlıklar (mm): {thickness_list_mm}")
    print(f"Event sayısı: {n_events}, event başına ort. track sayısı: {avg_tracks}")
    print(f"µ = {mu_per_mm} 1/mm, seed = {args.seed}\n")

    for d in thickness_list_mm:
        tracks = generate_tracks_for_thickness(
            thickness_mm=d,
            n_events=n_events,
            avg_tracks_per_event=avg_tracks,
            mu_per_mm=mu_per_mm,
            beam_radius_mm=10.0,
            z_sigma_mm=2.0,
            rng=rng,
        )
        d_int = int(round(d))
        filename = os.path.join(output_dir, f"shield_{d_int:03d}mm.root")
        write_root_file(filename, tracks)

    print("\nTüm dosyalar üretildi.")
    print("TTree adı: Hits")
    print("Branch'ler: eventID, trackID, parentID, stepIndex, particleID, volumeID, posX, posY, posZ, Edep")
    print("Bu dosyalar track bazlı görselleştirme yapan GUI ile uyumludur.")


if __name__ == "__main__":
    main()
