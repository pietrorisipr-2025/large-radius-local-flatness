#!/usr/bin/env python3
"""V1055B observable-radius calculation (illustrative, no fit).

Computes standard flat-LCDM distance scales from declared parameters:
Hubble radius, age, comoving CMB distance, particle horizon/observable
radius, observable diameter, and event horizon. This is a transparent
numerical illustration, not a parameter inference.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

C_KM_S = 299792.458
MPC_KM = 3.0856775814913673e19
SEC_PER_GYR = 31557600.0 * 1.0e9
MPC_TO_GLY = 3.261563777e-3  # 1 Mpc = 3.26156 million ly = 0.00326156 Gly

PARAMS = {
    "H0_km_s_Mpc": 67.4,
    "Omega_m": 0.315,
    "Omega_r": 9.2e-5,
}
PARAMS["Omega_Lambda"] = 1.0 - PARAMS["Omega_m"] - PARAMS["Omega_r"]
Z_CMB = 1089.9


def E_of_a(a: float) -> float:
    Om = PARAMS["Omega_m"]
    Or = PARAMS["Omega_r"]
    Ol = PARAMS["Omega_Lambda"]
    return math.sqrt(Or / a**4 + Om / a**3 + Ol)


def simpson(f, a: float, b: float, n: int = 200000) -> float:
    if n % 2:
        n += 1
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        x = a + i * h
        s += (4 if i % 2 else 2) * f(x)
    return s * h / 3.0


def integrate_split(f, intervals):
    total = 0.0
    for a, b, n in intervals:
        total += simpson(f, a, b, n)
    return total


def main() -> None:
    H0 = PARAMS["H0_km_s_Mpc"]
    RH_Mpc = C_KM_S / H0
    RH_Gly = RH_Mpc * MPC_TO_GLY
    H0_s = H0 / MPC_KM
    hubble_time_Gyr = 1.0 / H0_s / SEC_PER_GYR

    age_integral = integrate_split(lambda a: 1.0 / (a * E_of_a(a)), [
        (1e-8, 1e-4, 10000),
        (1e-4, 1e-2, 20000),
        (1e-2, 1.0, 200000),
    ])
    age_Gyr = hubble_time_Gyr * age_integral

    particle_integral = integrate_split(lambda a: 1.0 / (a * a * E_of_a(a)), [
        (1e-8, 1e-4, 10000),
        (1e-4, 1e-2, 20000),
        (1e-2, 1.0, 200000),
    ])
    particle_Gly = RH_Gly * particle_integral

    a_cmb = 1.0 / (1.0 + Z_CMB)
    cmb_integral = integrate_split(lambda a: 1.0 / (a * a * E_of_a(a)), [
        (a_cmb, 1e-2, 20000),
        (1e-2, 1.0, 200000),
    ])
    cmb_Gly = RH_Gly * cmb_integral

    event_integral = simpson(
        lambda u: 1.0 / math.sqrt(PARAMS["Omega_r"] * u**4 + PARAMS["Omega_m"] * u**3 + PARAMS["Omega_Lambda"]),
        0.0, 1.0, 200000
    )
    event_Gly = RH_Gly * event_integral

    values = [
        ("H0", H0, "km/s/Mpc"),
        ("Omega_m", PARAMS["Omega_m"], "dimensionless"),
        ("Omega_r", PARAMS["Omega_r"], "dimensionless"),
        ("Omega_Lambda", PARAMS["Omega_Lambda"], "dimensionless"),
        ("Hubble radius c/H0", RH_Gly, "Gly"),
        ("Integrated age", age_Gyr, "Gyr"),
        ("CMB comoving distance z=1089.9", cmb_Gly, "Gly"),
        ("Observable radius / particle horizon", particle_Gly, "Gly"),
        ("Observable diameter", 2.0 * particle_Gly, "Gly"),
        ("Event horizon illustrative", event_Gly, "Gly"),
        ("Horizon-chain one-dimensional union", 4.0 * particle_Gly, "Gly"),
    ]
    out = Path(__file__).resolve().parents[1] / "tables" / "v1055b_observable_radius_values.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value", "unit", "note"])
        for q, v, u in values:
            w.writerow([q, f"{v:.12g}", u, "illustrative flat LCDM; no fit"])
    print(out)
    for q, v, u in values:
        print(f"{q}: {v:.12g} {u}")

if __name__ == "__main__":
    main()
