from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_PATH = ROOT / "demo_data" / "simulated_alpha_gal_serology_demo.csv"
RANDOM_SEED = 20260529
N_PATIENTS = 1400


def reported_result(value: float) -> str:
    return "<0.10" if value < 0.10 else f"{value:.2f}"


def random_collection_date(rng: np.random.Generator) -> pd.Timestamp:
    month_weights = np.array(
        [0.070, 0.065, 0.070, 0.075, 0.080, 0.085, 0.105, 0.110, 0.110, 0.100, 0.070, 0.060]
    )
    month_weights = month_weights / month_weights.sum()
    month = int(rng.choice(np.arange(1, 13), p=month_weights))
    year = int(rng.integers(2018, 2026))
    first = pd.Timestamp(year=year, month=month, day=1)
    last = first + pd.offsets.MonthEnd(0)
    return first + pd.Timedelta(days=int(rng.integers(0, last.day)))


def simulate_alpha_gal_ige(rng: np.random.Generator, age: int, month: int) -> float:
    age_effect = np.clip((age - 35) / 60, 0, 1) * 0.12
    season_effect = 0.10 if month in {7, 8, 9, 10} else 0.0
    p_positive = np.clip(0.18 + age_effect + season_effect, 0.08, 0.48)
    band = rng.choice(
        ["negative", "borderline", "low", "moderate", "high"],
        p=[1 - p_positive - 0.10, 0.10, p_positive * 0.35, p_positive * 0.45, p_positive * 0.20],
    )
    if band == "negative":
        return 0.05
    if band == "borderline":
        return float(rng.uniform(0.10, 0.34))
    if band == "low":
        return float(rng.uniform(0.35, 0.69))
    if band == "moderate":
        return float(10 ** rng.uniform(np.log10(0.70), np.log10(3.49)))
    return float(10 ** rng.uniform(np.log10(3.50), np.log10(80.0)))


def simulate_meat_ige(rng: np.random.Generator, alpha_ige: float, scale: float) -> float:
    if alpha_ige < 0.10 and rng.random() < 0.72:
        return 0.05
    value = alpha_ige * scale * float(rng.lognormal(mean=-0.20, sigma=0.75))
    if rng.random() < 0.18:
        value = float(rng.uniform(0.05, 0.34))
    return float(np.clip(value, 0.05, 65.0))


def main() -> None:
    rng = np.random.default_rng(RANDOM_SEED)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    patients = pd.DataFrame(
        {
            "patient_id": [f"DEMO-P{idx:05d}" for idx in range(1, N_PATIENTS + 1)],
            "year_of_birth": rng.integers(1935, 2021, size=N_PATIENTS),
        }
    )

    rows = []
    accession_counter = 1
    for patient in patients.itertuples(index=False):
        repeat_count = int(rng.choice([1, 2, 3, 4], p=[0.70, 0.20, 0.075, 0.025]))
        dates = sorted({random_collection_date(rng) for _ in range(repeat_count + 1)})[:repeat_count]
        for visit_date in dates:
            age = int(visit_date.year - patient.year_of_birth)
            alpha_ige = simulate_alpha_gal_ige(rng, age, int(visit_date.month))
            tested = [("Alpha-Gal", "30039", alpha_ige)]

            if rng.random() < 0.55:
                tested.extend(
                    [
                        ("Beef", "50510S", simulate_meat_ige(rng, alpha_ige, 0.78)),
                        ("Pork", "52310S", simulate_meat_ige(rng, alpha_ige, 0.62)),
                        ("Lamb/Mutton", "54310S", simulate_meat_ige(rng, alpha_ige, 0.55)),
                    ]
                )
            else:
                for name, code, scale in [
                    ("Beef", "50510S", 0.78),
                    ("Pork", "52310S", 0.62),
                    ("Lamb/Mutton", "54310S", 0.55),
                ]:
                    if rng.random() < 0.18:
                        tested.append((name, code, simulate_meat_ige(rng, alpha_ige, scale)))

            for name, code, value in tested:
                rows.append(
                    {
                        "patient_id": patient.patient_id,
                        "accession_id": f"DEMO-A{accession_counter:07d}",
                        "year_of_birth": int(patient.year_of_birth),
                        "specimen_collection_date": visit_date.strftime("%Y-%m-%d"),
                        "specimen_collection_time": (
                            f"{int(rng.integers(7, 18)):02d}:"
                            f"{int(rng.choice([0, 10, 15, 20, 30, 40, 45, 50])):02d}"
                        ),
                        "analyte_name": name,
                        "result_code": code,
                        "reported_result": reported_result(value),
                        "ige_value_kU_L": round(value, 3),
                        "unit": "kU/L",
                    }
                )
                accession_counter += 1

    pd.DataFrame(rows).to_csv(OUTPUT_PATH, index=False)
    print(f"Wrote {len(rows)} synthetic records for {N_PATIENTS} fake patients to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
