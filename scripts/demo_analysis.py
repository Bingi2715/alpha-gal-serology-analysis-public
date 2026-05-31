from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "demo_data" / "simulated_alpha_gal_serology_demo.csv"
TABLE_DIR = ROOT / "outputs" / "tables"

ALPHA_GAL_CODE = "30039"
MEAT_ASSAYS = {
    "Beef": "50510S",
    "Pork": "52310S",
    "Lamb/Mutton": "54310S",
}


def classify_ige(value: float) -> str:
    if value < 0.10:
        return "negative"
    if value < 0.35:
        return "borderline"
    if value < 0.70:
        return "low positive"
    if value < 3.50:
        return "moderate positive"
    return "high positive"


def make_same_day_pairs(data: pd.DataFrame, meat_name: str, meat_code: str) -> pd.DataFrame:
    key_cols = ["patient_id", "specimen_collection_date"]
    alpha = (
        data.loc[data["result_code"] == ALPHA_GAL_CODE, key_cols + ["ige_value_kU_L"]]
        .rename(columns={"ige_value_kU_L": "alpha_gal_ige_kU_L"})
        .groupby(key_cols, as_index=False)
        .mean(numeric_only=True)
    )
    meat = (
        data.loc[data["result_code"] == meat_code, key_cols + ["ige_value_kU_L"]]
        .rename(columns={"ige_value_kU_L": f"{meat_name.lower().replace('/', '_')}_ige_kU_L"})
        .groupby(key_cols, as_index=False)
        .mean(numeric_only=True)
    )
    paired = alpha.merge(meat, on=key_cols, how="inner")
    paired.insert(2, "paired_assay", meat_name)
    return paired


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(DATA_PATH, dtype={"patient_id": str, "accession_id": str, "result_code": str})
    data["specimen_collection_date"] = pd.to_datetime(data["specimen_collection_date"]).dt.date
    data["ige_category"] = data["ige_value_kU_L"].apply(classify_ige)
    data["is_positive"] = data["ige_value_kU_L"] >= 0.35

    assay_counts = (
        data.groupby(["analyte_name", "result_code"], as_index=False)
        .size()
        .rename(columns={"size": "record_count"})
        .sort_values(["analyte_name", "result_code"])
    )

    alpha = data.loc[data["result_code"] == ALPHA_GAL_CODE].copy()
    alpha_summary = pd.DataFrame(
        [
            {
                "assay": "Alpha-Gal",
                "records": len(alpha),
                "unique_patients": alpha["patient_id"].nunique(),
                "positive_records": int(alpha["is_positive"].sum()),
                "positivity_percent": round(100 * alpha["is_positive"].mean(), 2),
                "median_ige_kU_L": round(float(alpha["ige_value_kU_L"].median()), 3),
            }
        ]
    )

    category_counts = (
        data.groupby(["analyte_name", "ige_category"], as_index=False)
        .size()
        .rename(columns={"size": "record_count"})
    )

    correlation_rows = []
    for meat_name, meat_code in MEAT_ASSAYS.items():
        paired = make_same_day_pairs(data, meat_name, meat_code)
        paired_file = TABLE_DIR / f"demo_same_day_alpha_gal_vs_{meat_name.lower().replace('/', '_')}.csv"
        paired.to_csv(paired_file, index=False)

        meat_value_col = [col for col in paired.columns if col.endswith("_ige_kU_L") and col != "alpha_gal_ige_kU_L"][0]
        if len(paired) >= 3:
            rho, p_value = spearmanr(
                np.log10(paired["alpha_gal_ige_kU_L"]),
                np.log10(paired[meat_value_col]),
            )
        else:
            rho, p_value = np.nan, np.nan

        correlation_rows.append(
            {
                "paired_assay": meat_name,
                "same_day_pairs": len(paired),
                "spearman_rho_log10_ige": round(float(rho), 4) if pd.notna(rho) else np.nan,
                "p_value": float(p_value) if pd.notna(p_value) else np.nan,
            }
        )

    correlation_summary = pd.DataFrame(correlation_rows)

    assay_counts.to_csv(TABLE_DIR / "demo_record_counts_by_assay.csv", index=False)
    alpha_summary.to_csv(TABLE_DIR / "demo_alpha_gal_positivity_summary.csv", index=False)
    category_counts.to_csv(TABLE_DIR / "demo_ige_category_counts.csv", index=False)
    correlation_summary.to_csv(TABLE_DIR / "demo_same_day_spearman_correlations.csv", index=False)

    print("\nRecord counts by assay")
    print(assay_counts.to_string(index=False))
    print("\nAlpha-gal positivity summary")
    print(alpha_summary.to_string(index=False))
    print("\nIgE category counts")
    print(category_counts.to_string(index=False))
    print("\nSame-day log10 IgE Spearman correlations")
    print(correlation_summary.to_string(index=False))
    print(f"\nDemo tables saved to: {TABLE_DIR}")


if __name__ == "__main__":
    main()
