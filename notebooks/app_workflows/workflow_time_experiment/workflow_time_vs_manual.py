"""
Run this script to compare the runtime of the workflow to the manual time.

To run:

python3 workflow_time_vs_manual.py


"""

import pandas as pd

# Load into DataFrame
df = pd.read_csv("runtime_experiment.csv")

# If your CSV has column “runtime” rename it to runtime_s
if "runtime" in df.columns and "runtime_s" not in df.columns:
    df = df.rename(columns={"runtime": "runtime_s"})

# Extract number of plasmids from the filename suffix, e.g. “…_10.ipynb” → 10
df["plasmids"] = df["notebook"].str.extract(r"_(\d+)\.ipynb$", expand=False).astype(int)

# Calculate manual times and conversions
df["manual_time_s"] = df["plasmids"] * 720  # 720 s per plasmid
df["runtime_h"] = df["runtime_s"] / 3600
df["manual_time_h"] = df["manual_time_s"] / 3600
df["time_saved_s"] = df["manual_time_s"] - df["runtime_s"]
df["speed_up_factor"] = df["manual_time_s"] / df["runtime_s"]

# Sort by notebook name
df = df.sort_values("notebook")

# Round values for clarity
df["runtime_s"] = df["runtime_s"].round(2)
df["runtime_h"] = df["runtime_h"].round(4)
df["manual_time_s"] = df["manual_time_s"].astype(int)
df["manual_time_h"] = df["manual_time_h"].round(1)
df["time_saved_s"] = df["time_saved_s"].round(2)
df["speed_up_factor"] = df["speed_up_factor"].round(1)

# Reorder columns
df = df[
    [
        "notebook",
        "plasmids",
        "runtime_s",
        "runtime_h",
        "manual_time_s",
        "manual_time_h",
        "time_saved_s",
        "speed_up_factor",
    ]
]

# Write out the comparison CSV
df.to_csv("workflow_time_vs_manual.csv", index=False)
# Write to excel
df.to_excel("workflow_time_vs_manual.xlsx", index=False)

print("✔ Written workflow_time_vs_manual.csv")
