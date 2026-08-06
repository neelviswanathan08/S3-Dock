
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

# =====================================================================
# 1. CONFIGURATION: Define your 3 folders here
# =====================================================================
FOLDERS = {
    "Positive Control": "/home/aicd3-2/Documents/RESULTS_PAPER_S3-Dock/CRIPT_POS_final_summary",
    "Negative Control": "/home/aicd3-2/Documents/RESULTS_PAPER_S3-Dock/CRIPT_NEGATIVE_final_summary",
    "Our Candidate": "/home/aicd3-2/Documents/RESULTS_PAPER_S3-Dock/CRIPT_OURS_final_summary"
}

OUTPUT_FILENAME = "Figure_Master_Comparison.png"

# =====================================================================
# 2. FILE IDENTIFICATION & DATA EXTRACTION
# =====================================================================
extracted_data = {}

for label, folder in FOLDERS.items():
    print(f"Scanning {label} in {folder}...")
    
    # Smart Identification based on your rules
    try:
        network_img = glob.glob(os.path.join(folder, "*network*.png"))[0]
        hotspot_img = glob.glob(os.path.join(folder, "*hotspot*.png"))[0]
        stability_img = glob.glob(os.path.join(folder, "*stability*.png"))[0]
        summary_csv = glob.glob(os.path.join(folder, "*.csv"))[0] # The only CSV
    except IndexError as e:
        print(f"ERROR: Missing a required file in {folder}. Check your filenames!")
        continue

    # Extract MM-GBSA Mean and SD from the summary CSV
    df = pd.read_csv(summary_csv)
    
    # Assuming the CSV has columns named 'Mean' and 'SD' or similar. 
    # Adjust these column names if your CSV headers are different!
    try:
        mean_val = df['Mean'].values[0] if 'Mean' in df.columns else df.iloc[:, 1].mean()
        sd_val = df['SD'].values[0] if 'SD' in df.columns else df.iloc[:, 1].std()
    except:
        # Fallback if CSV format is raw energies
        mean_val = df.select_dtypes(include=[np.number]).mean()[0]
        sd_val = df.select_dtypes(include=[np.number]).std()[0]

    extracted_data[label] = {
        "network": network_img,
        "hotspot": hotspot_img,
        "stability": stability_img,
        "mean": mean_val,
        "sd": sd_val
    }

# =====================================================================
# 3. BUILD THE PUBLICATION FIGURE (4 Rows, 3 Columns)
# =====================================================================
fig = plt.figure(figsize=(18, 20), dpi=300)
gs = GridSpec(4, 3, figure=fig, height_ratios=[1, 1, 1, 1.5], hspace=0.3, wspace=0.1)

labels = list(extracted_data.keys())

# Helper function to plot images cleanly
def plot_img(ax, img_path, title=""):
    img = mpimg.imread(img_path)
    ax.imshow(img)
    ax.axis('off')
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold', pad=10)

# Row 1: Stability (RMSD)
for i, label in enumerate(labels):
    ax = fig.add_subplot(gs[0, i])
    plot_img(ax, extracted_data[label]["stability"], title=f"{label}\nStability Profile")

# Row 2: Hotspots
for i, label in enumerate(labels):
    ax = fig.add_subplot(gs[1, i])
    plot_img(ax, extracted_data[label]["hotspot"], title="Binding Hotspots")

# Row 3: Interaction Networks
for i, label in enumerate(labels):
    ax = fig.add_subplot(gs[2, i])
    plot_img(ax, extracted_data[label]["network"], title="Interaction Network")

# Row 4: MM-GBSA Bar Chart (Spans across all 3 columns)
ax_bar = fig.add_subplot(gs[3, :])
means = [extracted_data[l]["mean"] for l in labels]
sds = [extracted_data[l]["sd"] for l in labels]
colors = ['#2b5c8f', '#d95f02', '#7570b3'] # Blue, Orange, Purple

bars = ax_bar.bar(labels, means, yerr=sds, capsize=8, color=colors, alpha=0.9, edgecolor='black', linewidth=1.5, width=0.5)

# Annotate the bars with the exact values
for bar, mean, sd in zip(bars, means, sds):
    y_pos = mean - sd - (abs(mean)*0.1) if mean < 0 else mean + sd + (abs(mean)*0.1)
    ax_bar.text(bar.get_x() + bar.get_width()/2.0, y_pos, 
            f"{mean:.2f} ± {sd:.2f} kcal/mol", 
            ha='center', va='bottom' if mean >= 0 else 'top', 
            fontsize=12, fontweight='bold')

ax_bar.set_ylabel("Binding Free Energy (MM-GBSA, kcal/mol)", fontsize=14, fontweight='bold')
ax_bar.set_title("Thermodynamic Validation Comparison", fontsize=16, fontweight='bold', pad=15)
ax_bar.axhline(0, color='black', linewidth=1.5, linestyle='--')
ax_bar.tick_params(axis='both', labelsize=12)

# Save the final masterpiece
plt.tight_layout()
plt.savefig(OUTPUT_FILENAME, bbox_inches='tight')
print(f"\n✅ SUCCESS: Master figure saved as {OUTPUT_FILENAME}")