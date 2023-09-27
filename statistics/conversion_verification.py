import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# Define pathing and aggregate data directories
ssd_index = 1
path_stub = f"/Volumes/data_ssd_{ssd_index}/bifrost_data/"
data_dirs = glob.glob(f"/Volumes/data_ssd_{ssd_index}/bifrost_data/*", recursive=True)

# Populate dictionary with file names of relevant binning flags
data_info = {}
for directory in data_dirs:
    index = directory.replace(path_stub, '')
    if index == "straycat":
        continue
    data_info[index] = {}
    seqids = glob.glob(directory + "/*")
    for seqid in seqids:
        obsid = seqid.replace(directory + '/', '')
        data_info[index][obsid] = []
        bin_path_stub = seqid + "/event_cl/"
        binning_files = glob.glob(bin_path_stub + "5000_det1convert_flag.txt")
        for file in binning_files:
            data_info[index][obsid].append(file.replace(bin_path_stub, ''))
        binning_files = glob.glob(bin_path_stub + "1000_det1convert_flag.txt")
        for file in binning_files:
            data_info[index][obsid].append(file.replace(bin_path_stub, ''))
        binning_files = glob.glob(bin_path_stub + "500_det1convert_flag.txt")
        for file in binning_files:
            data_info[index][obsid].append(file.replace(bin_path_stub, ''))
        evt_files = glob.glob(bin_path_stub + "*01_cl.evt")
        if len(evt_files) == 2:
            data_info[index][obsid].append("both_present")

# Perform population statistics
dirs = []
vals_fivek = []
vals_onek = []
vals_halfk = []
net_vals = []
rel_vals = []
for directory in data_info:
    dirs.append(directory)
    count_fivek = 0
    count_onek = 0
    count_halfk = 0
    count_rel = 0
    for seqid in data_info[directory]:
        if "both_present" in data_info[directory][seqid]:
            count_rel += 1
        if "5000_det1convert_flag.txt" in data_info[directory][seqid]:
            count_fivek += 1
        if "1000_det1convert_flag.txt" in data_info[directory][seqid]:
            count_onek += 1
        if "500_det1convert_flag.txt" in data_info[directory][seqid]:
            count_halfk += 1

    vals_fivek.append(count_fivek)
    vals_onek.append(count_onek)
    vals_halfk.append(count_halfk)
    
    net_vals.append(len(data_info[directory]))
    rel_vals.append(count_rel)


# Creating bar plot
fig, ax = plt.subplots()
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(10, 8)

height = 1
y_pos = np.arange(1, 3*len(dirs), 3)

for idx, val in enumerate(rel_vals):
    ax.barh(y_pos+height, rel_vals, height=3*height, color="white", edgecolor='k', linewidth=2)

ax.barh(y_pos, vals_fivek, height=height, label="5000", edgecolor='k', color="#3498DB")
ax.barh(y_pos+height, vals_onek, height=height, label="1000", edgecolor='k', color="#8E44AD")
ax.barh(y_pos+2*height, vals_halfk, height=height, label="500", edgecolor='k', color="#E74C3C")

percent_vals_fivek = list(np.round(100 * np.array(vals_fivek) / np.array(rel_vals), 1))
percent_vals_onek = list(np.round(100 * np.array(vals_onek) / np.array(rel_vals), 1))
percent_vals_halfk = list(np.round(100 * np.array(vals_halfk) / np.array(rel_vals), 1))

for i, v in enumerate(percent_vals_fivek):
    x = vals_fivek[i]
    ax.text(0, 3*i + 0.55, "5000 seconds", color="black")
    if v == 100.:
        ax.text(x-4, 3*i + 0.55, "100%", color="black")
    else:
        ax.text(x-4.5, 3*i + 0.55, str(v) + '%', color="black")
for i, v in enumerate(percent_vals_onek):
    x = vals_onek[i]
    ax.text(0, 3*i + 1.55, "1000 seconds", color="black")
    if v == 100.:
        ax.text(x-4, 3*i + 1.55, "100%", color="black")
    else:
        ax.text(x-4.5, 3*i + 1.55, str(v) + '%', color="black")
    
for i, v in enumerate(percent_vals_halfk):
    x = vals_halfk[i]
    ax.text(0, 3*i + 2.55, "500 seconds", color="black")
    if v == 100.:
        ax.text(x-4, 3*i + 2.55, "100%", color="black")
    else:
        ax.text(x-4.5, 3*i + 2.55, str(v) + '%', color="black")
    

plt.yticks(y_pos+height, dirs, fontsize=15)
ax.set_title("Merging stage summary")
ax.set_ylabel("Directory")
ax.set_xlabel("Number of Observations")
plt.xlim(0, 75)
plt.ylim(0.5, 3*len(dirs) + 0.5)

plt.savefig(f"summary_conversion_{ssd_index}.pdf", dpi=1000)
plt.close()
