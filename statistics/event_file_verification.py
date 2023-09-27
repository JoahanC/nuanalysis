import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# Define pathing and aggregate data directories
ssd_index = 2
path_stub = f"/Volumes/data_ssd_{ssd_index}/"
data_dirs = glob.glob(f"/Volumes/data_ssd_{ssd_index}/*", recursive=True)

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

        binning_files = glob.glob(bin_path_stub + "*01_cl.evt")
        for file in binning_files:
            data_info[index][obsid].append(file.replace(bin_path_stub, ''))

# Perform population statistics
dirs = []
vals_fivek = []
vals_onek = []
vals_halfk = []
net_vals = []
for directory in data_info:
    dirs.append(directory)
    count_fivek = 0
    count_onek = 0
    count_halfk = 0
    for seqid in data_info[directory]:
        if "5000_flag.txt" in data_info[directory][seqid]:
            count_fivek += 1
        if "1000_flag.txt" in data_info[directory][seqid]:
            count_onek += 1
        if len(data_info[directory][seqid]) == 2:
            count_halfk += 1
    vals_fivek.append(count_fivek)
    vals_onek.append(count_onek)
    vals_halfk.append(count_halfk)
    net_vals.append(len(data_info[directory]))


# Creating bar plot
fig, ax = plt.subplots()
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(10, 8)

height = 1
y_pos = np.arange(1, 3*len(dirs), 3)

for idx, val in enumerate(net_vals):
    ax.barh(y_pos+height, net_vals, height=3*height, color="white", edgecolor='k', linewidth=2)


ax.barh(y_pos+height, vals_halfk, height=3*height, label="500", edgecolor='k', color="#E74C3C")


percent_vals_halfk = list(np.round(100 * np.array(vals_halfk) / np.array(net_vals), 1))
   
for i, v in enumerate(percent_vals_halfk):
    x = vals_halfk[i]
    if v == 100.:
        ax.text(x-6, 3*i + 1.55, "100%", color="black", fontsize=15)
    else:
        ax.text(x-6.5, 3*i + 1.55, str(v) + '%', color="black", fontsize=15)
    

plt.yticks(y_pos+height, dirs, fontsize=15)
ax.set_title("Merging stage summary")
ax.set_ylabel("Directory")
ax.set_xlabel("Number of Observations")
plt.xlim(0, 75)
plt.ylim(0.5, 3*len(dirs) + 0.5)

plt.savefig(f"summary_evtfile_{ssd_index}.pdf", dpi=1000)
