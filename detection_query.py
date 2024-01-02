"""
This file generates basic statistics on all current detections listed under the ssdx format.
"""
from astropy.table import Table
import matplotlib.lines as mlines
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np


# Reading in data from ipac summary tables, this is for SSD1
basic_table_fivek_1 = Table.read("ssd1_basic_5000_final.tbl", format='ipac')
basic_table_onek_1 = Table.read("ssd1_basic_1000_final.tbl", format='ipac')
basic_table_halfk_1 = Table.read("ssd1_basic_500_final.tbl", format='ipac')
poisson_table_fivek_1 = Table.read("ssd1_poisson_5000_final.tbl", format='ipac')
poisson_table_onek_1 = Table.read("ssd1_poisson_1000_final.tbl", format='ipac')
poisson_table_halfk_1 = Table.read("ssd1_poisson_500_final.tbl", format='ipac')
# Reading in data from ipac summary tables, this is for SSD2
basic_table_fivek_2 = Table.read("ssd2_basic_5000_final.tbl", format='ipac')
basic_table_onek_2 = Table.read("ssd2_basic_1000_final.tbl", format='ipac')
basic_table_halfk_2 = Table.read("ssd2_basic_500_final.tbl", format='ipac')
poisson_table_fivek_2 = Table.read("ssd2_poisson_5000_final.tbl", format='ipac')
poisson_table_onek_2 = Table.read("ssd2_poisson_1000_final.tbl", format='ipac')
poisson_table_halfk_2 = Table.read("ssd2_poisson_500_final.tbl", format='ipac')

"""sep = 210
mask_basic_table_fivek_1 = basic_table_fivek_1['SEP'] > sep
mask_basic_table_onek_1 = basic_table_onek_1['SEP'] > sep
mask_basic_table_halfk_1 = basic_table_halfk_1['SEP'] > sep
mask_basic_table_fivek_2 = basic_table_fivek_2['SEP'] > sep
mask_basic_table_onek_2 = basic_table_onek_2['SEP'] > sep
mask_basic_table_halfk_2 = basic_table_halfk_2['SEP'] > sep
mask_poisson_table_fivek_1 = poisson_table_fivek_1['SEP'] > sep
mask_poisson_table_onek_1 = poisson_table_onek_1['SEP'] > sep
mask_poisson_table_halfk_1 = poisson_table_halfk_1['SEP'] > sep
mask_poisson_table_fivek_2 = poisson_table_fivek_2['SEP'] > sep
mask_poisson_table_onek_2 = poisson_table_onek_2['SEP'] > sep
mask_poisson_table_halfk_2 = poisson_table_halfk_2['SEP'] > sep

basic_table_fivek_1 = basic_table_fivek_1[mask_basic_table_fivek_1]
basic_table_onek_1 = basic_table_onek_1[mask_basic_table_onek_1]
basic_table_halfk_1 = basic_table_halfk_1[mask_basic_table_halfk_1]
basic_table_fivek_2 = basic_table_fivek_2[mask_basic_table_fivek_2]
basic_table_onek_2 = basic_table_onek_2[mask_basic_table_onek_2]
basic_table_halfk_2 = basic_table_halfk_2[mask_basic_table_halfk_2]
poisson_table_fivek_1 = poisson_table_fivek_1[mask_poisson_table_fivek_1]
poisson_table_onek_1 = poisson_table_onek_1[mask_poisson_table_onek_1]
poisson_table_halfk_1 = poisson_table_halfk_1[mask_poisson_table_halfk_1]
poisson_table_fivek_2 = poisson_table_fivek_2[mask_poisson_table_fivek_2]
poisson_table_onek_2 = poisson_table_onek_2[mask_poisson_table_onek_2]
poisson_table_halfk_2 = poisson_table_halfk_2[mask_poisson_table_halfk_2]"""

# Initialize lists for storing position information for the sources located in the summary tables for SSD1
basic_source_ra_fivek_1 = []
basic_source_dec_fivek_1 = []
basic_source_ra_onek_1 = []
basic_source_dec_onek_1 = []
basic_source_ra_halfk_1 = []
basic_source_dec_halfk_1 = []
poisson_source_ra_fivek_1 = []
poisson_source_dec_fivek_1 = []
poisson_source_ra_onek_1 = []
poisson_source_dec_onek_1 = []
poisson_source_ra_halfk_1 = []
poisson_source_dec_halfk_1 = []

# Initialize lists for storing position information for the sources located in the summary tables for SSD2
basic_source_ra_fivek_2 = []
basic_source_dec_fivek_2 = []
basic_source_ra_onek_2 = []
basic_source_dec_onek_2 = []
basic_source_ra_halfk_2 = []
basic_source_dec_halfk_2 = []
poisson_source_ra_fivek_2 = []
poisson_source_dec_fivek_2 = []
poisson_source_ra_onek_2 = []
poisson_source_dec_onek_2 = []
poisson_source_ra_halfk_2 = []
poisson_source_dec_halfk_2 = []

# Stores the number of targets listed among the summary tables, the [:-2] is to ensure that SEQID 
# sequences are properly counted and not overcounted for SSD1
basic_target_fivek_1 = []
for seqid in basic_table_fivek_1["SEQID"]:
    basic_target_fivek_1.append(seqid[:-2])
basic_target_onek_1 = []
for seqid in basic_table_onek_1["SEQID"]:
    basic_target_onek_1.append(seqid[:-2])
basic_target_halfk_1 = []
for seqid in basic_table_halfk_1["SEQID"]:
    basic_target_halfk_1.append(seqid[:-2])
poisson_target_fivek_1 = []
for seqid in poisson_table_fivek_1["SEQID"]:
    poisson_target_fivek_1.append(seqid[:-2])
poisson_target_onek_1 = []
for seqid in poisson_table_onek_1["SEQID"]:
    poisson_target_onek_1.append(seqid[:-2])
poisson_target_halfk_1 = []
for seqid in poisson_table_halfk_1["SEQID"]:
    poisson_target_halfk_1.append(seqid[:-2])

# Stores the number of targets listed among the summary tables, the [:-2] is to ensure that SEQID 
# sequences are properly counted and not overcounted for SSD2
basic_target_fivek_2 = []
for seqid in basic_table_fivek_2["SEQID"]:
    basic_target_fivek_2.append(seqid[:-2])
basic_target_onek_2 = []
for seqid in basic_table_onek_2["SEQID"]:
    basic_target_onek_2.append(seqid[:-2])
basic_target_halfk_2 = []
for seqid in basic_table_halfk_2["SEQID"]:
    basic_target_halfk_2.append(seqid[:-2])
poisson_target_fivek_2 = []
for seqid in poisson_table_fivek_2["SEQID"]:
    poisson_target_fivek_2.append(seqid[:-2])
poisson_target_onek_2 = []
for seqid in poisson_table_onek_2["SEQID"]:
    poisson_target_onek_2.append(seqid[:-2])
poisson_target_halfk_2 = []
for seqid in poisson_table_halfk_2["SEQID"]:
    poisson_target_halfk_2.append(seqid[:-2])

# Calculate the total number of SEQIDS and unique targets identified
basic_n_seqids_1 = len(set(basic_table_fivek_1["SEQID"]).union(set(basic_table_onek_1["SEQID"])).union(set(basic_table_halfk_1["SEQID"])))
basic_n_targets_1 = len(set(basic_target_fivek_1).union(set(basic_target_onek_1)).union(set(basic_target_halfk_1)))
basic_n_seqids_2 = len(set(basic_table_fivek_2["SEQID"]).union(set(basic_table_onek_2["SEQID"])).union(set(basic_table_halfk_2["SEQID"])))
basic_n_targets_2 = len(set(basic_target_fivek_2).union(set(basic_target_onek_2)).union(set(basic_target_halfk_2)))
basic_all_targets = set(basic_target_fivek_1).union(set(basic_target_onek_1)).union(set(basic_target_halfk_1)).union(set(basic_target_fivek_2)).union(set(basic_target_onek_2)).union(set(basic_target_halfk_2))
poisson_n_seqids_1 = len(set(poisson_table_fivek_1["SEQID"]).union(set(poisson_table_onek_1["SEQID"])).union(set(poisson_table_halfk_1["SEQID"])))
poisson_n_targets_1 = len(set(poisson_target_fivek_1).union(set(poisson_target_onek_1)).union(set(poisson_target_halfk_1)))
poisson_n_seqids_2 = len(set(poisson_table_fivek_2["SEQID"]).union(set(poisson_table_onek_2["SEQID"])).union(set(poisson_table_halfk_2["SEQID"])))
poisson_n_targets_2 = len(set(poisson_target_fivek_2).union(set(poisson_target_onek_2)).union(set(poisson_target_halfk_2)))
poisson_all_targets = set(poisson_target_fivek_1).union(set(poisson_target_onek_1)).union(set(poisson_target_halfk_1)).union(set(poisson_target_fivek_2)).union(set(poisson_target_onek_2)).union(set(poisson_target_halfk_2))

# Calculate the total number of targets per timescale class
basic_n_target_500 = len(set(basic_target_halfk_2).union(set(basic_target_halfk_1)))
basic_n_target_1000 = len(set(basic_target_onek_2).union(set(basic_target_onek_1)))
basic_n_target_5000 = len(set(basic_target_fivek_2).union(set(basic_target_fivek_1)))
poisson_n_target_500 = len(set(poisson_target_halfk_2).union(set(poisson_target_halfk_1)))
poisson_n_target_1000 = len(set(poisson_target_onek_2).union(set(poisson_target_onek_1)))
poisson_n_target_5000 = len(set(poisson_target_fivek_2).union(set(poisson_target_fivek_1)))

# Create list of source ras and sources decs from all detections
basic_added_set = []
basic_src_ras = []
basic_src_decs = []
for idx, seqid in enumerate(basic_table_fivek_1['SEQID']):
    if seqid[:-2] not in basic_added_set:
        basic_added_set.append(seqid[:-2])
        basic_src_ras.append(basic_table_fivek_1["RA"][idx])
        basic_src_decs.append(basic_table_fivek_1["DEC"][idx])
for idx, seqid in enumerate(basic_table_onek_1['SEQID']):
    if seqid[:-2] not in basic_added_set:
        basic_added_set.append(seqid[:-2])
        basic_src_ras.append(basic_table_onek_1["RA"][idx])
        basic_src_decs.append(basic_table_onek_1["DEC"][idx])
for idx, seqid in enumerate(basic_table_halfk_1['SEQID']):
    if seqid[:-2] not in basic_added_set:
        basic_added_set.append(seqid[:-2])
        basic_src_ras.append(basic_table_halfk_1["RA"][idx])
        basic_src_decs.append(basic_table_halfk_1["DEC"][idx])
for idx, seqid in enumerate(basic_table_fivek_2['SEQID']):
    if seqid[:-2] not in basic_added_set:
        basic_added_set.append(seqid[:-2])
        basic_src_ras.append(basic_table_fivek_2["RA"][idx])
        basic_src_decs.append(basic_table_fivek_2["DEC"][idx])
for idx, seqid in enumerate(basic_table_onek_2['SEQID']):
    if seqid[:-2] not in basic_added_set:
        basic_added_set.append(seqid[:-2])
        basic_src_ras.append(basic_table_onek_2["RA"][idx])
        basic_src_decs.append(basic_table_onek_2["DEC"][idx])
for idx, seqid in enumerate(basic_table_halfk_2['SEQID']):
    if seqid[:-2] not in basic_added_set:
        basic_added_set.append(seqid[:-2])
        basic_src_ras.append(basic_table_halfk_2["RA"][idx])
        basic_src_decs.append(basic_table_halfk_2["DEC"][idx])
poisson_added_set = []
poisson_src_ras = []
poisson_src_decs = []
for idx, seqid in enumerate(poisson_table_fivek_1['SEQID']):
    if seqid[:-2] not in poisson_added_set:
        poisson_added_set.append(seqid[:-2])
        poisson_src_ras.append(poisson_table_fivek_1["RA"][idx])
        poisson_src_decs.append(poisson_table_fivek_1["DEC"][idx])
for idx, seqid in enumerate(poisson_table_onek_1['SEQID']):
    if seqid[:-2] not in poisson_added_set:
        poisson_added_set.append(seqid[:-2])
        poisson_src_ras.append(poisson_table_onek_1["RA"][idx])
        poisson_src_decs.append(poisson_table_onek_1["DEC"][idx])
for idx, seqid in enumerate(poisson_table_halfk_1['SEQID']):
    if seqid[:-2] not in poisson_added_set:
        poisson_added_set.append(seqid[:-2])
        poisson_src_ras.append(poisson_table_halfk_1["RA"][idx])
        poisson_src_decs.append(poisson_table_halfk_1["DEC"][idx])
for idx, seqid in enumerate(poisson_table_fivek_2['SEQID']):
    if seqid[:-2] not in poisson_added_set:
        poisson_added_set.append(seqid[:-2])
        poisson_src_ras.append(poisson_table_fivek_2["RA"][idx])
        poisson_src_decs.append(poisson_table_fivek_2["DEC"][idx])
for idx, seqid in enumerate(poisson_table_onek_2['SEQID']):
    if seqid[:-2] not in poisson_added_set:
        poisson_added_set.append(seqid[:-2])
        poisson_src_ras.append(poisson_table_onek_2["RA"][idx])
        poisson_src_decs.append(poisson_table_onek_2["DEC"][idx])
for idx, seqid in enumerate(poisson_table_halfk_2['SEQID']):
    if seqid[:-2] not in poisson_added_set:
        poisson_added_set.append(seqid[:-2])
        poisson_src_ras.append(poisson_table_halfk_2["RA"][idx])
        poisson_src_decs.append(poisson_table_halfk_2["DEC"][idx])


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def sep_filter_detections(sep):
    mask_basic_table_fivek_1 = basic_table_fivek_1['SEP'] > sep
    mask_basic_table_onek_1 = basic_table_onek_1['SEP'] > sep
    mask_basic_table_halfk_1 = basic_table_halfk_1['SEP'] > sep
    mask_basic_table_fivek_2 = basic_table_fivek_2['SEP'] > sep
    mask_basic_table_onek_2 = basic_table_onek_2['SEP'] > sep
    mask_basic_table_halfk_2 = basic_table_halfk_2['SEP'] > sep

    sep_basic_table_fivek_1 = basic_table_fivek_1[mask_basic_table_fivek_1]
    sep_basic_table_onek_1 = basic_table_onek_1[mask_basic_table_onek_1]
    sep_basic_table_halfk_1 = basic_table_halfk_1[mask_basic_table_halfk_1]
    sep_basic_table_fivek_2 = basic_table_fivek_2[mask_basic_table_fivek_2]
    sep_basic_table_onek_2 = basic_table_onek_2[mask_basic_table_onek_2]
    sep_basic_table_halfk_2 = basic_table_halfk_2[mask_basic_table_halfk_2]
    
    
    # Create list of source ras and sources decs from all detections
    sep_basic_added_set = []
    sep_basic_src_ras = []
    sep_basic_src_decs = []
    for idx, seqid in enumerate(sep_basic_table_fivek_1['SEQID']):
        if seqid[:-2] not in sep_basic_added_set:
            sep_basic_added_set.append(seqid[:-2])
            sep_basic_src_ras.append(sep_basic_table_fivek_1["RA"][idx])
            sep_basic_src_decs.append(sep_basic_table_fivek_1["DEC"][idx])
    for idx, seqid in enumerate(sep_basic_table_onek_1['SEQID']):
        if seqid[:-2] not in sep_basic_added_set:
            sep_basic_added_set.append(seqid[:-2])
            sep_basic_src_ras.append(sep_basic_table_onek_1["RA"][idx])
            sep_basic_src_decs.append(sep_basic_table_onek_1["DEC"][idx])
    for idx, seqid in enumerate(sep_basic_table_halfk_1['SEQID']):
        if seqid[:-2] not in sep_basic_added_set:
            sep_basic_added_set.append(seqid[:-2])
            sep_basic_src_ras.append(sep_basic_table_halfk_1["RA"][idx])
            sep_basic_src_decs.append(sep_basic_table_halfk_1["DEC"][idx])
    for idx, seqid in enumerate(sep_basic_table_fivek_2['SEQID']):
        if seqid[:-2] not in sep_basic_added_set:
            sep_basic_added_set.append(seqid[:-2])
            sep_basic_src_ras.append(sep_basic_table_fivek_2["RA"][idx])
            sep_basic_src_decs.append(sep_basic_table_fivek_2["DEC"][idx])
    for idx, seqid in enumerate(sep_basic_table_onek_2['SEQID']):
        if seqid[:-2] not in sep_basic_added_set:
            sep_basic_added_set.append(seqid[:-2])
            sep_basic_src_ras.append(sep_basic_table_onek_2["RA"][idx])
            sep_basic_src_decs.append(sep_basic_table_onek_2["DEC"][idx])
    for idx, seqid in enumerate(sep_basic_table_halfk_2['SEQID']):
        if seqid[:-2] not in sep_basic_added_set:
            sep_basic_added_set.append(seqid[:-2])
            sep_basic_src_ras.append(sep_basic_table_halfk_2["RA"][idx])
            sep_basic_src_decs.append(sep_basic_table_halfk_2["DEC"][idx])
    
    im_range = np.arange(0, 800, 1)
    x_pix = np.array(sep_basic_table_halfk_1['XPIX'], dtype='f2')
    y_pix = np.array(sep_basic_table_halfk_1['YPIX'], dtype='f2')
    fig, ax = plt.subplots()
    pixelplot = ax.hist2d(x_pix, y_pix, bins = [im_range, im_range], cmap=discrete_cmap(5, "gnuplot2"))
    plt.xlim(250, 750)
    plt.ylim(250, 750)
    ax.set_xlabel("X Pixel")
    ax.set_ylabel("Y Pixel")
    ax.set_title("Net Detections per Pixel", loc="left")
    ax.set_title(f"Sep: {sep} arcsec", loc="right")
    cb = plt.colorbar(pixelplot[3], orientation='vertical')
    plt.savefig(f"filter_{sep}.pdf", dpi=1500)
    plt.close()

    return None


seps = [120, 150, 180, 210]
for sep in seps:
    sep_filter_detections(sep)

fig, ax = plt.subplots(figsize=(12, 7), subplot_kw=dict(projection="aitoff"))

plot_basic_ras = []
plot_basic_decs = []
for idx in range(len(basic_added_set)):
    sk = SkyCoord(f"{basic_src_ras[idx]} {basic_src_decs[idx]}", unit=(u.hourangle, u.deg), frame='fk5')
    plot_basic_ras.append(sk.galactic.l.wrap_at('180d').radian)
    plot_basic_decs.append(sk.galactic.b.radian)
ax.scatter(plot_basic_ras, plot_basic_decs, s=30, marker='o', color='red')
plot_poisson_ras = []
plot_poisson_decs = []
for idx in range(len(poisson_added_set)):
    sk = SkyCoord(f"{poisson_src_ras[idx]} {poisson_src_decs[idx]}", unit=(u.hourangle, u.deg), frame='fk5')
    plot_poisson_ras.append(sk.galactic.l.wrap_at('180d').radian)
    plot_poisson_decs.append(sk.galactic.b.radian)
ax.scatter(plot_poisson_ras, plot_poisson_decs, s=30, marker='o', color='blue')
ax.set_ylabel(r"$b$", color='white', fontsize=18)
ax.set_xlabel(r"$l$", color='white', fontsize=18)
#ax.grid()
#plt.show()
square_fivek = mlines.Line2D([], [], color='red', marker='o', linestyle='None',
                          markersize=10, label=f'Eliminated ximage detections')
diamond_onek = mlines.Line2D([], [], color='blue', marker='o', linestyle='None',
                          markersize=10, label=f'Remaining ximage detections')

plt.legend(bbox_to_anchor=(1.05, 1), handles=[square_fivek, diamond_onek], loc='upper right', prop={'size': 15}, borderaxespad=-2)
ax.grid()
plt.show()



im_range = np.arange(0, 800, 1)
x_pix = np.array(basic_table_halfk_1['XPIX'], dtype='f2')
y_pix = np.array(basic_table_halfk_1['YPIX'], dtype='f2')
fig, ax = plt.subplots()
pixelplot = ax.hist2d(x_pix, y_pix, bins = [im_range, im_range], cmap=discrete_cmap(6, "gnuplot2"))
plt.xlim(250, 750)
plt.ylim(250, 750)
ax.set_xlabel("X Pixel")
ax.set_ylabel("Y Pixel")
ax.set_title("Net Detections per Pixel")
cb = plt.colorbar(pixelplot[3], orientation='vertical')
plt.show()








"""
radius = 13 * u.arcsecond
for idx in tqdm(range(len(total_table['RA']))):
    pos = SkyCoord(f"{total_table['RA'][idx]} {total_table['DEC'][idx]}", unit=(u.hourangle, u.deg), frame='fk5')
    #print(float(total_table['PROB'][idx]))
    #result_table = Vizier.query_region(pos, radius=radius)
    #print(result_table)

prob_arr = np.array(total_table['PVAL'], dtype='f2')
fig, ax = plt.subplots()
MIN, MAX = 10e-7, 10e-3
ax.hist(prob_arr, bins = 10 ** np.linspace(np.log10(MIN), np.log10(MAX), 50))
ax.set_xscale('log')
ax.set_ylabel('Number')
ax.set_xlabel('Poisson Probability')
plt.show()


count_arr = np.array(total_table['COUNTS'], dtype='f2')
fig, ax = plt.subplots()
MIN, MAX = 10e-4, 1
ax.hist(count_arr, bins = 10 ** np.linspace(np.log10(MIN), np.log10(MAX), 50))
ax.set_xscale('log')
ax.set_ylabel('Number')
ax.set_xlabel('Count Rate')
plt.show()


fig, ax = plt.subplots()
yMIN, yMAX = 10e-8, 10e-4
xMIN, xMAX = 10e-4, 0.1
ax.hist2d(count_arr, prob_arr, bins = [10 ** np.linspace(np.log10(xMIN), np.log10(xMAX), 300), 10 ** np.linspace(np.log10(yMIN), np.log10(yMAX), 200)])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Count Rate')
ax.set_ylabel('Poisson Probability')
plt.show()

im_range = np.arange(0, 800, 1)
x_pix = np.array(total_table['XPIX'], dtype='f2')
y_pix = np.array(total_table['YPIX'], dtype='f2')
print(min(x_pix))
print(max(x_pix))
print(min(y_pix))
print(max(y_pix))
fig, ax = plt.subplots()
ax.hist2d(x_pix, y_pix, bins = [im_range, im_range])
plt.xlim(250, 750)
plt.ylim(250, 750)
plt.show()

skycoords = []
for i in range(len(total_table['XPIX'])):
    skycoords.append(SkyCoord(f"{total_table['RA'][i]} {total_table['DEC'][i]}", unit=(u.hourangle, u.deg)))
ras = [sk.ra.deg for sk in skycoords]
decs = [sk.dec.deg for sk in skycoords]





count_arra = np.array(total_table['ACOUNTS'], dtype='f2')
count_arrb = np.array(total_table['BCOUNTS'], dtype='f2')
fig, ax = plt.subplots()
MIN, MAX = 1, 10000
ax.hist(count_arra, bins = 10 ** np.linspace(np.log10(MIN), np.log10(MAX), 50))
ax.set_xscale('log')
ax.set_ylabel('Number')
ax.set_xlabel('Counts in FPMA')
plt.show()

fig, ax = plt.subplots()
MIN, MAX = 1, 10000
ax.hist(count_arrb, bins = 10 ** np.linspace(np.log10(MIN), np.log10(MAX), 50))
ax.set_xscale('log')
ax.set_ylabel('Number')
ax.set_xlabel('Counts in FPMB')
plt.show()

count_arr = count_arra + count_arrb
fig, ax = plt.subplots()
yMIN, yMAX = 10e-8, 10e-4
xMIN, xMAX = 1, 10000
ax.hist2d(count_arr, prob_arr, bins = [10 ** np.linspace(np.log10(xMIN), np.log10(xMAX), 10), 10 ** np.linspace(np.log10(yMIN), np.log10(yMAX), 10)])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Total Count Rate')
ax.set_ylabel('Poisson Probability')
plt.show()

fig, ax = plt.subplots()
yMIN, yMAX = 10e-8, 10e-4
xMIN, xMAX = 1, 10000
ax.hist2d(count_arra, prob_arr, bins = [10 ** np.linspace(np.log10(xMIN), np.log10(xMAX), 10), 10 ** np.linspace(np.log10(yMIN), np.log10(yMAX), 10)])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Total Counts in FPMA')
ax.set_ylabel('Poisson Probability')
plt.show()

fig, ax = plt.subplots()
yMIN, yMAX = 10e-8, 10e-4
xMIN, xMAX = 1, 10000
ax.hist2d(count_arrb, prob_arr, bins = [10 ** np.linspace(np.log10(xMIN), np.log10(xMAX), 10), 10 ** np.linspace(np.log10(yMIN), np.log10(yMAX), 10)])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Total Counts in FPMB')
ax.set_ylabel('Poisson Probability')
plt.show()
"""