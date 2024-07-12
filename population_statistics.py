"""
Helper script written to help understand detection population files to identify intersting targets.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u


dtimes = [500, 1000, 5000]
file_stub = ["ssd1_jan9_", "ssd2_jan9_"]
file_terminator = "_nss.tbl"
seqids = {}
plane_ids = set()
x_pixels = []
y_pixels = []
det1x = []
det1y = []
for stub in file_stub:
    for time in dtimes:
        filepath = stub + f"{time}" + file_terminator
        datatable = Table.read(filepath, format="ipac")
        
        for datum in datatable:
            if np.isnan(datum["DET1X"]) or np.isnan(datum["DET1Y"]):
                continue
            if datum["DET1X"] < 0 or datum["DET1Y"] < 0:
                continue
            x_pixels.append(datum["XPIX"])
            y_pixels.append(datum["YPIX"])
            det1x.append(datum["DET1X"])
            det1y.append(datum["DET1Y"])
            if datum["SEQID"] in seqids:
                seqids[datum["SEQID"]].append([datum["RA"], datum["DEC"]])
            if datum["SEQID"] in plane_ids:
                continue
            sk = SkyCoord(f"{datum['RA']} {datum['DEC']}", unit=(u.hourangle, u.deg), frame='fk5')
            if np.abs(sk.galactic.b.deg) < 15:
                plane_ids.add(datum["SEQID"])
                continue
            if datum["SEQID"] not in seqids:
                seqids[datum["SEQID"]] = [datum["RA"], datum["DEC"]]
            else:
                seqids[datum["SEQID"]].append([datum["RA"], datum["DEC"]])

# Histogram of detections in individual detections
counts = []
high_seqids = []
low_seqids = []
for seqid in seqids:
    counts.append(len(seqids[seqid]))
    if len(seqids[seqid]) > 100:
        high_seqids.append(seqid)
    else:
        low_seqids.append(seqid)

high_x = []
high_y = []
high_det1x = []
high_det1y = []
low_x = []
low_y = []
low_det1x = []
low_det1y = []
for datum in datatable:
    if datum["SEQID"] in high_seqids:
        high_x.append(datum["XPIX"])
        high_y.append(datum["YPIX"])
        high_det1x.append(datum["DET1X"])
        high_det1y.append(datum["DET1Y"])
    if datum["SEQID"] in low_seqids:
        low_x.append(datum["XPIX"])
        low_y.append(datum["YPIX"])
        low_det1x.append(datum["DET1X"])
        low_det1y.append(datum["DET1Y"])
    

fig, ax = plt.subplots()
ax.hist(counts, color="black", bins=np.logspace(np.log10(1),np.log10(10000), 50))
ax.set_xlabel("Number of detections")
ax.set_xscale("log")
ax.set_ylabel("Count")
ax.set_title(f"n = {len(seqids)}", loc="right")
#plt.show()
plt.close()

# Aitoff projection showing count histogram spatially

fig, ax = plt.subplots(figsize=(12, 7), subplot_kw=dict(projection="aitoff"))
for seqid in seqids:
    sk = SkyCoord(f"{seqids[seqid][0]} {seqids[seqid][1]}", unit=(u.hourangle, u.deg), frame='fk5')

    scale = len(seqids[seqid]) / np.max(counts) * 100
    ax.scatter(sk.galactic.l.wrap_at('180d').radian, sk.galactic.b.radian, s=scale, marker='o', color='red')
ax.grid()
plt.savefig("nss_aitoff.pdf", dpi=1000)
plt.close()

fig, ax = plt.subplots()
heatmap, _, _ = np.histogram2d(x_pixels, y_pixels, bins=800)
ax.imshow(heatmap)
ax.set_xlabel("Raw X Coordinate")
ax.set_ylabel("Raw Y Coordinate")
ax.set_title("Raw Locations (All)", loc="right")
plt.savefig("nss_raw_coords.pdf", dpi=1000)
plt.close()

fig, ax = plt.subplots()
heatmap, _, _ = np.histogram2d(high_x, high_y, bins=800)
ax.imshow(heatmap)
ax.set_xlabel("Raw X Coordinate")
ax.set_ylabel("Raw Y Coordinate")
ax.set_title("Raw Locations (Bright)", loc="right")
plt.savefig("nss_raw_coords_bright.pdf", dpi=1000)
plt.close()


fig, ax = plt.subplots()
heatmap, _, _ = np.histogram2d(low_x, low_y, bins=100)
ax.imshow(heatmap)
ax.set_xlabel("Raw X Coordinate")
ax.set_ylabel("Raw Y Coordinate")
ax.set_title("Raw Locations (Faint)", loc="right")
plt.savefig("nss_raw_coords_faint.pdf", dpi=1000)
plt.close()

fig, ax = plt.subplots()
heatmap, _, _ = np.histogram2d(det1x, det1y, bins=360)
ax.imshow(heatmap)
ax.set_xlabel("DET1X Coordinate")
ax.set_ylabel("DET1Y Coordinate")
ax.set_title("DET1 Locations (All)", loc="right")
plt.savefig("nss_det1_coords.pdf", dpi=1000)
plt.close()

fig, ax = plt.subplots()
heatmap, _, _ = np.histogram2d(high_det1x, high_det1y, bins=360)
ax.imshow(heatmap)
ax.set_xlabel("DET1X Coordinate")
ax.set_ylabel("DET1Y Coordinate")
ax.set_title("DET1 Locations (Bright)", loc="right")
plt.savefig("nss_det1_coords_bright.pdf", dpi=1000)
plt.close()

fig, ax = plt.subplots()
heatmap, _, _ = np.histogram2d(low_det1x, low_det1y, bins=360)
ax.imshow(heatmap)
ax.set_xlabel("DET1X Coordinate")
ax.set_ylabel("DET1Y Coordinate")
ax.set_title("DET1 Locations (Faint)", loc="right")
plt.savefig("nss_det1_coords_faint.pdf", dpi=1000)
plt.close()