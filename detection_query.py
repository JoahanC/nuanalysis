from astropy.table import QTable, Table, Column, vstack
#from astroquery.vizier import Vizier
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib
from astropy.coordinates import SkyCoord
from astropy import units as u
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt


table_fivek = Table.read("complete_5000.tbl", format='ipac')
table_onek = Table.read("complete_1000.tbl", format='ipac')
table_halfk = Table.read("complete_500.tbl", format='ipac')
ra_fivek = []
dec_fivek = []
source_ra_fivek = []
source_dec_fivek = []
pvals_fivek = []

ra_onek = []
dec_onek = []
source_ra_onek = []
source_dec_onek = []
pvals_onek = []

ra_halfk = []
dec_halfk = []
source_ra_halfk = []
source_dec_halfk = []
pvals_halfk = []

fig, ax = plt.subplots(figsize=(12, 7), subplot_kw=dict(projection="aitoff"))
fig.patch.set_facecolor('black')
#ax.set_facecolor('gray')
for idx in range(len(table_fivek['RA'])):
    sk = SkyCoord(f"{table_fivek['RA'][idx]} {table_fivek['DEC'][idx]}", unit=(u.hourangle, u.deg))
    ra_fivek.append(sk.galactic.l.deg)
    dec_fivek.append(sk.galactic.b.deg)
    mpos = SkyCoord(f"{table_fivek['SRCRA'][idx]} {table_fivek['SRCDEC'][idx]}", unit=(u.deg, u.deg))
    source_ra_fivek.append(mpos.galactic.l.deg)
    source_dec_fivek.append(mpos.galactic.b.deg)
    pvals_fivek.append(table_fivek['PVAL'][idx])
#ax.scatter(source_ra_fivek, source_dec_fivek, s=20, marker='o', color='green')
cords1 = ax.scatter(ra_fivek, dec_fivek, s=40, marker='s', c=pvals_fivek, cmap='cool', norm=matplotlib.colors.LogNorm(), label="5000 sec")

for idx in range(len(table_onek['RA'])):
    sk = SkyCoord(f"{table_onek['RA'][idx]} {table_onek['DEC'][idx]}", unit=(u.hourangle, u.deg))
    ra_onek.append(sk.galactic.l.deg)
    dec_onek.append(sk.galactic.b.deg)
    mpos = SkyCoord(f"{table_onek['SRCRA'][idx]} {table_onek['SRCDEC'][idx]}", unit=(u.deg, u.deg))
    source_ra_fivek.append(mpos.galactic.l.deg)
    source_dec_fivek.append(mpos.galactic.b.deg)
    pvals_onek.append(table_onek['PVAL'][idx])
#ax.scatter(source_ra_onek, source_dec_onek, s=20, marker='o', color='green')
cords2 = ax.scatter(ra_onek, dec_onek, s=40, marker='d', c=pvals_onek, cmap='cool', norm=matplotlib.colors.LogNorm(), label="1000 sec")

for idx in range(len(table_halfk['RA'])):
    sk = SkyCoord(f"{table_halfk['RA'][idx]} {table_halfk['DEC'][idx]}", unit=(u.hourangle, u.deg))
    ra_halfk.append(sk.galactic.l.deg)
    dec_halfk.append(sk.galactic.b.deg)
    mpos = SkyCoord(f"{table_halfk['SRCRA'][idx]} {table_halfk['SRCDEC'][idx]}", unit=(u.deg, u.deg))
    source_ra_fivek.append(mpos.galactic.l.deg)
    source_dec_fivek.append(mpos.galactic.b.deg)
    pvals_halfk.append(table_halfk['PVAL'][idx])
#ax.scatter(source_ra_halfk, source_dec_halfk, s=20, marker='o', color='green', label="Main Sources")
cords3 = ax.scatter(ra_halfk, dec_halfk, s=40, marker='^', c=pvals_halfk, cmap='cool', norm=matplotlib.colors.LogNorm(), label="500 sec")

cords3.axes.tick_params(color='white', labelcolor='white')

ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=15)
#lg = plt.legend()
cb = plt.colorbar(cords1, orientation="horizontal", pad=0.08, shrink=0.7)
cb.set_label('Probability', color='white', fontsize=15)
cb.ax.xaxis.set_tick_params(color='white')
cb.outline.set_edgecolor('white')

#leg = ax.get_legend()
#leg.legend_handles[0].set_color('black')
#leg.legend_handles[1].set_color('black')

square_fivek = mlines.Line2D([], [], color='black', marker='s', linestyle='None',
                          markersize=10, label=f'5000 seconds; n = {len(table_fivek["RA"])}')
diamond_onek = mlines.Line2D([], [], color='black', marker='d', linestyle='None',
                          markersize=10, label=f'1000 seconds; n = {len(table_onek["RA"])}')
triangle_halfk = mlines.Line2D([], [], color='black', marker='^', linestyle='None',
                          markersize=10, label=f'500 seconds; n = {len(table_halfk["RA"])}')

plt.legend(bbox_to_anchor=(1.05, 1), handles=[square_fivek, diamond_onek, triangle_halfk], loc='upper right', prop={'size': 15}, borderaxespad=-2)

plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='white', fontsize=15)
plt.setp(plt.getp(cords3.axes, 'xticklabels'), color='black')
#ax.set_title('Current Candidates', color='white', loc='right')
#ax.set_ylabel(r"$b$", color='white', fontsize=18)
#ax.set_xlabel(r"$l$", color='white', fontsize=18)
ax.grid()
plt.show()
plt.savefig("candidates.pdf")


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