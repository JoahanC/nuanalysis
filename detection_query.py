from astropy.table import QTable, Table, Column, vstack
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
from tqdm import tqdm


total_table = Table.read("complete.tbl", format='ipac')
print(len(total_table['RA']))
ra_vals = total_table['RA']
radius = 13 * u.arcsecond
for idx in tqdm(range(len(total_table['RA']))):
    pos = SkyCoord(f"{total_table['RA'][idx]} {total_table['DEC'][idx]}", unit=(u.hourangle, u.deg), frame='fk5')
    print(float(total_table['PROB'][idx]))
    #result_table = Vizier.query_region(pos, radius=radius)
    #print(result_table)