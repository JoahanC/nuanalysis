from astropy.table import QTable, Table, Column, vstack
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
from tqdm import tqdm


total_table = Table.read("complete.tbl", format='ipac')
ra_vals = total_table['RA']
radius = 32 * u.arcsecond
for idx in tqdm(range(len(total_table['RA']))):
    pos = SkyCoord(f"{total_table['RA'][idx]} {total_table['DEC'][idx]}", unit=(u.hourangle, u.deg))
    result_table = Vizier.query_region(pos, radius=radius)
    print(result_table)