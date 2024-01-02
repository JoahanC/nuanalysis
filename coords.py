import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

ra = 108 #266.41667
dec = 45 # -29.00778
ra2 = "07 14 22.645"
dec2 = "+45 38 27.766"
coord = SkyCoord(ra=ra, dec=dec, unit='deg', frame="fk5")
detect_position = SkyCoord(f"{ra2} {dec2}", unit=(u.hourangle, u.deg), frame='fk5')
gal_coord = coord.galactic
gal_coord2 = detect_position.galactic
print(gal_coord.l.deg, gal_coord.b.deg)
print(gal_coord2.l.deg, gal_coord2.b.deg)




fig, ax = plt.subplots(figsize=(12, 7), subplot_kw=dict(projection="aitoff"))
#ax.scatter(3, 5, s=30, marker='o', color='red')
ax.scatter(gal_coord.l.wrap_at('180d').radian, gal_coord.b.radian, s=30, marker='o', color='red')
plt.show()
