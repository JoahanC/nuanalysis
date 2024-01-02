import numpy as np
import matplotlib.pyplot as plt
from nuanalysis import *
from astropy.table import Table
from astropy import table



dataset = Table.read("poisson_500_raw_ssd2.tbl", format="ipac")
unique_by_name = table.unique(dataset, keys="SEQID")
print(unique_by_name)