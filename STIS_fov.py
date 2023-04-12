""" STIS_fov.py

Created on Apr 12, 2023
@author: Shin Satoh

Description:


Version
1.0.0 (Apr 12, 2023)
"""

import spiceypy as spice
import numpy as np

spice.furnsh('kernel/cassMetaK.txt')
csvdata = np.recfromcsv('doc/STIS_Jupiter_2014_2022.csv', encoding='utf_8_sig', skip_header=1, delimiter=',',
                        names=['sci_data_set_name', 'sci_targname', 'sci_ra', 'sci_dec', 'sci_refnum', 'sci_start_time', 'sci_stop_time', 'sci_actual_duration', 'sci_instrume', 'sci_aper_1234', 'sci_spec_1234', 'sci_central_wavelength', 'sci_pep_id', 'sci_pi_last_name', 'sci_release_date', 'sci_preview_name'])

utc0 = csvdata.sci_start_time

et_hst = spice.str2et(utc0)
