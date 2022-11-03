from astroquery.simbad import Simbad
import numpy as np
import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astropy.coordinates import Angle
from astroplan import FixedTarget
from astroplan import Observer
from astropy.time import Time, TimeDelta

import re
from mechanize import Browser
url = "http://simbad.cds.unistra.fr/simbad/sim-fbasic"

def text_to_sky(ra_t, dec_t):
    ra = ra_t.split(":")
    dec = dec_t.split(":")
    ra = ra[0] + 'h' + ra[1] + 'm' + ra[2] +'s'
    dec = dec[0] + 'd' + dec[1] + 'm' + dec[2] +'s'
    gc_c = SkyCoord(ra, dec) # galactic center
    return gc_c


def query_stars(qrstr): 
    result = Simbad.query_criteria(qrstr)
    star_pos = SkyCoord(result["RA"], result["DEC"], unit=(u.hourangle, u.deg))

    total_vstack = []
    for index, row in df.iterrows():
        ts = get_tell_star(row, result, star_pos)

        ts['target'] = [row["Name"]]
        total_vstack = vstack([total_vstack, result])
            
    return total_vstack


def get_tell_star(row, result, star_pos, obstime = None, 
                  tobs = None, data_path = None):
    lasilla = Observer.at_site('La Silla Observatory')
    if tobs is None:
        tobs = row.iloc[0]["Tobs"] 
    
    ra = Angle(row["ra_agn_hms"], unit = u.hourangle)
    dec = Angle(row["dec_agn_deg"], unit=u.deg)
    ha = Angle(tobs, unit = u.hourangle)
    lst = (ra + 0.5* ha).to_string(sep = ":")
    tra = ra.to_string(sep = ":")
    tdec = dec.to_string(sep = ":")
    
    ideal_pos = SkyCoord(lst, tdec, unit=(u.hourangle, u.deg))
    sep = star_pos.separation(ideal_pos)

    coordinates = SkyCoord(tra, tdec, unit=(u.hourangle, u.deg))
    targ = FixedTarget(coord=coordinates, name=None)
    
    if obstime is None: 
        time = Time('2022-09-09 21:00:00')
        obstime = lasilla.target_meridian_transit_time(time = time, target = targ, which = "nearest")
        #obstime = time
    else:
        tm = TimeDelta(0.5*tobs * u.hour)
        obstime = obstime + tm

    airmass = lasilla.altaz(obstime, targ).secz.value
    print ("Target airmass at ",obstime.iso, ":", airmass)
    
    tell_star = FixedTarget(coord=star_pos, name=None)
    tm = TimeDelta(0.5 * tobs * u.hour)
    star_obstime = obstime + tm# + row["Tobs"]*u.hour
    #star_obstime = lasilla.target_meridian_transit_time(time = time+row["Tobs"]*u.hour, 
    #                                                    target = tell_star, which = "nearest")
    airmass_star = lasilla.altaz(star_obstime, tell_star).secz.value
    d_airmass = np.abs(airmass_star-airmass)
    #print (star_obstime, tm)
    print ("Telluric observed at ", star_obstime.iso)
    
    result["sep"] = sep
    result["airmass"] = airmass_star

    result["d_airmass"] = d_airmass
    
    filtered = result[(sep < 45*u.degree) & (d_airmass < 0.1)] # change from 15 to 45
    filtered.sort([ "d_airmass", "sep"])
    mstable = filtered

    if data_path is not None:
        filtered.write(data_path + row["Name"]+"telluric_star.csv", overwrite=True) 
    
    filtered = result[(sep < 45*u.degree) & (d_airmass < 0.05)] # change from 15 to 45

    if len(filtered) == 0:
        filtered = result[(sep < 45*u.degree) & (d_airmass < 0.1)] # change from 15 to 45

    ii = 1
    while len(filtered) == 0:
        filtered = result[(sep < (45+(15*ii))*u.degree) & (d_airmass < 0.1)] # change from 15 to 45
        ii += 1

    filtered.sort(["sep", "d_airmass"])
    filtered['target'] = [row["Name"]] * len(filtered)

        
    return filtered[0], mstable, result

def load_stars():

    qrstr = "dec<25  &  (sptype = 'A0V' | sptype = 'B0V' | sptype = 'B1V' |  sptype = 'B2V' | sptype = 'B3V') &  Kmag < 8"

    result = Simbad.query_criteria(qrstr)

    star_pos = SkyCoord(result["RA"], result["DEC"], unit=(u.hourangle, u.deg))
    
    return result, star_pos

def load_target(target_file, sheet_name=0):
    df = pd.read_excel(target_file, sheet_name=sheet_name) # removed skiprows
    return df


def simbad_search(id_str): 
    br = Browser()
    br.open(url)
    br.select_form(nr = 0)
    br["Ident"] = id_str
    response = br.submit()
    df_tmp = pd.read_html(response.get_data(), header = [1])
    str_simbad = str(df_tmp[0].loc[0])
    info_end = str_simbad.find("Origin of")
    spt_st = str_simbad.find("Spectral type:")
    basic_info = (str_simbad[0:info_end])
    spt_type = str_simbad[spt_st:spt_st+25]
    return basic_info, spt_type

def telluric_search(name, obstime=None, tobs=None, show_max=5,  backup=False):
    global result, star_pos, url
    
    if backup:
        tdf =  back_df
    else:
        tdf = df
    
    select = (tdf["Name"] == name)
    row = tdf[select]
    print (row)
    
    ts, tb, rr = get_tell_star(row, result, star_pos, 
                           obstime = obstime, tobs = tobs, data_path = None)
    
    for row in tb[0:show_max]:
        try:
            basic_info, spt_type = simbad_search(row["MAIN_ID"])
            print ( basic_info, spt_type, "Sep:", row["sep"], "d_am:" , row["d_airmass"],  row["RA"], row["DEC"], row["airmass"])
        except:
            print (row["MAIN_ID"], "Sep", row["sep"], "d_am" , row["d_airmass"], row["RA"], row["DEC"], 
                  row["airmass"])            
    return 

result, star_pos = load_stars()
#target_file = "tabs/target_final_sept2022_combined.xlsx"
target_file = "tabs/target_tech_night.xlsx"
#sheet_name = 1
df = load_target(target_file)#, sheet_name=sheet_name)
back_df = load_target("tabs/target_final_sept2022_add_backup_excel.xlsx")
