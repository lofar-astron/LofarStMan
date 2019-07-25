# This script checks if LofarStMan and casacore's derivedmscal get the
# same values for UVWs of J2000 sources and of the sun.
# It uses the dataset L33277_SAP000_SB000_uv.MS/FIELD which has
# a phase center in J2000 coordinates.
# It sets the phase center to SUN to test its UVW coordinates.
# At the end it sets for each time slot the SUN's J2000 direction as
# a J2000 direction in the PHASE_DIR column to test if LofarStMan and
# DerivedMSCal give the same result as for a proper SUN direction.

from pyrap.tables import *
import numpy as np

# Set the phase dir in case it is incorrect.
print("Reset phase dir to J2000 ...")
t = table ('L33277_SAP000_SB000_uv.MS/FIELD', readonly=False, ack=False)
t.putcolkeyword ('PHASE_DIR', 'MEASINFO.Ref', 'J2000')
t.putcol ('PHASE_DIR', t.getcol('REFERENCE_DIR'))
t.close()

# Calculate the J2000 direction of the SUN for each time slot.
print("Calculate J2000 direction of SUN ...")
dird = taql('calc meas.j2000("SUN", [select unique TIME from L33277_SAP000_SB000_uv.MS] s)')
dirs = dird['0'];    # result is dict with one entry containing array with dirs

# Check if J2000 UVWs are fine.
# Note that LofarStMan calculates UVW.
print("Compare J2000 UVWs from LofarStMan and DerivedMSCal ...")
t = taql('select from L33277_SAP000_SB000_uv.MS where not all(near(UVW, mscal.uvw(), 1e-5))')
if t.nrows() > 0:
    print("***",t.nrows(),"rows mismatch between UVW and mscal.uvw()")
t.close()

# Find all rows for each time slot.
print("Find rows per time slot ...")
t = table ('L33277_SAP000_SB000_uv.MS', ack=False)
rowlist = []
for iter in t.iter('TIME'):
    rowlist.append (iter.rownumbers())
    iter.close()
t.close()

# Calculate the UVWs for the SUN using derivedmscal.
print("Compare mscal.uvw('SUN') and mscal.uvw() ...")
t = taql('select mscal.uvw("SUN") as uvw from L33277_SAP000_SB000_uv.MS')
uvwSun1 = t.getcol('uvw')
t.close()

# Put SUN in the reference type and get those UVWs again.
t = table ('L33277_SAP000_SB000_uv.MS/FIELD', readonly=False, ack=False)
t.putcolkeyword ('PHASE_DIR', 'MEASINFO.Ref', 'SUN')
t.close()
t = taql('select mscal.uvw() as uvw from L33277_SAP000_SB000_uv.MS')
uvwSun2 = t.getcol('uvw')
t.close()

# Check if equal.
res = np.where (abs(uvwSun1 - uvwSun2) > 0.00001)
if len(res[0]) > 0:
    print("*** diff between mscal() and mscal('SUN')")

# Check if SUN UVWs are fine.
# Note that LofarStMan calculates UVW.
print("Compare SUN UVWs from LofarStMan and DerivedMSCal ...")
t = taql('select from L33277_SAP000_SB000_uv.MS where not all(near(UVW, mscal.uvw(), 1e-5))')
if t.nrows() > 0:
    print("***",t.nrows(),"rows mismatch between UVW and mscal.uvw()")
t.close()

# Now loop through all timeslots.
print("Compare for each time slot UVW of SUN and of J2000's SUN ...")
for i in range(min(100,len(rowlist))):
    # Put the SUN J2000 dir in the PHASE_DIR column.
    # So for this time slot only we pretend the SUN has a J2000 direction.
    t = table ('L33277_SAP000_SB000_uv.MS/FIELD', readonly=False, ack=False)
    t.putcolkeyword ('PHASE_DIR', 'MEASINFO.Ref', 'J2000')
    t.putcell ("PHASE_DIR", 0, dirs[i:i+1,])
    t.close()
    # Get the table subset for the time slot.
    t = table ('L33277_SAP000_SB000_uv.MS', ack=False)
    t1 = t.selectrows(rowlist[i])
    t2 = t1.query ('not all(near(UVW, mscal.uvw(), 1e-5))')
    if t2.nrows() != 0:
        print("***",t2.nrows(),"rows mismatch between UVW and mscal.uvw() for time",i) 
    t2.close()
    t2 = t1.select('mscal.uvw() as uvw')   # get UVW using derivedmscal
    uvwMsCal = t2.getcol ('uvw')
    t2.close()
    t1.close()
    t.close()
    # The UVWs should match with the SUN ones calculated earlier.
    res =  np.where (abs(uvwSun1[rowlist[i][0]:rowlist[i][-1]+1,] - uvwMsCal)>0.00001)
    if len(res[0]) > 0:
        print("*** mscal difference for time",i)
