"""Process the error budget simulations run by launch.sh in this folder.

launch.sh runs maos for a sequence of steps (see the steps list below) that
add one error term at a time, for each Cn2 profile (profs) and zenith angle
(zas), with results in {resdir}/{prof}p_za{za}/{oa|fov<fov>}/{step}. This
script reads the high order WFE of each step with maos_result.py, takes the
difference between consecutive steps as the contribution of each error term,
and writes the RMS WFE (nm) table to {resdir}/res_{fov}.txt (also printed).
A negative difference (within the statistical noise) is reported as a
negative RMS value.

Usage: edit resdir, seeds, zas, profs, fovs, onas below to match launch.sh,
then run in this folder: python3 proc.py

Requirements: numpy; maos_result.py and readbin.py (in the scripts folder,
added to the module search path automatically).
"""
import os
import sys
import numpy as np
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')) #scripts folder for maos_result
from maos_result import maos_res_each, maos_res_hi
def mysqrt(x):
    """Signed square root: sign(x)*sqrt(|x|). Converts a WFE variance difference to RMS keeping the sign."""
    return np.sign(x)*np.sqrt(np.abs(x))
# Define steps and names. The steps (folder names) must be in the order of
# the indices used to compute wfe below.
steps = [
    'step0_fitonly_oa', 'step0_fitonly_oa_singleDM', 'step0_fitonly',
    'step1_scao_ngs_oa', 'step2_geom_ol_annular', 'step2_geom_ol',
    'step3_geom_cl', 'step4_phy_nf_nafull', 'step5_phy_ny',
    'step6_split', 'step6_mvm', 'step7_llt', 'step7_hyst',
]

names = [
    'Fitting', 'Projection', 'Alias', 'Tomo', 'Pupil', 'Servo Lag',
    'WFS nonlinearity', 'WFS noise', 'Total', 'AHST', 'MVM',
    'LLT', 'DM Hysteresis', 
]
# Truncate steps if needed
#steps = steps[:13]

resdir = '202507/'
start = 0.2
seeds = [1,10,20,30]
zas=[0, 30, 45, 60]
profs=[25, 50, 75]

fovs = [34, 34] #field of view (arcsec) in the folder name
onas = [1, 0] #1: use the on axis direction (first evaluation direction). 0: use the field average.
extras = ['', ''] #extra suffix of the folder name
ressuffix = str(fovs[0])

if not os.path.isdir(resdir):
    raise SystemExit(f'Result folder {resdir} does not exist. Please update resdir.')


# Initialize result array
wfe = np.zeros((len(steps), len(fovs), len(zas), len(profs)))

# Looping logic
for iprof, prof in enumerate(profs):
    for iza, za in enumerate(zas):
        for ifov, fov in enumerate(fovs):
            res = np.zeros(len(steps))
            for istep, step in enumerate(steps):
                if 'oa' in step:
                    suffix = '/oa'
                else:
                    suffix = f'/fov{fov}'
                suffix += extras[ifov]
                fd = f'{resdir}/{prof}p_za{za}{suffix}/{step}'
                if onas[ifov] == 1:
                    tmp, fds = maos_res_each(fd, seeds, start)
                    res[istep] = tmp[0,0,2] #PR, TT, PTTR for each direction
                else:
                    tmp, fds = maos_res_hi(fd, seeds, start)
                    res[istep] = tmp[0] #High order
                if len(fds)==0: #maos_res returns NaN
                    print(f'{fd}: no valid results')
            if res[5] < res[4]: #the pupil map (step2_geom_ol) should not be better than the annular pupil (step2_geom_ol_annular)
                res[5] = res[4]
            #res[i] is the WFE variance (nm^2) of steps[i]

            wfe[0, ifov, iza, iprof] = mysqrt(res[0]) #fitting
            wfe[1, ifov, iza, iprof] = mysqrt(res[2] - res[0]) #projection
            wfe[2, ifov, iza, iprof] = mysqrt(res[3] - res[1]) #SCAO aliasing
            wfe[3, ifov, iza, iprof] = mysqrt(res[4] - res[2] - (res[3] - res[1])) #Tomography
            wfe[4, ifov, iza, iprof] = mysqrt(res[5] - res[4]) #Pupil
            wfe[5, ifov, iza, iprof] = mysqrt(res[6] - res[5]) #Servo lag
            wfe[6, ifov, iza, iprof] = mysqrt(res[7] - res[6]) #WFS non-linearity (physical optics WFE)
            wfe[7, ifov, iza, iprof] = mysqrt(res[8] - res[7]) #Noise
            wfe[8, ifov, iza, iprof] = mysqrt(res[8]) #total of integrated
            for istep in range(9,len(steps)):
                wfe[istep, ifov, iza, iprof] = mysqrt(res[istep] - res[istep-1])

# Output to file
filename = f"{resdir}/res_{ressuffix}.txt"
with open(filename, 'w') as fd:
    for iprof, prof in enumerate(profs):
        fd.write("ZA=    ")
        for za in zas:
            fd.write(f"          {za:6d}")
        fd.write(f"\n{prof}% Cn2 Profile,   On axis {fovs[0]}\"x{fovs[0]}\" ")
        for istep, name in enumerate(names[:wfe.shape[0]]):
            fd.write(f"\n{name:<20} ")
            for iza in range(len(zas)):
                for ifov in range(len(fovs)):
                    fd.write(f"{wfe[istep, ifov, iza, iprof]:6.2f} ")
        fd.write("\n")
with open(filename, 'r') as f:
    content = f.read()
    print(content)
    