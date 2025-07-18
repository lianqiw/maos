import numpy as np
from maos_result import maos_res_each, maos_res_hi
def mysqrt(x):
    return np.sign(x)*np.sqrt(np.abs(x))
# Define steps and names
steps = [
    'step0_fitonly_oa', 'step0_fitonly_oa_singleDM', 'step0_fitonly',
    'step1_scao_ngs_oa', 'step2_geom_ol_annular', 'step2_geom_ol',
    'step3_geom_cl', 'step4_phy_nf_nafull', 'step5_phy_ny',
    'step6_split', 'step6_mvm', 'step7_llt', 'step7_hyst',
]

names = [
    'Fitting', 'Projection', 'Alias', 'Tomo', 'Pupil', 'Servo Lag',
    'WFS nonlinearity', 'WFS noise', 'Total', 'AHST', 'MVM',
    'LLT', 'DM Hystersis', 
]
# Truncate steps if needed
#steps = steps[:13]

resdir = '202507/'
start = 0.2
seeds = [1,10,20,30]
zas=[0, 30, 45, 60]
profs=[25, 50, 75]

fovs = [34, 34]
onas = [1, 0]
extras = ['', '']
ressuffix = str(fovs[0])


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
            if res[5] < res[4]:
                res[5] = res[4]

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
    