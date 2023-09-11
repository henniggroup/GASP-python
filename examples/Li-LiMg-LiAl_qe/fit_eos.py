from ase.io import read
from ase.eos import EquationOfState

from typing import List

from glob import glob

import numpy as np 


def FitEOS(output_file_lst, str = "sjeos",write_to_disk=True):
    es = []
    vs = []
    for file_name in output_file_lst:
        try:
            e = read(file_name).get_potential_energy()
        except: 
            continue
        v = read(file_name).get_volume()
        es.append(float(e))
        vs.append(float(v))

    f =  open("eos_fitting.txt", "w")    
    if len(vs) < 3:
        scale = 0
        if write_to_disk:
            f.write(str(scale))
            f.close()
        return scale 
    eos = EquationOfState(vs, es, eos_algorithm)
    v0, _, _ = eos.fit()
    minvol = min(vs)
    maxvol = max(vs)
    # fig, ax = plt.subplots()
    # if on_k8s:
    #     eos.plot(f"{eos_plot_name}_eos.png", ax=ax)
    # else:
    #     eos.plot(f"{eos_plot_name}_eos.png",ax=ax)
    # ax.set_xlabel(r"Volume ($\AA$)")
    # ax.set_ylabel(r"Energy (eV)")
    # plt.close(fig)
    if not (minvol < v0 < maxvol):
        scale = 0
        if write_to_disk:
            f.write(str(scale))
            f.close()
        return scale
    else:
        try:
            scale_1_atoms_file_name = glob("*_1.pwi.pwo")[0]
            # scale_1_atoms_file_name =  np.array(output_file_lst)[np.array(scale_lst)==1][0]
            atoms=read(scale_1_atoms_file_name)
        except:
            scale_1_atoms_file_name = glob("*_1.pwi")[0]
            atoms=read(scale_1_atoms_file_name)
        v = atoms.get_volume()
        scale = (v0 / v)**(1/3)
        if write_to_disk:
            f.write(str(scale))
            f.close()

        return float(scale)
    
if __name__ == "__main__":
    output_file_lst = sorted(glob("*.pwo"))
    scale = FitEOS(output_file_lst)