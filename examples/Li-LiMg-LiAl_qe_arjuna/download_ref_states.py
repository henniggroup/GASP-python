from pymatgen.ext.matproj import MPRester

ID=1210753
api_key='mui9dtio4ndhyHOSaVmyCHa5Vs6wXVPP'
a = MPRester(api_key)
results = a.materials.search(material_ids=[f"mp-{ID}"], fields=["structure","formula_pretty"])
results = results[0]
initial_structures = results.structure #.make_supercell([2, 2, 2]) only use this when the number of atoms is less than 2
formula = results.formula_pretty
initial_structures.to(fmt='POSCAR',filename=f'./ref_states/POSCAR.{formula}')
print(len(initial_structures))
print(formula)
print(initial_structures.lattice.abc)
