massfile = open("masslist.txt","r")

masses = massfile.readlines()

massfile.close()

elementlist = []
masslist = []

for line in masses:
    if len(line.split()[2])>3:
        element = line.split()[3]
        massf = line.split()[4]
        massf2 = massf.split("(")[0]


    else:
        element = line.split()[2]
        massf = line.split()[3]
        massf2 = massf.split("(")[0]


    mass = massf2.split(']')[0].replace('[','')
    mass = float(mass)
    elementlist.append(element)
    masslist.append(mass)

massdict = dict(zip(elementlist,masslist))

import re

def mass_from_atom(in_atom, multiplier=1):
    "Calcualtes mass from atoms with coefficients (e.g. Cl2)"
    coef_re = re.compile('[\d]+')
    atom_re = re.compile('[A-Z][a-z]*')
    atom = atom_re.search(in_atom).group() if atom_re.search(in_atom) else in_atom
    coefficient = float(coef_re.search(in_atom).group()) if coef_re.search(in_atom) else 1.0
    mass_total = massdict[atom] * multiplier * coefficient
    return mass_total

def split_atoms(in_atom, multiplier=1):
    "Splits apart adjacent atoms (e.g. ClO)"
    elem_coef = re.compile('([A-Z][a-z]*[\d]*)')
    atom_split = elem_coef.split(in_atom)
    atom_split = [i for i in atom_split if i != '']
    atom_masses = 0.0
    if atom_split:
        for i in atom_split:
            atom_masses += mass_from_atom(i, multiplier)
        return atom_masses
    else:
        return mass_from_atom(in_atom, multiplier)
        
def split_polyatom(brack_polyatom):
    "Separates polyatom of form (CO3)2 into CO3 and mass multiplier 2"
    poly_coef = re.compile('(\)[\d]*)')
    poly_split = poly_coef.split(brack_polyatom)
    poly_split = [i for i in poly_split if i != '']
    outside_coefficient, polyatom = poly_split[1], poly_split[0].replace('(','')
    outside_coefficient = float(outside_coefficient.replace(')',''))
    return polyatom, outside_coefficient

def molarmass(formula):
    "Returns molar mass of formula string (e.g. Ca(NO3)2)"
    mass_total = 0.0
    poly = re.compile('(\([\w]+\)[\d]*)')
    elem = re.compile('([A-Z][a-z]*[\d]*)')
    poly_from_elem = poly.split(formula)
    poly_from_elem = [i for i in poly_from_elem if i != '']
    for part in poly_from_elem:
        if poly.search(part):
            polyatom_func = split_polyatom(part)
            mass_total += split_atoms(polyatom_func[0], polyatom_func[1])
        else:
            elem_split = elem.split(part)
            elem_split = [z for z in elem_split if z != '']
            for i in elem_split:
                mass_total += split_atoms(i)
            
    return mass_total
