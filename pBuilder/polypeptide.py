import mbuild as mb
from .aminos import *
import sys
import numpy as np

aminos = {}
aminos[('alanine','ala','A')] = 'CC(C(=O)O)N'
aminos[('arginine','arg','R')] = 'C(C[C@@H](C(=O)O)N)CN=C(N)N'
aminos[('asparagine' ,'asn', 'N')] = 'C([C@@H](C(=O)O)N)C(=O)N'
aminos[('aspartic acid', 'asp','D')] = 'C([C@@H](C(=O)O)N)C(=O)O'
aminos[('cysteine', 'cys', 'C')] = 'C([C@@H](C(=O)O)N)S'
aminos[('glutamine', 'gln', 'Q')] = 'C(CC(=O)N)[C@@H](C(=O)O)N'
aminos[('glutamic acid', 'glu', 'E')] = 'C(CC(=O)O)C(C(=O)O)N'
aminos[('glycine', 'gly', 'G')] = 'C(C(=O)O)N'
aminos[('histidine', 'his', 'H')] = 'C1=C(NC=N1)C[C@@H](C(=O)O)N'
aminos[('isoleucine', 'ile', 'I')] = 'CC[C@H](C)[C@@H](C(=O)O)N'
aminos[('leucine', 'leu', 'L')] = 'CC(C)C[C@@H](C(=O)O)N'
aminos[('lysine', 'lys', 'K')] = 'C(CCN)C[C@@H](C(=O)O)N'
aminos[('methionine', 'met', 'M')] = 'CSCC[C@@H](C(=O)O)N'
aminos[('phenylalanine', 'phe', 'F')] = 'C1=CC=C(C=C1)C[C@@H](C(=O)O)N'
aminos[('proline', 'pro', 'P')] = 'C1C[C@H](NC1)C(=O)O'
aminos[('serine', 'ser', 'S')] = 'C([C@@H](C(=O)O)N)O'
aminos[('threonine', 'thr', 'T')] = 'C[C@H]([C@@H](C(=O)O)N)O'
aminos[('tryptophan', 'trp', 'W')] = 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N'
aminos[('tyrosine', 'tyr', 'Y')] = 'C1=CC(=CC=C1C[C@@H](C(=O)O)N)O'
aminos[('valine', 'val', 'V')] = 'CC(C)[C@@H](C(=O)O)N'

just_name = {k[0]:aminos[k] for k in aminos.keys()}
just_symbols = {k[2]:aminos[k] for k in aminos.keys()}

def str_to_class(classname):
    return getattr(sys.modules[__name__], classname)

symbol_class = {}
aas = []
for full, triple, single in aminos.keys():
    if ' ' in full:
        f, s = full.split()
        c_name = f[0].upper() + f[1:] + '_' + s[0].upper()+s[1:]
    else:
        c_name = full[0].upper() + full[1:]

    symbol_class[single] = str_to_class(c_name)

class Protein(mb.Compound):
    def __init__(self):
        super().__init__()
        
    def build(self, seq, n = 1):
        start_mol = symbol_class[seq[0]]()
        N_terminal = self.remove_C_terminal(start_mol)
        previous_acid = N_terminal
        self.add(previous_acid,label='N-terminal')
        
        for i, letter in enumerate(seq[1:-1]):
            current_acid = symbol_class[letter]()
            current_acid = self.remove_both_terminals(current_acid)
            N_port = mb.Port(anchor=current_acid.amine, 
                             orientation=[0, 1, -1], 
                             separation=0.1)
            C_port = mb.Port(anchor=current_acid.carboxyl, 
                             orientation=[1, 0, 0], 
                             separation=0.1)
            current_acid.add(N_port, label='head')
            current_acid.add(C_port, label='tail')
            # current_acid['head'].rotate(around=[-10,0,1], theta=(i+1)*np.pi/4)
            mb.force_overlap(move_this=current_acid,
                             from_positions=current_acid['head'],
                             to_positions=previous_acid['tail'])
            self.add(current_acid, label='aa_'+str(i+1))
            previous_acid = current_acid
            
        final_mol = symbol_class[seq[-1]]()
        C_terminal = self.remove_N_terminal(final_mol)

        mb.force_overlap(move_this=C_terminal,
                         from_positions=C_terminal['head'],
                         to_positions=previous_acid['tail'])
        self.add(C_terminal,label='C-terminal')

    def remove_C_terminal(self, mol):
        amine_h, carboxyl_o = mol.indices
        mol.remove(carboxyl_o)
        mol.add(mb.Port(anchor=mol.carboxyl,
                               orientation=[0, 1, 0], 
                               separation=0.1), 
                       label='tail')
        mol.spin(around=[0,1,0], theta=np.pi/2)
        return mol

    def remove_N_terminal(self, mol):
        amine_h, carboxyl_o = mol.indices
        mol.remove(amine_h[0])
        mol.add(mb.Port(anchor=mol.amine, 
                               orientation=[-1, .7, -.2], 
                               separation=0.05), 
                       label='head')
        return mol
    
    def remove_both_terminals(self, mol):
        amine_h, carboxyl_o = mol.indices
        mol.remove(amine_h[0])
        mol.remove(carboxyl_o)
        return mol

if __name__ == '__main__':
    for i in np.random.choice(list(just_symbols.keys()), size=(15,15)):
        chain = Protein()
        chain.build(''.join(list(i)))
        print(''.join(list(i))+ ' Works')
        # chain.visualize(backend='nglview',show_ports=True)
