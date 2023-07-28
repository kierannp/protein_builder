import mbuild as mb

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


class Alanine(mb.Compound):
    def __init__(self):
        super(Alanine, self).__init__()
        alanine = mb.load(just_name['alanine'],smiles=True)
        self.add(alanine)
        self.name = "Alanine"
        self.amine = alanine[5]
        self.carboxyl = alanine[2]
        amine_h = [alanine[11], alanine[12]]
        carboxyl_o = alanine[3]
        self.indices = [amine_h, carboxyl_o]
        
class Arginine(mb.Compound):
    def __init__(self):
        super(Arginine, self).__init__()
        arginine = mb.load(just_name['arginine'],smiles=True)
        self.add(arginine)
        self.name = "Arginine"
        self.amine = arginine[8]
        self.carboxyl = arginine[9]
        amine_h = [arginine[18], arginine[19]]
        carboxyl_o = arginine[4]
        self.indices = [amine_h, carboxyl_o]
        
class Asparagine(mb.Compound):
    def __init__(self):
        super(Asparagine, self).__init__()
        asparagine = mb.load(just_name['asparagine'],smiles=True)
        self.add(asparagine)
        self.name = "Asparagine"
        self.amine = asparagine[5]
        self.carboxyl = asparagine[2]
        amine_h = [asparagine[13], asparagine[14]]
        carboxyl_o = asparagine[3]
        self.indices = [amine_h, carboxyl_o]
        

class Aspartic_Acid(mb.Compound):
    def __init__(self):
        super(Aspartic_Acid, self).__init__()
        aspartic_acid = mb.load(just_name['aspartic acid'],smiles=True)
        self.add(aspartic_acid)
        self.name = "Aspartic_Acid"
        self.amine = aspartic_acid[5]
        self.carboxyl = aspartic_acid[6]
        amine_h = [aspartic_acid[13], aspartic_acid[14]]
        carboxyl_o = aspartic_acid[7]
        self.indices = [amine_h, carboxyl_o]       

class Cysteine(mb.Compound):
    def __init__(self):
        super(Cysteine, self).__init__()
        cysteine = mb.load(just_name['cysteine'],smiles=True)
        self.add(cysteine)
        self.name = "Cysteine"
        self.amine = cysteine[5]
        self.carboxyl = cysteine[2]
        amine_h = [cysteine[11],cysteine[12]]
        carboxyl_o = cysteine[3]
        self.indices = [amine_h, carboxyl_o]

class Glutamine(mb.Compound):
    def __init__(self):
        super(Glutamine, self).__init__()
        glutamine = mb.load(just_name['glutamine'],smiles=True)
        self.add(glutamine)
        self.name = "Glutamine"
        self.amine = glutamine[9]
        self.carboxyl = glutamine[6]
        amine_h = [glutamine[18],glutamine[19]]
        carboxyl_o = glutamine[7]
        self.indices = [amine_h, carboxyl_o]

class Glutamic_Acid(mb.Compound):
    def __init__(self):
        super(Glutamic_Acid, self).__init__()
        glutamic_acid = mb.load(just_name['glutamic acid'],smiles=True)
        self.add(glutamic_acid)
        self.name = "Glutamic_Acid"
        self.amine = glutamic_acid[9]
        self.carboxyl = glutamic_acid[6]
        amine_h = [glutamic_acid[17], glutamic_acid[18]]
        carboxyl_o = glutamic_acid[7]
        self.indices = [amine_h, carboxyl_o]

        
class Glycine(mb.Compound):
    def __init__(self):
        super(Glycine, self).__init__()
        glycine = mb.load(just_name['glycine'],smiles=True)
        self.add(glycine)
        self.name = "Glycine"
        self.amine = glycine[4]
        self.carboxyl = glycine[1]
        amine_h = [glycine[8], glycine[9]]
        carboxyl_o = [glycine[3], glycine[7]]
        self.indices = [amine_h, carboxyl_o]
 
class Histidine(mb.Compound):
    def __init__(self):
        super(Histidine, self).__init__()
        histidine = mb.load(just_name['histidine'],smiles=True)
        self.add(histidine)
        self.name = "Histidine"
        self.amine = histidine[10]
        self.carboxyl = histidine[7]
        amine_h = [histidine[18], histidine[19]]
        carboxyl_o = histidine[8]
        self.indices = [amine_h, carboxyl_o]

class Isoleucine(mb.Compound):
    def __init__(self):
        super(Isoleucine, self).__init__()
        isoleucine = mb.load(just_name['isoleucine'],smiles=True)
        self.add(isoleucine)
        self.name = "Isoleucine"
        self.amine = isoleucine[8]
        self.carboxyl = isoleucine[5]
        amine_h = [isoleucine[20], isoleucine[21]]
        carboxyl_o = isoleucine[6]
        self.indices = [amine_h, carboxyl_o]

class Leucine(mb.Compound):
    def __init__(self):
        super(Leucine, self).__init__()
        leucine = mb.load(just_name['leucine'],smiles=True)
        self.add(leucine)
        self.name = "Leucine"
        self.amine = leucine[8]
        self.carboxyl = leucine[5]
        amine_h = [leucine[20], leucine[21]]
        carboxyl_o = leucine[6]
        self.indices = [amine_h, carboxyl_o]

class Lysine(mb.Compound):
    def __init__(self):
        super(Lysine, self).__init__()
        lysine = mb.load(just_name['lysine'],smiles=True)
        self.add(lysine)
        self.name = "Lysine"
        self.amine = lysine[9]
        self.carboxyl = lysine[6]
        amine_h = [lysine[22], lysine[23]]
        carboxyl_o = lysine[7]
        self.indices = [amine_h, carboxyl_o]
 
class Methionine(mb.Compound):
    def __init__(self):
        super(Methionine, self).__init__()
        methionine = mb.load(just_name['methionine'],smiles=True)
        self.add(methionine)
        self.name = "Methionine"
        self.amine = methionine[8]
        self.carboxyl = methionine[5]
        amine_h = [methionine[18], methionine[19]]
        carboxyl_o = methionine[6]
        self.indices = [amine_h, carboxyl_o]
   
class Phenylalanine(mb.Compound):
    def __init__(self):
        super(Phenylalanine, self).__init__()
        phenylalanine = mb.load(just_name['phenylalanine'],smiles=True)
        self.add(phenylalanine)
        self.name = "Phenylalanine"
        self.amine = phenylalanine[11]
        self.carboxyl = phenylalanine[8]
        amine_h = [phenylalanine[21], phenylalanine[22]]
        carboxyl_o = phenylalanine[9]
        self.indices = [amine_h, carboxyl_o]
 
class Proline(mb.Compound):
    def __init__(self):
        super(Proline, self).__init__()
        proline = mb.load(just_name['proline'],smiles=True)
        self.add(proline)
        self.name = "Proline"
        self.amine = proline[3]
        self.carboxyl = proline[5]
        amine_h = proline[13]
        carboxyl_o = proline[6]
        self.indices = [amine_h, carboxyl_o]
  
class Serine(mb.Compound):
    def __init__(self):
        super(Serine, self).__init__()
        serine = mb.load(just_name['serine'],smiles=True)
        self.add(serine)
        self.name = "Serine"
        self.amine = serine[5]
        self.carboxyl = serine[2]
        amine_h = [serine[11], serine[12]]
        carboxyl_o = serine[3]
        self.indices = [amine_h, carboxyl_o]

class Threonine(mb.Compound):
    def __init__(self):
        super(Threonine, self).__init__()
        threonine = mb.load(just_name['threonine'],smiles=True)
        self.add(threonine)
        self.name = "Threonine"
        self.amine = threonine[6]
        self.carboxyl = threonine[3]
        amine_h = [threonine[14], threonine[15]]
        carboxyl_o = threonine[4]
        self.indices = [amine_h, carboxyl_o]
 
class Tryptophan(mb.Compound):
    def __init__(self):
        super(Tryptophan, self).__init__()
        tryptophan = mb.load(just_name['tryptophan'],smiles=True)
        self.add(tryptophan)
        self.name = "Tryptophan"
        self.amine = tryptophan[14]
        self.carboxyl = tryptophan[11]
        amine_h = [tryptophan[25], tryptophan[26]]
        carboxyl_o = tryptophan[12]
        self.indices = [amine_h, carboxyl_o]

class Tyrosine(mb.Compound):
    def __init__(self):
        super(Tyrosine, self).__init__()
        tyrosine = mb.load(just_name['tyrosine'],smiles=True)
        self.add(tyrosine)
        self.name = "Tyrosine"
        self.amine = tyrosine[11]
        self.carboxyl = tyrosine[8]
        amine_h = [tyrosine[21], tyrosine[22]]
        carboxyl_o = tyrosine[9]
        self.indices = [amine_h, carboxyl_o]

class Valine(mb.Compound):
    def __init__(self):
        super(Valine, self).__init__()
        valine = mb.load(just_name['valine'],smiles=True)
        self.add(valine)
        self.name = "Valine"
        self.amine = valine[7]
        self.carboxyl = valine[4]
        amine_h = [valine[17], valine[18]]
        carboxyl_o = valine[5]
        self.indices = [amine_h, carboxyl_o]

if __name__ == '__main__':
    import sys
    def str_to_class(classname):
        return getattr(sys.modules[__name__], classname)
    aa = [i[0].upper()+i[1:] for i in list(just_name.keys())]
    aas = []
    for i in aa:
        if ' ' in i:
            f, s = i.split()
            aas.append(f + '_' + s[0].upper()+s[1:])
        else:
            aas.append(i)
    aa = aas

    for n in aa:
        try:
            mol = str_to_class(n)()
            if n =='Proline':
                if mol.amine.element.symbol + ' ' +  mol.carboxyl.element.symbol != 'N C':
                    raise Exception('error with initialize of '+ n)
                if mol.indices[0].element.symbol + ' ' + mol.indices[1].element.symbol != 'H O':
                    raise Exception('error with initialize of '+ n)
                continue
            if mol.amine.element.symbol + ' ' +  mol.carboxyl.element.symbol != 'N C':
                raise Exception('error with initialize of '+ n)
            if mol.indices[0][0].element.symbol + ' ' + mol.indices[0][1].element.symbol + ' ' + mol.indices[1].element.symbol != 'H H O':
                raise Exception('error with initialize of '+ n)
            print(n + ' Works')
        except:
            raise Exception('error with initialize of '+ n)
