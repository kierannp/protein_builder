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