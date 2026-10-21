import sys
from collections import defaultdict

cif_path = r"C:\Users\user\.gemini\antigravity\brain\306d5cd7-295e-4c4a-a823-72f75a5942be\6X18\6X18.cif"

chains = defaultdict(lambda: {'auth_seqs': set(), 'comp_ids': {}, 'label_asym': set(), 'entity_id': set()})

with open(cif_path, 'r', encoding='utf-8') as f:
    for line in f:
        if line.startswith('ATOM '):
            parts = line.split()
            label_asym = parts[6]
            entity_id = parts[7]
            auth_seq = int(parts[16])
            auth_comp = parts[17]
            auth_asym = parts[18]
            chains[auth_asym]['auth_seqs'].add(auth_seq)
            chains[auth_asym]['comp_ids'][auth_seq] = auth_comp
            chains[auth_asym]['label_asym'].add(label_asym)
            chains[auth_asym]['entity_id'].add(entity_id)

three_to_one = {
    'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',
    'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
    'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',
    'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
}

for auth_asym in sorted(chains.keys()):
    d = chains[auth_asym]
    seqs = sorted(d['auth_seqs'])
    one_letter = "".join([three_to_one.get(d['comp_ids'][s], 'X') for s in seqs])
    print(f"Chain {auth_asym}: Entity {list(d['entity_id'])}, Label {list(d['label_asym'])}, Length: {len(seqs)}, Range: {seqs[0]}-{seqs[-1]}")
    print(f"  Sequence: {one_letter}")
    print()
