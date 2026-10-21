import math
from collections import defaultdict

cif_path = r"C:\Users\user\.gemini\antigravity\brain\306d5cd7-295e-4c4a-a823-72f75a5942be\6X18\6X18.cif"

atoms = []
with open(cif_path, 'r', encoding='utf-8') as f:
    for line in f:
        if line.startswith('ATOM '):
            parts = line.split()
            atom_name = parts[3]
            auth_comp = parts[17]
            auth_asym = parts[18]
            auth_seq = int(parts[16])
            x, y, z = float(parts[10]), float(parts[11]), float(parts[12])
            atoms.append({
                'name': atom_name, 'res': auth_comp, 'chain': auth_asym,
                'seq': auth_seq, 'coord': (x, y, z)
            })

p_atoms = [a for a in atoms if a['chain'] == 'P']
r_atoms = [a for a in atoms if a['chain'] == 'R']

print(f"P atoms: {len(p_atoms)}, R atoms: {len(r_atoms)}")

def dist(c1, c2):
    return math.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2)

contacts = defaultdict(list)
for pa in p_atoms:
    for ra in r_atoms:
        d = dist(pa['coord'], ra['coord'])
        if d <= 4.0:
            contacts[(pa['seq'], pa['res'])].append((ra['seq'], ra['res'], d, pa['name'], ra['name']))

print("=== GLP-1 (Chain P) Residue Contacts with GLP1R (Chain R) (cutoff <= 4.0 A) ===")
for (p_seq, p_res) in sorted(contacts.keys()):
    r_res_list = sorted(set([(r_seq, r_res) for (r_seq, r_res, _, _, _) in contacts[(p_seq, p_res)]]))
    min_d = min([d for (_, _, d, _, _) in contacts[(p_seq, p_res)]])
    r_str = ", ".join([f"{r_res}{r_seq}" for (r_seq, r_res) in r_res_list])
    domain = "TM/Bundle (signaling)" if p_seq < 16 else ("Hinge/Stalk" if p_seq < 22 else "ECD groove (affinity)")
    print(f"GLP-1 {p_res}{p_seq:2d} [{domain:24s}]: min_dist={min_d:.2f}A -> interacts with GLP1R: {r_str}")
