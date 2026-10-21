import math

cif_path = r"C:\Users\user\.gemini\antigravity\brain\306d5cd7-295e-4c4a-a823-72f75a5942be\6X18\6X18.cif"

ca_coords = {}
with open(cif_path, 'r', encoding='utf-8') as f:
    for line in f:
        if line.startswith('ATOM '):
            parts = line.split()
            atom_name = parts[3]
            auth_comp = parts[17]
            auth_asym = parts[18]
            auth_seq = int(parts[16])
            if atom_name == 'CA':
                ca_coords[(auth_asym, auth_seq)] = (float(parts[10]), float(parts[11]), float(parts[12]), auth_comp)

def dist(c1, c2):
    return math.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2)

print("=== Key Distances in 6X18 ===")
h7 = ca_coords.get(('P', 7))
d15 = ca_coords.get(('P', 15))
v16 = ca_coords.get(('P', 16))
e21 = ca_coords.get(('P', 21))
a24 = ca_coords.get(('P', 24))
r36 = ca_coords.get(('P', 36))

print(f"His7 CA: {h7[:3]}")
print(f"Asp15 CA: {d15[:3]}")
print(f"Val16 CA: {v16[:3]}")
print(f"Glu21 CA: {e21[:3]}")
print(f"Arg36 CA: {r36[:3]}")

print(f"Distance His7 -> Asp15: {dist(h7, d15):.2f} A (N-terminal signaling stretch, 9 aa)")
print(f"Distance Asp15 -> Glu21: {dist(d15, e21):.2f} A (Hinge/stalk span, 6 aa)")
print(f"Distance Glu21 -> Arg36: {dist(e21, r36):.2f} A (ECD-bound alpha-helix, 15 aa)")
print(f"Direct end-to-end His7 -> Arg36: {dist(h7, r36):.2f} A")

# ECD binding pocket center
ecd_res = [35, 39, 68, 69, 90, 91, 121]
ecd_cas = [ca_coords[('R', r)] for r in ecd_res if ('R', r) in ca_coords]
ecd_center = (
    sum(c[0] for c in ecd_cas)/len(ecd_cas),
    sum(c[1] for c in ecd_cas)/len(ecd_cas),
    sum(c[2] for c in ecd_cas)/len(ecd_cas)
)
print(f"ECD Pocket Centroid: {ecd_center}")
print(f"Distance His7 (TM trigger) to ECD Centroid: {dist(h7, ecd_center):.2f} A")
print(f"Distance Glu21 (junction) to ECD Centroid: {dist(e21, ecd_center):.2f} A")
print(f"Distance Arg36 (ECD tip) to ECD Centroid: {dist(r36, ecd_center):.2f} A")
