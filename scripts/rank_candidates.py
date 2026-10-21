import csv

with open(r"C:\Users\user\.gemini\antigravity\scratch\glp1r_agonist_design\glp1r_candidate_library_40.csv", "r", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

# Multi-parameter triage scoring function simulating the Claude/Adaptyv automated filter:
# 1. Agonist trigger integrity (+15 if free N-terminal His)
# 2. Solubility / GRAVY (+10 if in sweet spot -0.7 to -0.3; -30 if >0)
# 3. pI distance from 7.4 (+10 if outside 6.5-8.0 precipitation band; -25 if in band)
# 4. Length feasibility (ideal 55-90 aa for high-yield cell-free synthesis; -10 if >130)
# 5. Cysteine safety (0 unpaired cysteines; -25 if non-canonical Cys present)
# 6. Architectural diversity bonus (+5 for top representative of each architecture)

scored = []
for r in rows:
    score = 70.0
    length = int(r['length'])
    if 55 <= length <= 90:
        score += 10
    elif length > 130:
        score -= 5
    
    gravy = float(r['gravy'])
    if gravy > 0:
        score -= 30
    elif -0.75 <= gravy <= -0.35:
        score += 10
    
    pi = float(r['pI'])
    if 6.8 <= pi <= 7.8:
        score -= 25
    else:
        score += 10
        
    cys = int(r['cys_count'])
    if 'VHH' not in r['architecture'] and cys > 0:
        score -= 25
    elif 'VHH' in r['architecture'] and cys == 2:
        score += 5 # canonical framework disulfide
        
    if r['free_n_his'] == 'True':
        score += 15
        
    scored.append((score, r))

# Sort by score descending
scored.sort(key=lambda x: x[0], reverse=True)

print("=== TOP 12 PRIORITIZED 'HERO' CANDIDATES FOR ADAPTYV BIO-FOUNDRY ===")
for i, (s, r) in enumerate(scored[:12], 1):
    print(f"{i:2d}. {r['id']:14s} | Score: {s:4.0f} | Arch: {r['architecture']:26s} | Len: {r['length']:3s} | pI: {r['pI']:4s} | GRAVY: {r['gravy']:6s} | {r['notes']}")

# Export ranked version
with open(r"C:\Users\user\.gemini\antigravity\scratch\glp1r_agonist_design\glp1r_candidate_library_ranked.csv", "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=['rank', 'triage_score', 'id', 'architecture', 'length', 'pI', 'net_charge_pH7.4', 'gravy', 'cys_count', 'free_n_his', 'notes', 'sequence'])
    writer.writeheader()
    for rank, (score, r) in enumerate(scored, 1):
        row_out = {'rank': rank, 'triage_score': int(score)}
        row_out.update(r)
        writer.writerow(row_out)

print("\nSuccessfully saved ranked library to glp1r_candidate_library_ranked.csv")
