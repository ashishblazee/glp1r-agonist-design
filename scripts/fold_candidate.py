"""
Interactive 3D Structure Prediction script for undergraduates & researchers.
Sends a candidate sequence to the Meta ESMFold API, downloads the 3D atomic coordinates (PDB),
and analyzes the per-residue pLDDT confidence scores.
"""

import urllib.request
import os
import sys

def fold_sequence(candidate_id, sequence, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    pdb_path = os.path.join(output_dir, f"{candidate_id}.pdb")
    
    print(f"\n--- Folding {candidate_id} ({len(sequence)} amino acids) ---")
    print(f"Sequence: {sequence}")
    print("Contacting Meta ESMFold API (this takes ~4-8 seconds)...")
    
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    req = urllib.request.Request(
        url,
        data=sequence.encode("utf-8"),
        headers={"Content-Type": "text/plain"}
    )
    
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            pdb_content = resp.read().decode("utf-8")
            
        with open(pdb_path, "w", encoding="utf-8") as f:
            f.write(pdb_content)
            
        print(f"Success! Saved 3D coordinate file to: {pdb_path}")
        
        # Parse pLDDT from PDB (stored in the B-factor column, characters 60-66)
        plddt_scores = []
        for line in pdb_content.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                score = float(line[60:66].strip())
                plddt_scores.append(score)
                
        if plddt_scores:
            avg_plddt = sum(plddt_scores) / len(plddt_scores)
            high_conf = sum(1 for s in plddt_scores if s >= 70) / len(plddt_scores) * 100
            print(f"Average pLDDT Confidence: {avg_plddt:.1f} / 100")
            print(f"Residues with high confidence (>=70): {high_conf:.1f}%")
            
            # Breakdown: N-term warhead, linker, scaffold
            warhead_plddt = sum(plddt_scores[:15]) / 15
            print(f"  - N-terminal Warhead (first 15 aa) pLDDT: {warhead_plddt:.1f}")
            if len(plddt_scores) > 25:
                scaffold_plddt = sum(plddt_scores[25:]) / len(plddt_scores[25:])
                print(f"  - De Novo Scaffold Body pLDDT: {scaffold_plddt:.1f}")
                
        return pdb_path
        
    except Exception as e:
        print(f"Error during folding: {e}")
        return None

if __name__ == "__main__":
    # Test on our Rank 1 candidate: GLP1R_DES_A01
    rank1_id = "GLP1R_DES_A01"
    rank1_seq = "HAEGTFTSDVSSYLEGGGGSGGGGSGQAAKEFIAWLVKGRGGSGGEEAKKIAEEAAKKIAEEAAKKGE"
    
    output_directory = r"C:\Users\user\.gemini\antigravity\scratch\glp1r_agonist_design\predicted_structures"
    fold_sequence(rank1_id, rank1_seq, output_directory)
