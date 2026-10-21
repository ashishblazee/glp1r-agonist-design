"""
Verifies that de novo warheads have NO matching substring >= 7 aa against known GLP1R binders.
Known binders:
- Human GLP-1 (7-36): HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR
- Exendin-4 (1-39):   HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS
- Semaglutide:        HAEGTFTSDVSSYLEGQAAKEFIAWLVRGRG
- Tirzepatide:        YXEGTFTSDYSIYLDKIAQKAFVQWLIAGGPSSGAPPPS
- Human Glucagon:     HSQGTFTSDYSKYLDSRRAQDFVQWLMNT
"""

known_binders = [
    "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR",
    "HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS",
    "HAEGTFTSDVSSYLEGQAAKEFIAWLVRGRG",
    "YXEGTFTSDYSIYLDKIAQKAFVQWLIAGGPSSGAPPPS",
    "HSQGTFTSDYSKYLDSRRAQDFVQWLMNT"
]

def check_novelty(seq, threshold=7):
    # Returns max match length and offending substring
    max_len = 0
    max_sub = ""
    for ref in known_binders:
        for k in range(threshold, len(ref) + 1):
            for i in range(len(ref) - k + 1):
                sub = ref[i:i+k]
                if sub in seq:
                    if k > max_len:
                        max_len = k
                        max_sub = sub
    return max_len, max_sub

# Test our candidate de novo warheads:
# Essential conserved: His1, Glu3, Phe6, Asp9
# Conservative mutations breaking any 7-aa or 8-aa window:
test_warheads = [
    ("DW1", "HVEGSFTADIASYIE", "Val2 + Ser5 + Ala8 + Ile10 + Ala11 + Ile14"),
    ("DW2", "HAEGSFTADVALYIE", "Ser5 + Ala8 + Ala11 + Leu12 + Ile14"),
    ("DW3", "HVEGTFTADIAAYLE", "Val2 + Ala8 + Ile10 + Ala11"),
    ("DW4", "HGEGSFTADIASYLE", "Gly2 + Ser5 + Ala8 + Ile10 + Ala11"),
    ("DW5", "HVEGTFTSDIASYIE", "Val2 + Ile10 + Ala11 + Ile14"),
    ("DW6", "HAEGTFTADISAYIE", "Ala8 + Ala11 + Ile14"),
    ("DW7", "HVEGSFTSDVALYIE", "Val2 + Ser5 + Ala11 + Leu12 + Ile14"),
    ("DW8", "HGEGTFTADIAAYIE", "Gly2 + Ala8 + Ile10 + Ala11 + Ile14"),
]

print("=== CHECKING DE NOVO WARHEAD SUBSTRING MATCHES (THRESHOLD >= 7) ===")
for wid, wseq, notes in test_warheads:
    mlen, msub = check_novelty(wseq, threshold=6)
    status = "PASS (NO MATCH >= 6 aa)" if mlen < 6 else f"FAIL ({mlen} aa: {msub})"
    print(f"{wid:4s}: {wseq} | {status} | {notes}")
