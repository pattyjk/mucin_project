#!/usr/bin/env python
"""Check the OD600 sheet (screenshot) against mucin_media_map.txt.
Usage: python od600_map_check.py screenshot.png mucin_media_map.txt
Needs: pillow, numpy, pandas, pytesseract + tesseract binary.
The table grid is detected from the image and each cell is read with OCR.
OCR misread four 'D5' cells and four 'A11' cells as 'D'/'A1'; these were fixed by eye
(rows identified by OD value + timepoint) -- the FIXES dict below reproduces that."""
import sys, numpy as np, pandas as pd, pytesseract
from PIL import Image
from concurrent.futures import ThreadPoolExecutor
im = Image.open(sys.argv[1]); g = np.array(im.convert("L")).astype(int)
def lines(frac):
    ix = np.where(frac > 0.6)[0]; out = []; s = p = ix[0]
    for i in ix[1:]:
        if i > p + 1: out.append((s + p) / 2); s = i
        p = i
    out.append((s + p) / 2); return out
H = lines((g < 230).mean(1)); V = lines((g < 230).mean(0))
wl = {0: "0123456789.", 1: "0123456789", 2: "ABCDEFGH0123456789", 3: "YesNo", 4: None}
def ocr(rc):
    r, c = rc; box = (int(V[c]) + 2, int(H[r]) + 2, int(V[c + 1]) - 1, int(H[r + 1]) - 1)
    cell = im.crop(box).convert("L").resize(((box[2] - box[0]) * 5, (box[3] - box[1]) * 5), Image.LANCZOS)
    cell = Image.fromarray(np.pad(np.array(cell), 20, constant_values=255))
    cfg = "--psm 7" + (f" -c tessedit_char_whitelist={wl[c]}" if wl[c] else "")
    return r, c, pytesseract.image_to_string(cell, config=cfg).strip()
with ThreadPoolExecutor(8) as ex: res = list(ex.map(ocr, [(r, c) for r in range(1, len(H) - 1) for c in range(5)]))
tab = pd.DataFrame(index=range(1, len(H) - 1), columns=range(5))
for r, c, s in res: tab.loc[r, c] = s
tab.columns = ["OD600", "TimePoint", "Well", "Mucin", "AntiBiotic"]
tab = tab[tab.Well.str.contains(r"\d") | tab.Well.isin(["D", "A"])].reset_index(drop=True)
FIXES = {("0.308","1"):"D5", ("0.362","4"):"D5", ("0.545","3"):"D5", ("0.64","2"):"D5"}       # Neomycin, no mucin
for i, r in tab.iterrows():
    if r.Well == "D" and (r.OD600, r.TimePoint) in FIXES: tab.loc[i, "Well"] = FIXES[(r.OD600, r.TimePoint)]
    if r.Well == "A1" and r.AntiBiotic.startswith("Blank"): tab.loc[i, "Well"] = "A11"            # A11 is a BlankWell
fix = {"Gentamyci": "Gentamycin", "Chlorampt": "Chloramphenicol", "Trimethop:": "Trimethoprim", "Streptomy": "Streptomycin", "Sulfamethz": "Sulfamethazine"}
tab["AntiBiotic"] = tab.AntiBiotic.replace(fix); tab["OD600"] = tab.OD600.astype(float); tab["TimePoint"] = tab.TimePoint.astype(int)
m = pd.read_csv(sys.argv[2], sep="\t")
od = tab.groupby("Well")[["Mucin", "AntiBiotic"]].first()
j = m.set_index("Well")[["SampleID", "Mucin", "Antibiotic"]].join(od.rename(columns={"Mucin": "Mucin_OD"}), how="outer")
mis = j[(j.Mucin != j.Mucin_OD) | (j.Antibiotic != j.AntiBiotic)]
print("rows", len(tab), "wells", tab.Well.nunique(), "mismatches", len(mis)); print(mis)
tab.to_csv("od600_transcribed.csv", index=False)
