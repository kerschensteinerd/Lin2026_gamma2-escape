"""Physically remove A-E/I from the existing PDF, retaining F-H exactly.

No hidden obsolete model panels are retained. Experimental pixels must match
at 288 dpi in independent PDFium renders before this base is accepted.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pypdfium2 as pdfium
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tmp/fig4_v9/python_deps'))
import pymupdf as fitz
SRC=ROOT/'figures/main/PDF/Fig4_v8.pdf';OUT=ROOT/'figures/main/assets/Fig4_v9'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    source_sha=hashlib.sha256(SRC.read_bytes()).hexdigest()
    doc=fitz.open(SRC);page=doc[0]
    assert len(doc)==1 and abs(page.rect.width-493.228)<.01 and abs(page.rect.height-257.509)<.01
    for rect in [(0,0,246,257.509),(407,165,493.228,257.509)]:
        page.add_redact_annot(fitz.Rect(rect),fill=False,cross_out=False)
    page.apply_redactions(images=1,graphics=2,text=0)
    # Illustrator otherwise ignores modified PDF operators and opens its stale
    # embedded native payload. Remove that original-model payload and preview.
    for key in ['PieceInfo','Thumb']:
        doc.xref_set_key(page.xref,key,'null')
    dest=OUT/'Fig4_experimental_F-H.pdf';doc.save(dest,garbage=4,deflate=True);doc.close()
    scale=4
    a=np.array(pdfium.PdfDocument(SRC)[0].render(scale=scale).to_pil().convert('RGB'))
    b=np.array(pdfium.PdfDocument(dest)[0].render(scale=scale).to_pil().convert('RGB'))
    regions={'F':(250,0,493.228,165),'G-H':(250,165,407,257.509)};checks={}
    for key,box in regions.items():
        x0,y0,x1,y1=[int(v*scale) for v in box]
        aa=a[y0:y1,x0:x1];bb=b[y0:y1,x0:x1]
        checks[key]={'pixel_identical':bool(np.array_equal(aa,bb)),'different_pixels':int(np.any(aa!=bb,axis=2).sum()),'max_channel_difference':int(np.abs(aa.astype(int)-bb.astype(int)).max()),'crop_points':box}
        assert checks[key]['pixel_identical'],(key,checks[key])
    text=fitz.open(dest)[0].get_text()
    assert 'GABA removed' not in text and 'Spatial attenuation' not in text and 'OFF comparison' not in text
    assert 'OFF sublamina' in text and 'ON sublamina' in text
    clean=fitz.open(dest)
    assert clean.xref_get_key(clean[0].xref,'PieceInfo')[0]=='null'
    clean.close()
    assert source_sha==hashlib.sha256(SRC.read_bytes()).hexdigest()
    (OUT/'experimental_preservation.json').write_text(json.dumps({'source':str(SRC.relative_to(ROOT)),'source_sha256':source_sha,'output':str(dest.relative_to(ROOT)),'output_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'method':'True PDF redaction removes painting operators for A-E and I; no masking or clipped obsolete model objects.','pymupdf_version':fitz.__version__,'render_engine':'PDFium','resolution_dpi':288,'checks':checks,'obsolete_model_text_removed':True,'source_unchanged':True},indent=2)+'\n')
    print(json.dumps(checks,indent=2))
if __name__=='__main__':main()
