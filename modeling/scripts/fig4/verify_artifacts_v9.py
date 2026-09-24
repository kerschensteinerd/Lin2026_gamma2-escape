"""Verify final figures, native AI editability, preserved data, and PDF appearance."""
from pathlib import Path
import hashlib,json
import numpy as np,pypdfium2 as pdfium
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[3]
TMP=ROOT/'tmp/fig4_v9'

def render(path,scale=4):
    d=pdfium.PdfDocument(path);assert len(d)==1
    return np.asarray(d[0].render(scale=scale).to_pil().convert('RGB'))
def diff(a,b):
    assert a.shape==b.shape
    delta=np.abs(a.astype(int)-b.astype(int))
    return {'mean_absolute_RGB_difference':float(delta.mean()),'fraction_pixels_any_channel_gt20':float(np.mean(np.any(delta>20,axis=2))),'different_pixels':int(np.any(delta!=0,axis=2).sum()),'pixel_identical':bool(np.array_equal(a,b))}
def crop(a,box,s=4):
    x0,y0,x1,y1=[int(v*s) for v in box];return a[y0:y1,x0:x1]
result={}
for fig,old,pdf,ai,out in [('Figure4','figures/main/PDF/Fig4_v8.pdf','figures/main/PDF/Fig4_v9.pdf','figures/main/AI/Fig4_v9.ai','figures/main/assets/Fig4_v9'),('Figure7','output/pdf/Figure7_v7.pdf','output/pdf/Figure7_v8.pdf','figures/main/AI/Fig7_v8.ai','figures/main/assets/Fig7_v8')]:
    a,b,native=[render(ROOT/p) for p in [old,pdf,ai]]
    qa=json.loads((ROOT/out/'native_illustrator_qa.json').read_text());assert qa['native_save_reopened'] and qa['after']['external_links']==0
    assert b.shape==native.shape
    nd=diff(b,native);assert nd['mean_absolute_RGB_difference']<.6 and nd['fraction_pixels_any_channel_gt20']<.02
    entry={'native_vs_PDF':nd,'native_editability':qa['after'],'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in [old,pdf,ai]},'render_dpi':288,'single_page':True}
    if fig=='Figure4':
        assert '/PieceInfo' not in PdfReader(ROOT/pdf).pages[0], 'Stale native Illustrator payload could override PDF content on import.'
        for path in [pdf,ai]:
            txt=PdfReader(ROOT/path).pages[0].extract_text()
            assert '25% retained' in txt and 'Lines: fixed-parameter model' in txt
            assert 'GABA removed' not in txt and 'Lines: calibrated model' not in txt
        entry['F-H_preservation']={}
        entry['native_F-H_preservation']={}
        for k,box in {'F':[250,0,493.228,165],'G-H':[250,165,407,257.509]}.items():
            dd=diff(crop(a,box),crop(b,box));assert dd['pixel_identical'];entry['F-H_preservation'][k]=dd
            dd_native=diff(crop(a,box),crop(native,box));assert dd_native['pixel_identical'];entry['native_F-H_preservation'][k]=dd_native
        entry['updated_model_labels_verified_in_PDF_and_native']=True
        checks=json.loads((ROOT/'vg3_cid3_rf_simulation-main/outputs/fig4_v9_retained_inhibition_20260916/validation.json').read_text());assert all(q['passed'] for q in checks)
        entry['model_validation_checks_passed']=len(checks)
        assembly=json.loads((ROOT/out/'assembly_manifest.json').read_text())
        for p,h in assembly['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    else:
        mask=np.ones(a.shape[:2],dtype=bool);mask[int(464*4):int(478*4),int(148*4):int(193*4)]=False
        assert np.array_equal(a[mask],b[mask]);entry['outside_model_marker_pixel_identical']=True;entry['allowed_change_box_points']=[148,464,193,478]
        marker=json.loads((ROOT/out/'marker_manifest.json').read_text())
        oldstream=PdfReader(ROOT/old).pages[0].get_contents().get_data();newstream=PdfReader(ROOT/pdf).pages[0].get_contents().get_data()
        assert newstream.replace(marker['new_block'].encode(),marker['old_block'].encode())==oldstream
        entry['all_other_decoded_content_bytes_identical']=True
        for p,h in marker['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    # Editable Illustrator private data, not a PDF renamed to AI.
    raw=(ROOT/ai).read_bytes();assert b'Adobe Illustrator' in raw and b'/Illustrator' in raw
    entry['Illustrator_private_data_present']=True
    pdfium.PdfDocument(ROOT/ai)[0].render(scale=3).to_pil().save(TMP/f'{fig}_native.png')
    (ROOT/out/'final_verification.json').write_text(json.dumps(entry,indent=2)+'\n');result[fig]=entry
print(json.dumps({k:{'native_vs_PDF':v['native_vs_PDF'],'native_editability':v['native_editability']} for k,v in result.items()},indent=2))
