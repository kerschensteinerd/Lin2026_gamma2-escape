"""Assemble new model A-E/I onto the exactly preserved current F-H base."""
from pathlib import Path
from copy import deepcopy
import json,hashlib
from pypdf import PdfReader,PdfWriter,Transformation
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'figures/main/assets/Fig4_v9';DEST=ROOT/'figures/main/PDF/Fig4_v9.pdf'
DATA=ROOT/'vg3_cid3_rf_simulation-main/outputs/fig4_v9_retained_inhibition_20260916'

def main():
    checks=json.loads((DATA/'validation.json').read_text());assert all(x['passed'] for x in checks)
    assert not DEST.exists(),'New version already exists; refusing to overwrite.'
    baseline=ROOT/'figures/main/PDF/Fig4_v8.pdf'
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [baseline,ROOT/'figures/main/AI/Fig4_v8.ai',ROOT/'Data_Figures 9.xlsx']}
    layout_source=ROOT/'figures/main/assets/Fig3_v8/assembly_manifest.json'
    historical_layout=json.loads(layout_source.read_text())
    layout={k:v for k,v in historical_layout.items() if k not in ['experimental_sha256','model_sha256','output_sha256','source_blocks','layout_target','transform']}
    layout['layout_provenance']=str(layout_source.relative_to(ROOT))
    layout['layout_target']='Existing Figure 4 v8 geometry retained: new model A-E at left, original experimental F-H and new model I at right.'
    layout['transform']='Inherited uniform scale and translations only; new model panels are generated from the nominal 25% retained-inhibition outputs.'
    scale=layout['uniform_scale'];model_path=OUT/'Fig4_modeling_panels.pdf';ipath=OUT/'Fig4_panel_I.pdf';base=OUT/'Fig4_experimental_F-H.pdf'
    writer=PdfWriter();page=writer.add_page(PdfReader(base).pages[0])
    for key in ['/PieceInfo','/Thumb']:
        if key in page:del page[key]
    model=PdfReader(model_path).pages[0]
    assert abs(float(model.mediabox.width)-layout['source_model_page_points'][0])<1e-6
    page.merge_transformed_page(model,Transformation().scale(scale).translate(*layout['model_offset_points']))
    panel_i=deepcopy(PdfReader(ipath).pages[0]);width=float(panel_i.mediabox.width)
    panel_i.cropbox.lower_left=(2*width/3,0);panel_i.cropbox.upper_right=(width,130)
    page.merge_transformed_page(panel_i,Transformation().scale(scale).translate(*layout['experimental_offset_points']))
    page.compress_content_streams()
    writer.add_metadata({'/Title':'Figure 4 v9: dendritic calcium with 25% effective inhibition retained','/Subject':'Updated model A-E/I with fixed prior calcium calibration; experimental F-H preserved exactly.','/Author':'Kerschensteiner laboratory'})
    with DEST.open('wb') as f:writer.write(f)
    for path,digest in hashes.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    files=[base,model_path,ipath,DEST]
    manifest={'sources':hashes,'outputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},'model_source':str((DATA/'nominal25').relative_to(ROOT)),'model_checks_passed':len(checks),'model_calcium_refitted':False,'retained_effective_inhibition':.25,'layout':layout,'page_points':[float(page.mediabox.width),float(page.mediabox.height)],'F-H_preservation':'Original current PDF F-H retained via physical removal of obsolete A-E/I; new model content merged at original layout transforms.','model_I_label':'Lines: fixed-parameter model','condition_label':'25% retained'}
    (OUT/'assembly_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(DEST)
if __name__=='__main__':main()
