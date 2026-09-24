"""Change only the visible Figure 7G compartmental-model marker and label."""
from pathlib import Path
import json,hashlib,math
from pypdf import PdfReader,PdfWriter
from pypdf.generic import DecodedStreamObject,NameObject
ROOT=Path(__file__).resolve().parents[3]
SRC=ROOT/'output/pdf/Figure7_v7.pdf';DEST=ROOT/'output/pdf/Figure7_v8.pdf';OUT=ROOT/'figures/main/assets/Fig7_v8'
MODEL=ROOT/'analysis/fig4_results_v9_20260916/model_shift_v9.json'

def main():
    assert not DEST.exists(),'New PDF version exists; refusing to overwrite.'
    protected=[SRC,ROOT/'figures/main/AI/Fig7_v7.ai',ROOT/'Data_Figures 9.xlsx',ROOT/'data_vos.xlsx']
    before={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
    model=json.loads(MODEL.read_text());row=next(r for r in model['sensitivity'] if r['retained_effective_inhibition_fraction']==.25)
    value=row['descriptive_marker_arithmetic_mean'];assert row['n_valid_readout_contrast_matches']==8
    x=92+math.log(value/.8)/math.log(13/.8)*140
    old=(ROOT/'analysis/fig4_inhibition25_20260916/figure7_visible_marker_original_stream.txt').read_bytes()
    new=old.replace(b'152.9057 63 4 4 re',f'{x-2:.8f} 63 4 4 re'.encode()).replace(b'160.9057 62.8 cm',f'{x+6:.8f} 62.8 cm'.encode()).replace(b'(\\0012\\0568)',f'(\\001{value:.1f})'.replace('.',r'\056').encode())
    reader=PdfReader(SRC);data=reader.pages[0].get_contents().get_data();assert data.count(old)==1
    updated=data.replace(old,new);assert updated.replace(new,old)==data
    writer=PdfWriter();writer.clone_document_from_reader(reader)
    stream=DecodedStreamObject();stream.set_data(updated);writer.pages[0][NameObject('/Contents')]=stream;writer.pages[0].compress_content_streams()
    writer.add_metadata({'/Title':'Figure 7 v8: updated fixed-parameter model marker','/Subject':'Only visible model marker and its label changed; behavioral data, fitted estimates, intervals, layout, and OMR observations preserved.'})
    with DEST.open('wb') as f:writer.write(f)
    result=PdfReader(DEST).pages[0].get_contents().get_data();assert result==updated
    for p,h in before.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'marker_manifest.json').write_text(json.dumps({'sources':before,'model_source':str(MODEL.relative_to(ROOT)),'model_source_sha256':hashlib.sha256(MODEL.read_bytes()).hexdigest(),'new_pdf':str(DEST.relative_to(ROOT)),'new_pdf_sha256':hashlib.sha256(DEST.read_bytes()).hexdigest(),'old_marker_display':2.8,'new_marker_arithmetic_mean':value,'new_marker_display_label':f'≈{value:.1f}','new_marker_center_points':[x,470],'new_marker_n_matches':8,'marker_has_uncertainty_interval':False,'only_visible_marker_block_changed':True,'remaining_decoded_content_bytes_identical':True,'old_block':old.decode(),'new_block':new.decode(),'behavioral_data_fits_and_layout_unchanged':True,'sources_unchanged':True},indent=2)+'\n')
    print(DEST);print('marker',value,'x',x,'label',f'≈{value:.1f}')
if __name__=='__main__':main()
