#!/usr/bin/env python3
"""Assemble experimental and model Figure 3 artwork at native size."""
from pathlib import Path
import hashlib
import json
from pypdf import PdfReader, PdfWriter, Transformation
PROJECT = Path(__file__).resolve().parents[3]
OUT = PROJECT / 'figures/main/assets/Fig3_v3'
DEST = PROJECT / 'figures/main/PDF/Fig3_v3.pdf'


def main():
    experimental_path = OUT/'Fig3_experimental_panels.pdf'
    model_path = OUT/'Fig3_modeling_panels.pdf'
    ep = PdfReader(experimental_path).pages[0]
    mp = PdfReader(model_path).pages[0]
    ew, eh = map(float, (ep.mediabox.width, ep.mediabox.height))
    mw, mh = map(float, (mp.mediabox.width, mp.mediabox.height))
    assert abs(ew-mw) < 1e-6
    margin, gap = 12., 10.
    width, height = mw+2*margin, mh+eh+gap+2*margin
    writer = PdfWriter()
    dest = writer.add_blank_page(width=width, height=height)
    dest.merge_transformed_page(ep, Transformation().translate(margin, margin+mh+gap))
    dest.merge_transformed_page(mp, Transformation().translate(margin, margin))
    dest.compress_content_streams()
    writer.add_metadata({'/Title': 'Figure 3: VG3 spatial integration and calcium contrast gain',
                         '/Subject': 'Imaging and compartmental modeling; common vector figure style',
                         '/Author': 'Kerschensteiner laboratory'})
    with DEST.open('wb') as stream:
        writer.write(stream)
    (OUT/'assembly_manifest.json').write_text(json.dumps(dict(
        experimental_sha256=hashlib.sha256(experimental_path.read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest(),
        output_sha256=hashlib.sha256(DEST.read_bytes()).hexdigest(),
        original_source_ai_sha256=hashlib.sha256((PROJECT/'figures/main/AI/Fig3_v1.ai').read_bytes()).hexdigest(),
        output_page_points=[width,height], experimental_offset_points=[margin,margin+mh+gap],
        model_offset_points=[margin,margin],
        transform='Translation only for both native-size panel blocks; no font or marker scaling'), indent=2)+'\n')
    print(DEST)


if __name__ == '__main__':
    main()
