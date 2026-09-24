#!/usr/bin/env python3
"""Place the v6 model and vertically aligned imaging blocks side by side."""
from pathlib import Path
import hashlib
import json
from pypdf import PdfReader, PdfWriter, Transformation

PROJECT = Path(__file__).resolve().parents[3]
OUT = PROJECT / 'figures/main/assets/Fig3_v8'
DEST = PROJECT / 'figures/main/PDF/Fig3_v8.pdf'


def main():
    model_path = OUT/'Fig3_modeling_panels.pdf'
    experimental_path = OUT/'Fig3_experimental_panels.pdf'
    model = PdfReader(model_path).pages[0]
    experimental = PdfReader(experimental_path).pages[0]
    mw, mh = map(float, (model.mediabox.width, model.mediabox.height))
    ew, eh = map(float, (experimental.mediabox.width, experimental.mediabox.height))
    margin, gutter = 8., 16.
    assert abs(mw-ew) < 1e-6
    width = 174*72/25.4
    half_width = (width-2*margin-gutter)/2
    scale = half_width/mw
    height = max(mh, eh)*scale+2*margin
    model_offset = [margin, height-margin-mh*scale]
    experimental_offset = [margin+half_width+gutter, height-margin-eh*scale]
    writer = PdfWriter()
    page = writer.add_blank_page(width=width, height=height)
    page.merge_transformed_page(model, Transformation().scale(scale).translate(*model_offset))
    page.merge_transformed_page(experimental, Transformation().scale(scale).translate(*experimental_offset))
    page.compress_content_streams()
    writer.add_metadata({
        '/Title': 'Figure 3: VG3 spatial integration and calcium contrast gain',
        '/Subject': 'v6 A-E left and F-I right; taller imaging traces align the lower plots',
        '/Author': 'Kerschensteiner laboratory',
    })
    with DEST.open('wb') as stream:
        writer.write(stream)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = {
        'experimental_sha256': sha(experimental_path),
        'model_sha256': sha(model_path),
        'output_sha256': sha(DEST),
        'output_page_points': [width, height],
        'output_page_mm': [width*25.4/72, height*25.4/72],
        'layout_target': '174-mm width: v6 A-E at left and F-I at right, with aligned top and bottom plots',
        'source_blocks': 'Model block byte-identical to v6; imaging retains v6 horizontal layout with taller traces and translated summaries',
        'model_offset_points': model_offset,
        'experimental_offset_points': experimental_offset,
        'source_model_page_points': [mw, mh],
        'source_experimental_page_points': [ew, eh],
        'half_width_points': half_width,
        'central_gutter_points': gutter,
        'uniform_scale': scale,
        'scaled_marker_diameter_pt': 3.3*scale,
        'scaled_marker_outline_pt': .2*scale,
        'scaled_axis_font_pt': 7.3*scale,
        'transform': 'Identical uniform scale for both blocks plus translation; right trace height and vertical spacing adjusted before assembly',
    }
    (OUT/'assembly_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(DEST)


if __name__ == '__main__':
    main()
