"""Pair the fixed-matrix v0.3 evaluation with the indexed-222 unweighted evaluation."""
import hashlib
import json
from pathlib import Path
import sys

import fitz


def build(left, right, dest):
    dest.mkdir(exist_ok=True, parents=True)
    source = json.loads((left / 'evaluation.json').read_text())
    current = json.loads((right / 'unweighted/fit_report.json').read_text())
    left_n, right_n = source['evaluated_events'], current['fit_events']
    for path, n in [(left,left_n), (right,right_n)]:
        selection = json.loads((path/'selection.json').read_text())
        assert selection['native_bin_cap'] is None and selection['fit_events'] == n
        assert all(v['fitted_each_case'] * 3 == n for v in selection['foil_counts'].values())
    sides = [left, right]
    hashes = [hashlib.sha256((p/'unweighted/matrices/nps_shms_newfit_25521_unweighted.dat').read_bytes()).hexdigest() for p in sides]
    assert hashes == [source['matrix_sha256'], json.loads((right/'validation.json').read_text())['cases']['unweighted']['matrix_sha256']]
    plots = Path('unweighted/plots')
    new = plots/'new'
    delta = new/'25521_new_delta_vs_ytar.pdf'
    summary = plots/'fit_opt_matrix_25521_unweighted_new_matrix_diff.pdf'
    foils = ['foil_0_z_plus10', 'foil_1_z_0', 'foil_2_z_minus10']
    full = [(delta, 'Reconstructed Ytar versus momentum deviation', None),
            (summary, 'Target residuals | Gaussian fits', None)]
    for pdf in sorted((right/new).glob('*.pdf')):
        if pdf.name != delta.name:
            full.append((new/pdf.name, pdf.stem.removeprefix('25521_new_').replace('_', ' '), None))
    for pdf in sorted((right/new/'nps_legacy').glob('25521*.pdf')):
        full.append((new/'nps_legacy'/pdf.name, 'NPS | '+pdf.stem.removeprefix('25521_new_').replace('_', ' '), None))
    quick = full[:2].copy()
    intervals = ['(-15, -6)', '[-6, -2)', '[-2, 2)', '[2, 6)', '[6, 10)', '[10, 22)']
    for i, foil in enumerate(foils):
        files = sorted((right/new/foil).glob('delta_*.pdf'))
        assert len(files) == 6
        for j, pdf in enumerate(files):
            title = f'Foil {i} | delta {intervals[j]}%'
            full.append((new/foil/pdf.name, title, None))
            if j == 4:
                quick.extend([(new/foil/pdf.name, title+' | target residuals', [0]),
                              (new/foil/pdf.name, title+' | sieve reprojection', [4]),
                              (new/'nps_legacy/25521_new_xpdiff.pdf', title+' | Xptar by hole', [6*i+j])])
    labels = ['v0.3 clusters  /  original repository unweighted new matrix',
              '222 indexed clusters  /  unweighted SVD fitted to the updated data']
    plots_used = []
    counts = {}
    for name, sections in [('unweighted_complete_atlas', full), ('unweighted_highlights', quick)]:
        book = fitz.open()
        toc = []
        for rel, title, indices in sections:
            with fitz.open(left/rel) as a, fitz.open(right/rel) as b:
                assert len(a) == len(b) and len(a) > 0, rel
                pages = range(len(a)) if indices is None else indices
                toc.append([1, title, len(book)+1])
                if name == 'unweighted_complete_atlas':
                    plots_used.append(dict(source=rel.as_posix(), pages=len(a), first_atlas_page=len(book)+1))
                for k in pages:
                    page = book.new_page(width=1440, height=620)
                    page.insert_text((18,22), f'Run 25521   |   {title}   |   {k+1}/{len(a)}', fontsize=14)
                    for x, doc, label, pop, color in [(18,a,labels[0],f'{left_n:,} evaluated  /  {left_n//3:,} per foil',(.13,.28,.46)),
                                                     (728,b,labels[1],f'{right_n:,} evaluated  /  {right_n//3:,} per foil',(.05,.38,.32))]:
                        page.insert_text((x,46), label, fontsize=15, color=color)
                        page.insert_text((x,64), pop, fontsize=10, color=color)
                        # Flatten the rotated, offset ROOT CropBox before placing the vector page.
                        with fitz.open('pdf', doc.convert_to_pdf(from_page=k, to_page=k)) as visible:
                            page.show_pdf_page(fitz.Rect(x,73,x+692,596), visible, 0)
                    page.insert_text((18,612), 'Each panel uses its own foil/hole targets and selected events. Native axis ranges and color scales are retained.', fontsize=9)
        book.set_toc(toc)
        book.set_metadata({'title':'Run 25521 | v0.3 and indexed-222 unweighted matrix comparison'})
        file = dest/(name+'.pdf')
        book.save(file, garbage=4, deflate=True)
        book.close()
        with fitz.open(file) as check:
            assert len(check.get_toc()) == len(sections)
            for page in check:
                assert len(page.get_xobjects()) >= 2 and all(s in page.get_text() for s in labels)
            counts[name] = len(check)
    report = dict(left=source, right=dict(clusters=222, matrix_sha256=hashes[1], evaluated_events=right_n,
                  events_per_foil=right_n//3, rmse_xptar_mrad_ytar_mm_yptar_mrad=current['new_at_same_geometric_xtar']['rmse_mrad_mm_mrad']),
                  atlas_pages=counts, sections=plots_used)
    (dest/'atlas_sources.json').write_text(json.dumps(report, indent=2)+'\n')
    print(counts)


if __name__ == '__main__':
    build(*(Path(arg).resolve() for arg in sys.argv[1:]))
