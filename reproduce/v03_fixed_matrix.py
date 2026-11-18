"""Make the left atlas panels. Run with native ROOT in PATH and PyROOT installed.

Arguments: upstream checkout, v0.3 data_1.npz, empty working directory.
For a sparse checkout, argument four can supply the unmodified run25521 tree.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import uproot


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


repo, release, work = [Path(value).resolve() for value in sys.argv[1:4]]
historical = Path(sys.argv[4]).resolve() if len(sys.argv) > 4 else repo / 'run25521'
work.mkdir(exist_ok=True, parents=True)
out = work / 'run25521'
out.mkdir(exist_ok=False)
for name in ['config', 'diagnostics', 'spectrometer_root.h', 'spectrometer_config.h', 'spectrometer_profiles.def']:
    src = repo / name
    if src.is_dir():
        shutil.copytree(src, work / name)
    else:
        shutil.copy2(src, work / name)
scripts = out / 'reproduce'
shutil.copytree(historical / 'reproduce', scripts)
shutil.copyfile(repo / 'run25521/reproduce/prepare.py', scripts / 'prepare.py')
meta = out / 'config/SHMS_run25521'
meta.mkdir(parents=True)
matrix = out / 'unweighted/matrices/nps_shms_newfit_25521_unweighted.dat'
matrix.parent.mkdir(parents=True)
shutil.copyfile(historical / matrix.relative_to(out), matrix)
assert digest(matrix) == '5e8ba5d81ba6ae75d8de5e0e50047822e83cd2a186d2dc83201a0df2e07d577c'
assert digest(release) == 'de16b6577a8ef2b08d0659f91a455267c624de7486043eb928f9d532c2c6d9ff'
keys = ('run source_full_entry cluster_id foil_id hole_col hole_row fp_x_cm fp_y_cm fp_xp_rad '
        'fp_yp_rad delta_percent beam_x_cm beam_y_cm reco_xtar_cm reco_xptar_rad reco_yptar_rad reco_ytar_cm').split()
with np.load(release, allow_pickle=False) as src:
    inputs = {key: src[key] for key in keys}
    inputs['z_reco_cm'] = src['z_tar_reco_cm']
    inputs['index_review_required'] = np.zeros(len(inputs['run']), dtype=bool)
assert len(np.unique(inputs['cluster_id'])) == 219
np.savez(work / 'data_1.npz', **inputs)
subprocess.run([sys.executable, scripts / 'prepare.py', work, work], check=True, cwd=work)
selection = json.loads((out / 'selection.json').read_text())
assert selection['native_bin_cap'] is None
with uproot.open(work / 'SHMS_run25521/balanced/Optics_25521_-1_fit_tree_gmm.root') as file:
    events = file['TFit'].arrays(library='np')
sys.path.insert(0, str(work / 'diagnostics'))
from make_unweighted_plots import evaluate

pred = evaluate(matrix, events)
targets = np.column_stack([events['xptarT'], events['ytarT']/100, events['yptarT']])
residuals = (pred-targets)*1000
n = selection['fit_events']
assert len(pred) == n and np.isfinite(residuals).all()
table = {k: events[k].astype(float) for k in
         'ys xs ysT xsT xtar ztar xptarT yptarT ytarT ztarT delta xpfp ypfp xfp yfp'.split()}
table.update(xptar=pred[:, 0], ytar=pred[:, 1]*100, yptar=pred[:, 2])
tree = work / 'SHMS_run25521/unweighted_new.root'
with uproot.recreate(tree) as file:
    file.mktree('TFit', table)
with uproot.open(tree) as file:
    assert file['TFit'].num_entries == n
dest = out / 'unweighted/plots/new'
dest.mkdir(parents=True)
for axis in ['xptar', 'ytar', 'yptar']:
    macro = work / 'diagnostics/fit' / f'plot_{axis}_residuals.C'
    cmd = f'{macro}(25521,"{tree}","{meta}/optics_metadata.dat","{dest}","25521_new")'
    result = subprocess.run(['root', '-l', '-b', '-q', cmd], cwd=work, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, check=True)
    assert f'start loop {n}' in result.stdout and 'Error in <' not in result.stdout, result.stdout[-3000:]

# Limit the pinned plotting programs to this matrix and derive count checks from v0.3.
for name in ['foil_slice_diagnostics.py', 'run_nps_legacy.py', 'plot_delta_vs_ytar.py']:
    path = scripts / name
    code = path.read_text().replace("['unweighted','weighted']", "['unweighted']").replace("['old','new']", "['new']")
    code = code.replace(",('weighted','Supplied weights')", "").replace('Unweighted legacy fit', 'Published unweighted matrix')
    code = code.replace('73188', str(n)).replace('24396', str(n//3))
    if name == 'plot_delta_vs_ytar.py':
        code = '\n'.join(line for line in code.splitlines() if not line.startswith("assert np.array_equal(arrays['unweighted']"))
    if name == 'run_nps_legacy.py':
        code = code.split("lines=['#")[0]
    path.write_text(code)
    subprocess.run([sys.executable, path, work], cwd=work, check=True)

# Same three histograms, bin limits, Gaussian fits and pad order as fit_opt_matrix_gmm.C.
import ROOT
ROOT.gROOT.SetBatch(True)
canvas = ROOT.TCanvas('cnewdiff', 'Newfit diff target', 800, 800)
canvas.Divide(2, 2)
histograms = []
for pad, (name, title, values, bins, low, high) in enumerate([
    ('hytarnewdiff', 'ytar new diff(cm)', residuals[:, 1]/10, 70, -5, 5),
    ('hyptarnewdiff', 'yptar new diff (mr)', residuals[:, 2], 100, -10, 10),
    ('hxptarnewdiff', 'xptar new diff (mr)', residuals[:, 0], 100, -15, 15),
], start=1):
    canvas.cd(pad)
    hist = ROOT.TH1F(name, title, bins, low, high)
    hist.FillN(n, np.ascontiguousarray(values), np.ones(n))
    assert hist.GetEntries() == n
    hist.Draw()
    hist.Fit('gaus', 'Q')
    histograms.append(hist)
canvas.Print(str(dest.parent / 'fit_opt_matrix_25521_unweighted_new_matrix_diff.pdf'))
audit = dict(source_sha256=digest(release), matrix_sha256=digest(matrix), clusters=219,
             available_events=selection['events'], evaluated_events=n, events_per_foil=n//3,
             rmse_xptar_mrad_ytar_mm_yptar_mrad=np.sqrt(np.mean(residuals**2,axis=0)).tolist())
(out / 'evaluation.json').write_text(json.dumps(audit, indent=2) + '\n')
