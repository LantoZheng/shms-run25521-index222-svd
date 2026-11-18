"""Compare both unweighted matrices on shared events under each target convention.

Arguments: patched upstream checkout, v0.3 work directory, indexed-222 work
directory, published v0.3-side matrix, fitted indexed-222 matrix, output CSV.
"""
import csv
from pathlib import Path
import sys

import numpy as np
import uproot

repo, left, right, original, refit, dest = map(Path, sys.argv[1:])
sys.path.insert(0, str(repo.resolve() / 'diagnostics'))
from make_unweighted_plots import evaluate

def read(work):
    path = work / 'SHMS_run25521/balanced/Optics_25521_-1_fit_tree_gmm.root'
    keys = ('source_full_entry foil_id hole_col hole_row xfp xpfp yfp ypfp '
            'xtar xptarT ytarT yptarT').split()
    with uproot.open(path) as file:
        data = file['TFit'].arrays(keys, library='np')
    assert np.all(np.diff(data['source_full_entry']) > 0)
    return data

a, b = read(left), read(right)
ids, ia, ib = np.intersect1d(a['source_full_entry'], b['source_full_entry'], return_indices=True)
a = {k:v[ia] for k,v in a.items()}
b = {k:v[ib] for k,v in b.items()}
assert len(ids) and all(np.array_equal(a[k], b[k]) for k in ('source_full_entry','xfp','xpfp','yfp','ypfp'))
same = (a['foil_id']==b['foil_id']) & (a['hole_col']==b['hole_col']) & (a['hole_row']==b['hole_row'])
assert same.any() and all(np.array_equal(a[k][same], b[k][same]) for k in ('xtar','xptarT','ytarT','yptarT'))
with dest.open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['target_index','matrix','population','events','xptar_mrad','ytar_mm','yptar_mrad'])
    for target, data in [('v0.3',a), ('indexed222',b)]:
        truth = np.column_stack((data['xptarT'],data['ytarT']/100,data['yptarT']))
        for label, matrix in [('upstream_unweighted',original), ('indexed222_unweighted',refit)]:
            residual = (evaluate(matrix,data)-truth)*1000
            for population,mask in [('all',np.ones(len(ids),bool)), ('unchanged_labels',same)]:
                writer.writerow([target,label,population,int(mask.sum()),*np.sqrt(np.mean(residual[mask]**2,axis=0))])
