from pathlib import Path
import sys
import numpy as np
import uproot

source, out = map(Path, sys.argv[1:])
out.mkdir(parents=True, exist_ok=True)
keys = """run source_full_entry foil_id hole_col hole_row cluster_id fp_x_cm fp_y_cm
fp_xp_rad fp_yp_rad delta_percent beam_x_cm beam_y_cm reco_xtar_cm reco_xptar_rad
reco_yptar_rad reco_ytar_cm index_review_required""".split()
with np.load(source, allow_pickle=False) as src:
    data = {key: src[key] for key in keys}
    data["z_reco_cm"] = src["z_tar_reco_cm"]
assert len(data["run"]) == 106058 and data["index_review_required"].sum() == 147
np.savez_compressed(out / "data_1.npz", **data)
with uproot.recreate(out / "data_1.root") as file:
    file.mktree("train", {k: v.dtype for k, v in data.items()})
    file["train"].extend(data)
