from .functions import DECONVO, pre_marker_select
import pandas as pd
import collections
import numpy as np
import os


def ARIC(
    mix_path,
    ref_path,
    save_path=None,
    marker_path=None,
    selected_marker=False,
    scale=0.1,
    delcol_factor=10,
    iter_num=10,
    confidence=0.75,
    w_thresh=10,
    unknown=False,
    is_methylation=False,
):
    """
     This function is used for deconvolute DNA methylation or gene expression data.

     ARIC(mix_path, ref_path, save_path=None, marker_path=None, selected_marker=False, scale=0.1, delcol_factor=10,
         iter_num=10, confidence=0.75, w_thresh=10, unknown=False, is_methylation=False)
    {P}arameters:
        mix_path: Path to mixture data, must be an csv file with colnames and rownames.
        ref_path: Path to reference data, must be an csv file with colnames and rownames.
        save_path: Where to save the deconvolution results. Default: mix_path_prefix_prop.csv.
        marker_path: Path to the user specificed markers. Must be an csv file.
        selected_marker: Output selected marker for every sample. Marker files will be saved in a folder named "sample_marker.csv".
        scale: Used for controlling the convergence of SVR. A smaller value makes the convergence much faster. Default: 0.1.
        delcol_factor: Used for controlling the extent of removing collinearity. Default: 10.
        iter_num: Iterative numbers of outliers detection. Default: 10.
        confidence: Ratio of remained markers in each outlier detection loop. Default: 0.75.
        w_thresh: Threshold to cut the weights designer. Default: 10.
        unknown: Whether to estimate unknown content proportion.
        is_methylation: Whether the data type belongs to methylation data. If true, preliminary marker selection will be performed.
    """

    print("---------------------------------------------")
    print("--------------WELCOME TO ARIC----------------")
    print("---------------------------------------------")

    ref = pd.read_csv(ref_path, index_col=0)
    mix = pd.read_csv(mix_path, index_col=0)

    if save_path is None:
        save_path = os.path.splitext(mix_path)[0] + "_prop.csv"

    if selected_marker:
        tmppath = os.path.splitext(mix_path)[0]
        os.makedirs(tmppath, exist_ok=True)
    else:
        tmppath = None

    raw_markers = list(ref.index)
    cell_type = ref.columns.values
    samples = mix.columns.values
    prop = collections.OrderedDict()
    prop["cell types"] = cell_type
    
    if marker_path is not None:
        reference = []
        mixture = []
        markers = pd.read_csv(marker_path, index_col=0)
        markers = markers.index.values

        # valid marker
        markers = list(set(markers).intersection(set(mix.index)))
        markers = list(set(markers).intersection(set(ref.index)))

        for i in range(len(markers)):
            reference.append(ref.loc[markers[i]])
            mixture.append(mix.loc[markers[i]])
    else:
        reference, mixture = ref, mix

    reference = np.asarray(reference)
    mixture = np.asarray(mixture)

    if is_methylation:
        reference, mixture = pre_marker_select(reference, mixture)

    print("Data reading finished!")

    print("ARIC Engines Start, Please Wait......")

    prop_predict = DECONVO(
        scale * reference,
        scale * mixture,
        raw_markers=raw_markers,
        raw_ref=np.asarray(ref),
        delcol_factor=delcol_factor,
        iter_num=iter_num,
        confidence=confidence,
        w_thresh=w_thresh,
        unknown=unknown,
        celltypes=cell_type,
        samples=samples,
        tmppath=tmppath,
    )

    print("Deconvo Results Saving!")

    for i in range(len(samples)):
        prop[samples[i]] = []
        for j in range(len(cell_type)):
            prop[samples[i]].append(prop_predict[j, i])
    prop = pd.DataFrame(prop)

    prop.to_csv(save_path, index=False)
    print("Finished!")
