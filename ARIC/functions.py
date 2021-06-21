import numpy as np
from tqdm import tqdm
from sklearn.svm import NuSVR
import statsmodels.api as sm
import pandas as pd
import collections


def DECONVO(
    ref,
    mix,
    raw_markers,
    raw_ref,
    delcol_factor=10,
    iter_num=10,
    confidence=0.75,
    w_thresh=10,
    unknown=False,
    celltypes=None,
    samples=None,
    tmppath=None,
):
    """
    :param ref: reference data
    :param mix: mixture data
    :param delcol_factor: control the extent of removing collinearity
    :param iter_num: iterative numbers of outliers detection
    :param confidence: ratio of remained markers in each outlier detection loop
    :param w_thresh: threshold to cut the weights designer
    :param unknown: if there is unknown content
    :return: proportion data
    """
    reference, mixtureData = filt_zeros(ref, mix)
    markerNum = np.size(reference, 0)
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(mixtureData, 1)
    bestReference = []
    bestMarker = []
    conditionNumberHistory = []
    bestNumberHistory = []
    proportionDeconvolution = (
        np.zeros([tissueNum, numOfSamples]) if unknown == False else np.zeros([tissueNum + 1, numOfSamples])
    )
    for i in tqdm(range(numOfSamples)):

        # get sample name
        sample_name = samples[i]

        selectArea = np.arange(markerNum)
        selectMixture = mixtureData[selectArea, i]
        selectReference = reference[selectArea, :]
        minimumConditionNumber = 10 ** (100)
        endNumber = np.size(selectReference, 0)
        for selectedNumber in range(int(endNumber / delcol_factor)):
            minDistances = 10 ** (50)
            for j in range(tissueNum):
                for k in range(j + 1, tissueNum):
                    distances = selectReference[:, j] - selectReference[:, k]
                    distances = np.sqrt(np.sum(np.multiply(distances, distances)))
                    if distances < minDistances:
                        minDistances = distances
                        closetJ = j
                        closetK = k
            sumData = selectReference[:, closetJ] + selectReference[:, closetK]
            area = sumData == 0
            sumData[area] = 10 ** (-100)
            collinearity = np.abs(selectReference[:, closetJ] - selectReference[:, closetK]) / (sumData)
            collinearityIndex = np.argsort(collinearity)
            area = np.ones(np.size(selectReference, 0))
            area = area.astype(np.bool)
            area[collinearityIndex[0]] = False
            selectArea = np.arange(np.size(selectReference, 0))
            selectArea = selectArea[area]
            selectMixture = selectMixture[selectArea]
            selectReference = selectReference[selectArea, :]
            ConditionNumber = cmpConditionNumber(selectReference, selectMixture)
            conditionNumberHistory.append(ConditionNumber)
            if ConditionNumber < minimumConditionNumber:
                minimumConditionNumber = ConditionNumber
                bestReference = selectReference
                bestMarker = np.zeros([np.size(selectReference, 0), 1])
                bestMarker[:, 0] = selectMixture
                bestNumber = selectedNumber

        t = RobustSVR(
            bestReference,
            bestMarker,
            raw_markers=raw_markers,
            raw_ref=raw_ref,
            iter_num=iter_num,
            confidence=confidence,
            w_thresh=w_thresh,
            unknown=unknown,
            sample_name=sample_name,
            celltypes=celltypes,
            tmppath=tmppath,
        )
        bestNumberHistory.append(bestNumber)
        proportionDeconvolution[:, i] = t[:, 0]
    return proportionDeconvolution


def cmpConditionNumber(referenceSelect, mixtureSelect):
    """
    compute the componentwise condition number for each reference and mixture data
    :param referenceSelect:
    :param markerSelect:
    :return:
    """
    pinvReferenceSelect = np.linalg.pinv(referenceSelect)
    bNorm = np.linalg.norm(mixtureSelect)
    maxConditionNumber = 0
    tissueNumber = np.size(referenceSelect, 1)
    for i in range(tissueNumber):
        tq = pinvReferenceSelect[i, :]
        conditionNumber = (bNorm / np.abs(np.dot(tq, mixtureSelect))) * np.linalg.norm(tq)
        if conditionNumber > maxConditionNumber:
            maxConditionNumber = conditionNumber
    return maxConditionNumber


def stdRes(res, d, H):
    """
    compute the residuals
    :param res:
    :param d:
    :param H:
    :return:
    """
    res_std = np.zeros([np.size(res, 0), np.size(res, 1)])
    s = np.sqrt(np.sum(np.power(res, 2)) / d)
    for i in range(np.size(res, 0)):
        res_std[i, 0] = res[i, 0] / (s * (1 - H[i, i]))
    return res_std


def RobustSVR(
    reference,
    mixtureData,
    raw_markers,
    raw_ref,
    iter_num=10,
    confidence=0.75,
    w_thresh=10,
    unknown=False,
    sample_name=None,
    celltypes=None,
    tmppath=None,
):
    """
    outlier detection, weights designer and SVR deconvolution
    :param reference:
    :param mixtureData:
    :param iter_num:
    :param confidence:
    :param w_thresh:
    :param unknown:
    :return:
    """
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(mixtureData, 1)
    markerNumber = np.size(reference, 0)
    proportionDeconvolution = (
        np.zeros([tissueNum, numOfSamples]) if unknown == False else np.zeros([tissueNum + 1, numOfSamples])
    )
    for i in range(numOfSamples):

        iterReference = reference
        itermixtureData = mixtureData[:, i]
        mixture = sm.RLM(itermixtureData, iterReference).fit()
        test = mixture.params
        t = test / np.sum(test) if unknown == False else test
        c1 = np.zeros([np.size(iterReference, 0), 1])
        c1[:, 0] = itermixtureData[:]
        t1 = np.zeros([tissueNum, 1])
        t1[:, 0] = t
        c2 = np.dot(iterReference, t1)
        res = c1 - c2
        s_2 = np.sum(np.power(res, 2)) / (np.size(iterReference, 0) - np.size(iterReference, 1))
        res_std = np.abs(res / np.sqrt(s_2))
        res_Sort = np.sort(res_std[:, 0])
        T = res_Sort[int(confidence * np.size(res_Sort))]
        memRef = np.zeros([np.size(iterReference, 0), np.size(iterReference, 1)])
        memRef[:, :] = iterReference[:, :]
        memMix = np.zeros(np.size(itermixtureData))
        memMix[:] = itermixtureData[:]
        for j in range(iter_num):
            mixture = sm.RLM(itermixtureData, iterReference).fit()
            test = mixture.params
            t = test / np.sum(test) if unknown == False else test
            c1 = np.zeros([np.size(iterReference, 0), 1])
            c1[:, 0] = itermixtureData[:]
            t1 = np.zeros([tissueNum, 1])
            t1[:, 0] = t
            c2 = np.dot(iterReference, t1)
            res = c1 - c2
            s_2 = np.sum(np.power(res, 2)) / (np.size(iterReference, 0) - np.size(iterReference, 1))
            res_std = res / np.sqrt(s_2)
            iterSelected = np.arange(np.size(iterReference, 0))
            area = np.abs(res_std[:, 0]) <= T
            iterSelected = iterSelected[area]
            iterReference = iterReference[iterSelected, :]
            itermixtureData = itermixtureData[iterSelected]
            if np.size(iterReference, 0) < int(tissueNum):
                iterReference = memRef
                itermixtureData = memMix
                break
            if np.size(iterReference, 0) < int(0.5 * markerNumber):
                break
            memRef = np.zeros([np.size(iterReference, 0), np.size(iterReference, 1)])
            memRef[:, :] = iterReference[:, :]
            memMix = np.zeros(np.size(itermixtureData))
            memMix[:] = itermixtureData[:]

        # get selected markers
        if tmppath is not None:
            marker_name = ["marker" + str(i) for i in range(iterReference.shape[0])]

            # output single sample deconvolution information
            # mix_new = pd.DataFrame(data=itermixtureData, index=marker_name, columns=[sample_name])
            # ref_new = pd.DataFrame(data=iterReference, index=marker_name, columns=celltypes)
            # mix_new.to_csv(tmppath + sample_name + "_mix.csv", index=True)
            # ref_new.to_csv(tmppath + sample_name + "_ref.csv", index=True)

            cos_M = compute_similarity(iterReference, raw_ref)
            markers_output = collections.OrderedDict()
            markers_output["markers"], markers_output["confidence"] = [], []
            for c in range(np.size(cos_M, 0)):
                idx = np.argsort(-cos_M[c, :])
                sim_c = cos_M[c, idx[0]]
                markers_output["markers"].append(raw_markers[idx[0]])
                markers_output["confidence"].append(sim_c)
            markers_output = pd.DataFrame(markers_output)
            markers_output.to_csv(tmppath + "/" + sample_name + "_markers.csv", index=True)

        weights = weightsDesigner(iterReference, itermixtureData, w_thresh=w_thresh)

        t = nuSVR(iterReference, itermixtureData.reshape([-1, 1]), weights, unknown=unknown)
        t = t[:, 0]
        if unknown == False:
            proportionDeconvolution[:, i] = t
        else:
            proportionDeconvolution[0:-1, i] = t
            proportionDeconvolution[-1, i] = max(1 - np.sum(t), 0)
    return proportionDeconvolution


def compute_similarity(iter_ref, raw_ref):
    m_dot = np.dot(iter_ref, np.transpose(raw_ref))
    iter_L = np.sqrt(np.sum(np.power(iter_ref, 2), axis=1))
    raw_L = np.sqrt(np.sum(np.power(raw_ref, 2), axis=1))
    cos_M = m_dot / np.dot(iter_L.reshape((-1, 1)), np.transpose(raw_L.reshape(-1, 1)))
    return cos_M


def weightsDesigner(ref, mix, w_thresh=10):
    mixture = sm.RLM(mix, ref).fit()
    test = mixture.params
    x_pre = test / np.sum(test)
    weights = np.abs(np.dot(ref, x_pre))
    for i in range(np.size(weights)):
        weights[i] = 1 / weights[i]
    weights = weights / np.min(weights)
    for i in range(np.size(weights)):
        if weights[i] > w_thresh:
            weights[i] = w_thresh
    return weights / np.mean(weights)


def nuSVR(reference, mixtureData, weights, unknown=False):
    nu = [0.25, 0.50, 0.75]
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(mixtureData, 1)
    proportionDeconvolution = np.zeros([tissueNum, numOfSamples])
    nuHistory = []
    p0 = np.zeros([3, tissueNum, numOfSamples])
    for i in range(0, 3, 1):
        nuI = nu[i]
        clf = NuSVR(nu=nuI, kernel="linear")
        for j in range(numOfSamples):
            clf.fit(reference, mixtureData[:, j], sample_weight=weights)
            t = clf.coef_
            t1 = np.zeros(tissueNum)
            t1[:] = t[0, :]
            area = t1 < 0
            t1[area] = 0
            t1 = t1 / np.sum(t1) if unknown == False else t1
            p0[i, :, j] = t1[:]
    for i in range(numOfSamples):
        minRMSE = 10 ** (50)
        truth = np.zeros([np.size(reference, 0), 1])
        truth[:, 0] = mixtureData[:, i]
        for k in range(0, 3, 1):
            pVector = np.zeros([tissueNum, 1])
            pVector[:, 0] = p0[k, :, i]
            temp = np.dot(reference, pVector)
            error = truth - temp
            error = np.sqrt(np.mean(np.multiply(error, error)))
            if error < minRMSE:
                minRMSE = error
                proportionDeconvolution[:, i] = p0[k, :, i]
                bestNu = k
        nuHistory.append(nu[bestNu])
    return proportionDeconvolution


def pre_marker_select(reference, mixtureData):
    cellTypeNumber = np.size(reference, 1)
    markerNumber = np.size(reference, 0)
    selectedCpG = np.arange(markerNumber)
    area = np.zeros(markerNumber)
    area = area.astype(np.bool)
    for i in range(cellTypeNumber):
        for j in range(i + 1, cellTypeNumber):
            temp = reference[:, [i, j]]
            tempSum = np.sum(temp, axis=1)
            tempSum1 = np.zeros([markerNumber, 2])
            tempSum1[:, 0] = tempSum
            tempSum1[:, 1] = tempSum
            temp = temp / tempSum1
            pairSortIndexIncrease = np.argsort(temp[:, 0])
            pairSortIndexDecrease = np.argsort(temp[:, 1])
            area[pairSortIndexIncrease[0:100]] = True
            area[pairSortIndexDecrease[0:100]] = True
    selectedCpG = selectedCpG[area]
    ref = reference[selectedCpG, :]
    mix = mixtureData[selectedCpG, :]
    return ref, mix


def filt_zeros(ref, mix):
    ref1 = []
    mix1 = []
    for i in range(np.size(ref, 0)):
        if np.max(ref[i, :]) > 0:
            ref1.append(ref[i, :])
            mix1.append(mix[i, :])
    return np.asarray(ref1), np.asarray(mix1)
