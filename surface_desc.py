import numpy as np


def politzer(area, potential, suffix=None):
    S = area
    V = potential

    SV = S * V
    where_neg = V < 0
    where_pos = V >= 0

    V_min = np.min(V)
    V_max = np.max(V)

    S_tot = np.sum(S)
    S_neg = np.sum(S[where_neg])
    S_pos = np.sum(S[where_pos])

    V_avg = np.sum(SV) / S_tot

    if S_neg != 0:
        V_avg_neg = np.sum(SV[where_neg]) / S_neg
        sigma2_neg = np.sum((S[where_neg] * (V[where_neg] - V_avg_neg)) ** 2) / S_neg
    else:
        V_avg_neg = 0.0
        sigma2_neg = 0.0

    if S_pos != 0:
        V_avg_pos = np.sum(SV[where_pos]) / S_pos
        sigma2_pos = np.sum((S[where_pos] * (V[where_pos] - V_avg_pos)) ** 2) / S_pos
    else:
        V_avg_pos = 0.0
        sigma2_pos = 0.0

    sigma2_tot = sigma2_neg + sigma2_pos
    mui = sigma2_neg * sigma2_pos / sigma2_tot ** 2
    pi = np.sum(S * np.abs(V - V_avg)) / S_tot

    d = {
        'S_neg_ratio': S_neg / S_tot, 'S_pos_ratio': S_pos / S_tot,
        'V_min': V_min, 'V_max': V_max, 'V_avg': V_avg, 'V_avg_neg': V_avg_neg, 'V_avg_pos': V_avg_pos,
        'sigma2_tot': sigma2_tot, 'sigma2_neg': sigma2_neg, 'sigma2_pos': sigma2_pos,
        'mui': mui, 'pi': pi, 'mui_sigma2': mui * sigma2_tot
    }

    if suffix:
        for k in [x for x in d]:
            d[k + suffix] = d.pop(k)
    return d
