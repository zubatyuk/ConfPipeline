import numpy as np

def politzer(area, potential):
    S = area
    V = potential

    SV = S * V

    where_neg = V < 0
    where_pos = V >= 0

    V_min = np.min(V)
    V_max = np.max(V)
    V_avg = np.mean(V)
    V_avg_neg = np.mean(V[where_neg])
    V_avg_pos = np.mean(V[where_pos])

    S_tot = np.sum(S)
    S_neg = np.sum(S[where_neg])
    S_pos = np.sum(S[where_pos])

    sigma2_neg = np.sum((S[where_neg] * (V[where_neg] - V_avg_neg)) ** 2) / S_neg
    sigma2_pos = np.sum((S[where_pos] * (V[where_pos] - V_avg_pos)) ** 2) / S_pos
    sigma2_tot = sigma2_neg + sigma2_pos
    mui = sigma2_neg * sigma2_pos / sigma2_tot ** 2

    pi = np.sum(S * np.abs(V - V_avg)) / S_tot

    mui_sigma2 = mui * sigma2_tot

    return dict((name, eval(name)) for name in ['V_min', 'V_max', 'V_avg', 'V_avg_neg', 'V_avg_pos',
                                                'S_tot', 'S_neg', 'S_pos',
                                                'sigma2_neg', 'sigma2_pos', 'sigma2_tot', 'mui', 'pi', 'mui_sigma2'])



