"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Gauss con il Pivoting
"""

import numpy as np


def risolvi_tridiagonale(e, f, g, r):
    """
    e : array (diagonale inferiore/sub-diagonale)
    f : array (diagonale principale)
    g : array (diagonale superiore/sovra-diagonale)
    r : array (termini noti)
    """
    n = len(f)

    # Creiamo copie per non modificare i dati originali
    ee = e.copy()
    ff = f.copy()
    rr = r.copy()
    x = np.zeros(n)

    # Decomposizione
    for k in range(1, n):
        # ek = ek / f_{k-1}
        factor = ee[k] / ff[k - 1]
        ee[k] = factor  # Memorizziamo il moltiplicatore

        # fk = fk - ek * g_{k-1}
        ff[k] = ff[k] - factor * g[k - 1]

    # Sostituzione in avanti
    for k in range(1, n):
        # rk = rk - ek * r_{k-1}
        rr[k] = rr[k] - ee[k] * rr[k - 1]

    # Sostituzione all'indietro (Back substitution)
    # xn = rn / fn
    x[n - 1] = rr[n - 1] / ff[n - 1]

    for k in range(n - 2, -1, -1):
        # xk = (rk - gk * x_{k+1}) / fk
        x[k] = (rr[k] - g[k] * x[k + 1]) / ff[k]

    return x