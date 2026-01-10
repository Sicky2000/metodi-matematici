"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Thomas per i Sistemi Tridiagonali
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

# --- MAIN DI ESEMPIO ---
if __name__ == "__main__":
    # 1. DEFINIZIONE DEI DATI
    # Sistema:
    # [ 2 -1  0  0 ]
    # [-1  2 -1  0 ]
    # [ 0 -1  2 -1 ]
    # [ 0  0 -1  2 ]

    n = 4

    # Diagonale principale (f)
    f_in = np.array([2.0, 2.0, 2.0, 2.0])

    # Diagonale inferiore (e) - il primo elemento è ignorato (0.0)
    e_in = np.array([0.0, -1.0, -1.0, -1.0])

    # Diagonale superiore (g) - l'ultimo elemento è ignorato (0.0)
    g_in = np.array([-1.0, -1.0, -1.0, 0.0])

    # Termini noti (r)
    # Se la soluzione attesa è [1, 1, 1, 1], allora r deve essere [1, 0, 0, 1]
    r_in = np.array([1.0, 0.0, 0.0, 1.0])

    print(f"Dimensione sistema: {n}x{n}")
    print(f"Termini noti: {r_in}")

    # 2. CHIAMATA ALLA FUNZIONE
    soluzione = risolvi_tridiagonale(e_in, f_in, g_in, r_in)

    # 3. STAMPA DEI RISULTATI
    print("\nRisultato calcolato (x):")
    for i, val in enumerate(soluzione):
        print(f"x[{i + 1}] = {val:.4f}")