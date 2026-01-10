"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Interpolazione di Newton con i nodi di Chebyshev
"""

import numpy as np

def newton_chebyshev(x, y, xi_target):
    """
    x: Array dei punti x (dati noti)
    y: Array dei punti y (dati noti)
    xi_target: Il valore x per cui vogliamo stimare y (punto di interpolazione)
    """
    n = len(x)

    # inizializza la tabella delle differenze divise (matrice n x n)
    # fdd sta per "Finite Divided Differences"
    fdd = np.zeros((n, n))

    # riempie la prima colonna con i valori y iniziali
    for i in range(n):
        fdd[i][0] = y[i]

    # calcola le Differenze Divise (La struttura piramidale)
    # j rappresenta la colonna (l'ordine della differenza)
    for j in range(1, n):
        for i in range(n - j):
            numeratore = fdd[i + 1][j - 1] - fdd[i][j - 1]
            denominatore = x[i + j] - x[i]
            fdd[i][j] = numeratore / denominatore

    # valutazione del polinomio (interpolazione) nel punto xi_target
    xterm = 1
    yint = fdd[0][0]  # Il primo termine è b0 (che è y[0])

    for ordine in range(1, n):
        # calcola il prodotto dei termini (x - x0)(x - x1)...
        xterm = xterm * (xi_target - x[ordine - 1])

        # calcola il valore del nuovo termine da aggiungere
        termine_aggiuntivo = fdd[0][ordine] * xterm

        # aggiorna la stima totale
        yint_precedente = yint
        yint += termine_aggiuntivo

        # calcolo dell'errore approssimato (ea)
        ea = yint - yint_precedente

        print(f"Stima all'Ordine {ordine}: {yint:.12f} (Errore stimato: {ea:.4e})")

    return yint

# Calcolo nodi di Chebyshev
def calcola_nodi_chebyshev(a, b, n):
    """
    a : float -> Estremo sinistro dell'intervallo
    b : float -> Estremo destro dell'intervallo
    n : int   -> Numero di nodi da generare
    """

    # Generazione degli indici k
    k = np.arange(n)

    # Calcolo dei nodi nell'intervallo standard [-1, 1]
    t_k = np.cos((2 * k + 1) / (2 * n) * np.pi)

    # Mappatura (Trasformazione Affine) nell'intervallo [a, b]
    x_k = (a + b) / 2 + (b - a) / 2 * t_k

    return x_k


# --- MAIN DI ESEMPIO ---
# 1. DEFINIAMO LA FUNZIONE DA INTERPOLARE
def funzione(x):
    return 1 / (1 + 25 * x ** 2)

if __name__ == "__main__":
    # 2. PARAMETRI
    a, b = -1, 1  # Intervallo
    n_nodi = 10  # Numero di nodi (grado polinomio = n-1)
    xi_target = 0.5  # Punto dove vogliamo stimare il valore (Test singolo)

    # 3. GENERAZIONE DEI DATI USANDO CHEBYSHEV
    x_cheb = calcola_nodi_chebyshev(a, b, n_nodi)
    y_cheb = funzione(x_cheb)

    print(f"Nodi di Chebyshev calcolati: \n{np.round(x_cheb, 4)}\n")

    # 4. ESECUZIONE INTERPOLAZIONE NEL PUNTO XI_TARGET
    print(f"--- Interpolazione nel punto x = {xi_target} ---")
    valore_stimato = newton_chebyshev(x_cheb, y_cheb, xi_target)
    valore_vero = funzione(xi_target)

    print("-" * 50)
    print(f"Risultato Interpolazione: {valore_stimato:.6f}")
    print(f"Valore Reale (Runge):     {valore_vero:.6f}")
    print(f"Errore Assoluto:          {abs(valore_stimato - valore_vero):.6e}")