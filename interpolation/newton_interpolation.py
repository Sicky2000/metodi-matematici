"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Interpolazione di Newton
"""

import numpy as np


def interpolazione_newton(x, y, xi_target):
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

    print("-" * 50)
    print(f"Stima all'Ordine 0: {yint:.12f}")

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

# --- MAIN DI ESEMPIO ---
if __name__ == "__main__":
    # 1. DEFINIZIONE DELLA FUNZIONE
    def funzione(x):
        return np.log(x)

    # 2. DEFINIAMO SOLO I NODI X (I PUNTI DI CAMPIONAMENTO)
    x_dati = [1, 4, 6]

    # 3. CALCOLIAMO LE Y DINAMICAMENTE CHIAMANDO LA FUNZIONE
    #    Questo garantisce la massima precisione float
    y_dati = [funzione(x) for x in x_dati]

    # Punto in cui vogliamo interpolare
    target = 2

    print(f"Calcolo interpolazione per x = {target}")
    print(f"Nodi X utilizzati: {x_dati}")
    print(f"Valori Y calcolati: {[round(y, 6) for y in y_dati]}")  # Stampa arrotondata per leggibilità
    print("-" * 50)

    risultato = interpolazione_newton(x_dati, y_dati, target)

    # Calcolo dell'errore vero
    valore_vero = funzione(target)
    errore_assoluto = abs(risultato - valore_vero)

    print(f"Risultato Interpolato: {risultato:.5f}")
    print(f"Valore Vero (ln({target})): {valore_vero:.5f}")
    print(f"Errore Assoluto: {errore_assoluto:.5e}")