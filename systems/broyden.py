"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Broyden
"""

import numpy as np


def broyden(func, x0, tol, max_iter, b0=None):
    """
        func (callable): funzione che accetta x e restituisce un array numpy.
        x0 (array_like): vettore tentativo iniziale (initial guess).
        b0 (array_like, opzionale): approssimazione iniziale dello Jacobiano.
                                    se None, usa la matrice Identità.
        tol (float): tolleranza per la convergenza (norma del passo).
        max_iter (int): numero massimo di iterazioni.
    """

    # configurazione e conversione dei tipi
    # assicura che x0 sia un array numpy 1D di float
    x0 = np.array(x0, dtype=float).flatten()
    n = len(x0)

    # inizializza b0 (Approssimazione dello Jacobiano)
    if b0 is None:
        b0 = np.eye(n) # matrice Identità se non fornita
    else:
        b0 = np.array(b0, dtype=float)

    # Valuta la funzione al punto iniziale
    fx0 = np.array(func(x0)).flatten()

    difv = [] # inizializza la lista per lo storico della convergenza
    diff = tol + 1.0 # imposta un valore iniziale > tol per garantire l'ingresso nel ciclo while
    i = 0 # contatore delle iterazioni

    # inizializzazione preventiva per evitare errori se il loop non dovesse partire

    x1 = x0
    fx1 = fx0

    # ciclo Principale
    while diff >= tol and i < max_iter:
        i += 1

        # risolve il sistema lineare: B * deltax = -f(x)
        try:
            deltax = np.linalg.solve(b0, -fx0)
        except np.linalg.LinAlgError:
            print("Error: Jacobian approximation became singular.")
            break

        # calcola il nuovo punto x
        x1 = x0 + deltax
        fx1 = np.array(func(x1)).flatten()

        # aggiornamento di Broyden
        # formula: B_new = B_old + (df * deltax') / (deltax' * deltax)
        delta_norm_sq = np.dot(deltax, deltax)

        # evita la divisione per zero se il passo è troppo piccolo
        if delta_norm_sq > 1e-16:
            b0 = b0 + np.outer(fx1, deltax) / delta_norm_sq

        # aggiorna le metriche di convergenza
        diff = np.linalg.norm(deltax)
        difv.append(diff)

        # aggiorna le variabili per la prossima iterazione
        x0 = x1
        fx0 = fx1

    # finalizzazione dei risultati
    zero = x1
    res = np.linalg.norm(fx1)

    # controllo se ho raggiunto il numero massimo di iterazioni (FALLIMENTO)
    if i == max_iter and diff > tol:
        print("Broyden terminato: Numero massimo di iterazioni raggiunto senza convergenza.")

    return zero, res, i, np.array(difv)