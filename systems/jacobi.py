"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Iterazione di Jacobi
"""

import numpy as np

def jacobi(A, b, x0, tolleranza, max_iter):
    """
    A: Matrice dei coefficienti (nxn)
    b: Vettore dei termini noti (n)
    x0: Ipotesi iniziale (n)
    """
    n = len(b)
    x_old = x0.copy()
    x_new = np.zeros(n)

    for k in range(max_iter):
        # Iteriamo su ogni riga i
        for i in range(n):
            sigma = 0
            # Somma su j diverso da i
            for j in range(n):
                if j != i:
                    sigma += A[i][j] * x_old[j]

            # Formula di aggiornamento
            x_new[i] = (b[i] - sigma) / A[i][i]

        # Controllo convergenza (norma della differenza tra due iterazioni)
        if np.linalg.norm(x_new - x_old) < tolleranza:
            print(f"Convergenza raggiunta in {k + 1} iterazioni.")
            return x_new

        # Aggiorniamo x_old per la prossima iterazione
        x_old = x_new.copy()

    print("Numero massimo di iterazioni raggiunto.")
    return x_new