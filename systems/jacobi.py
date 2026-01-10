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

# --- MAIN DI ESEMPIO ---
if __name__ == "__main__":
    # 1. DEFINIZIONE DEL SISTEMA
    # Usiamo una matrice a diagonale dominante per garantire la convergenza.
    # Sistema:
    # 10x + 2y + z = 7
    # x + 5y + z = -8
    # 2x + 3y + 10z = 6

    A = np.array([
        [10.0, 2.0, 1.0],
        [1.0, 5.0, 1.0],
        [2.0, 3.0, 10.0]
    ])

    b = np.array([7.0, -8.0, 6.0])

    # 2. PARAMETRI INIZIALI
    # Ipotesi iniziale (vettore di zeri)
    x0 = np.zeros(len(b))
    tolleranza = 1e-6
    max_iter = 100

    # 3. CHIAMATA ALLA FUNZIONE
    soluzione_jacobi = jacobi(A, b, x0, tolleranza, max_iter)

    # 4. OUTPUT DEI RISULTATI
    print("-" * 30)
    print("Soluzione trovata con Jacobi:")
    print(soluzione_jacobi)
    print("-" * 30)

    # 5. VERIFICA (OPZIONALE) CON LA SOLUZIONE ESATTA DI NUMPY
    soluzione_esatta = np.linalg.solve(A, b)
    print("Soluzione esatta (numpy.linalg.solve):")
    print(soluzione_esatta)

    # Calcolo dell'errore
    errore = np.linalg.norm(soluzione_jacobi - soluzione_esatta)
    print(f"Errore rispetto alla soluzione esatta: {errore:.2e}")