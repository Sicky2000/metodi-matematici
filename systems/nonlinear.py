"""
Modulo per la risoluzione di Sistemi di Equazioni Non Lineari.

Autore:      Sicky2005
Corso:       Metodi Numerici per l'Ingegneria
"""

import numpy as np

import numpy as np


def newton_system(F, J, x0, tol=1e-6, max_iter=50):
    """
    Metodo di Newton-Raphson per sistemi di equazioni non lineari.
    Risolve J(x) * delta = -F(x)

    Args:
        F (callable): Funzione che accetta un vettore x e restituisce un vettore F(x).
        J (callable): Funzione che accetta un vettore x e restituisce la matrice Jacobiana J(x).
        x0 (list/array): Vettore delle stime iniziali.
        tol (float): Tolleranza per l'arresto (norma dell'errore).
        max_iter (int): Numero massimo di iterazioni.

    Returns:
        numpy.ndarray: Il vettore soluzione x.
    Raises:
        ValueError: Se la matrice Jacobiana è singolare o il metodo non converge.
    """
    x = np.array(x0, dtype=float)

    for k in range(max_iter):
        # 1. Valutiamo le funzioni e lo Jacobiano nel punto corrente
        Fx = np.array(F(x))
        Jx = np.array(J(x))

        # 2. Controllo convergenza (Norma del residuo)
        if np.linalg.norm(Fx) < tol:
            return x

        # 3. Risolviamo il sistema lineare J * delta = -F
        try:
            delta = np.linalg.solve(Jx, -Fx)
        except np.linalg.LinAlgError:
            raise ValueError("La matrice Jacobiana è singolare (determinante 0). Il metodo fallisce.")

        # 4. Aggiorniamo la soluzione
        x = x + delta

        # 5. Controllo convergenza sullo step
        if np.linalg.norm(delta) < tol:
            return x

    print("Attenzione: Numero massimo di iterazioni raggiunto.")
    return x

def broyden(f, x0, tol=1e-6, max_iter=100, B0=None):
    """
    Risolve un sistema di equazioni non lineari f(x) = 0 usando il metodo di Broyden.
    È un metodo "Quasi-Newton" che aggiorna un'approssimazione dello Jacobiano
    invece di ricalcolarlo ad ogni passo.

    Args:
        f (callable): Funzione che accetta un array x e restituisce un array f(x).
        x0 (list/array): Vettore tentativo iniziale.
        tol (float): Tolleranza per la norma del passo (o residuo).
        max_iter (int): Numero massimo di iterazioni.
        B0 (list/array, optional): Stima iniziale dello Jacobiano.
                                   Se None, usa la matrice Identità.

    Returns:
        np.array: Il vettore soluzione x.

    Raises:
        np.linalg.LinAlgError: Se la matrice B diventa singolare.
        RuntimeError: Se il metodo non converge.
    """
    x = np.array(x0, dtype=float).flatten()
    n = len(x)

    # Inizializzazione Jacobiano approssimato
    if B0 is None:
        B = np.eye(n)
    else:
        B = np.array(B0, dtype=float)

    # Valutazione iniziale
    fx = np.array(f(x), dtype=float).flatten()

    for i in range(max_iter):
        # 1. Risolvi il sistema lineare B * delta_x = -f(x)
        try:
            delta_x = np.linalg.solve(B, -fx)
        except np.linalg.LinAlgError:
            raise np.linalg.LinAlgError("L'approssimazione dello Jacobiano è singolare. Il metodo fallisce.")

        # 2. Aggiorna x
        x_new = x + delta_x
        fx_new = np.array(f(x_new), dtype=float).flatten()

        # 3. Controllo convergenza (sulla norma del passo o del residuo)
        if np.linalg.norm(delta_x) < tol:
            return x_new

        # 4. Aggiornamento di Broyden (Formula ottimizzata)
        # B_new = B + (f(x_new) * delta_x^T) / ||delta_x||^2
        # Nota: Usiamo fx_new perché (y - B*s) si semplifica in fx_new
        delta_norm_sq = np.dot(delta_x, delta_x)

        if delta_norm_sq > 1e-16:
            B += np.outer(fx_new, delta_x) / delta_norm_sq
        else:
            # Se il passo è minuscolo ma non siamo ancora a convergenza,
            # potremmo resettare B all'identità o uscire. Qui continuiamo.
            pass

        # Aggiornamento variabili per il prossimo giro
        x = x_new
        fx = fx_new

    raise RuntimeError(f"Il metodo di Broyden non ha convertito dopo {max_iter} iterazioni.")