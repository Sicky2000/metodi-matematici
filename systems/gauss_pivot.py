"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Gauss con il Pivoting
"""

import numpy as np

def pivot(a, b, s, n, k):
    """
      a: Matrice dei coefficienti (viene modificata in-place scambiando le righe).
      b: Vettore dei termini noti (viene modificato in-place).
      s: Vettore dei fattori di scala (viene modificato in-place).
      n: Dimensione del sistema (numero di equazioni).
      k: Indice della colonna corrente su cui stiamo lavorando.
    """
    # usiamo riga_pivot per indicare l'indice della riga migliore trovata
    riga_pivot = k

    # big = rapporto massimo trovato finora
    big = abs(a[k, k] / s[k])

    # cerca la riga migliore tra quelle sottostanti (i da k+1 a n)
    for i in range(k + 1, n):
        dummy = abs(a[i, k] / s[i])

        if dummy > big:
            big = dummy
            riga_pivot = i

    # se la riga migliore (riga_pivot) è diversa da quella corrente (k), scambiamo
    if riga_pivot != k:
        # scambia righe in a
        a[[riga_pivot, k]] = a[[k, riga_pivot]]
        # scambia elementi in b
        b[riga_pivot], b[k] = b[k], b[riga_pivot]
        # scambia fattori di scala s
        s[riga_pivot], s[k] = s[k], s[riga_pivot]


def eliminazione(a, s, n, b, tol):
    """
      a: Matrice dei coefficienti.
      s: Vettore dei fattori di scala.
      n: Dimensione del sistema.
      b: Vettore dei termini noti.
      tol: Valore di tolleranza per determinare la singolarità.
    """
    # ciclo k da 0 a n-2
    for k in range(0, n - 1):
        pivot(a, b, s, n, k)

        # controllo tolleranza: se il pivot è quasi zero, ea = -1.0
        if abs(a[k, k] / s[k]) < tol:
            return -1.0

        # eliminazione sulle righe sottostanti
        for i in range(k + 1, n):
            factor = a[i, k] / a[k, k]

            for j in range(k + 1, n):
                a[i, j] = a[i, j] - factor * a[k, j]

            b[i] = b[i] - factor * b[k]

    # controllo finale sull'ultimo elemento
    if abs(a[n - 1, n - 1] / s[n - 1]) < tol:
        return -1.0

    return 0.0


def sostituzione(a, n, b):
    """
      a: Matrice triangolare superiore (dopo l'eliminazione).
      n: Dimensione del sistema.
      b: Vettore dei termini noti (modificato dall'eliminazione).
    """
    xr = np.zeros(n)

    # calcolo ultimo elemento
    xr[n - 1] = b[n - 1] / a[n - 1, n - 1]

    # ciclo all'indietro (da n-2 a 0)
    for i in range(n - 2, -1, -1):
        sum_val = 0
        for j in range(i + 1, n):
            sum_val = sum_val + a[i, j] * xr[j]

        xr[i] = (b[i] - sum_val) / a[i, i]

    return xr


def gauss_pivot(matrice_in, b_in, tol=1e-6):
    """
      matrice_in : Matrice dei coefficienti originale (lista di liste o array).
      b_in       : Vettore dei termini noti originale.
      tol        : Tolleranza per il pivot (default 1e-6).
    """
    # copia dei dati
    a = np.array(matrice_in, dtype=float)
    b = np.array(b_in, dtype=float)
    n = len(b)
    s = np.zeros(n)

    # calcolo fattori di scala (s)
    for i in range(n):
        s[i] = abs(a[i, 0])
        for j in range(1, n):
            if abs(a[i, j]) > s[i]:
                s[i] = abs(a[i, j])

    # eliminazione
    ea = eliminazione(a, s, n, b, tol)

    if ea != -1.0:
        # STAMPA ITERAZIONI (Passi di eliminazione)
        print(f"Iterazioni (passi di eliminazione): {n - 1}")

        # sostituzione
        xr = sostituzione(a, n, b)
        return xr
    else:
        print("Errore: Il sistema è singolare (ea = -1.0).")
        return None

# --- MAIN DI ESEMPIO ---
if __name__ == "__main__":
    # SISTEMA
    # 3x - 0.1y - 0.2z = 7.85
    # 0.1x + 7y - 0.3z = -19.3
    # 0.3x - 0.2y + 10z = 71.4

    A_input = [
        [3.0, -0.1, -0.2],
        [0.1, 7.0, -0.3],
        [0.3, -0.2, 10.0]
    ]

    b_input = [7.85, -19.3, 71.4]

    print("--- Inizio Calcolo ---")
    risultato = gauss_pivot(A_input, b_input)

    if risultato is not None:
        print("\nIl vettore soluzione x è:")
        print(risultato)