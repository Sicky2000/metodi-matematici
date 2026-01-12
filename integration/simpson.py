"""
Modulo per l'Integrazione Numerica (Metodo di Simpson).

Implementazione basata sull'algoritmo 'SimpInt' che gestisce dinamicamente
il numero di intervalli pari o dispari combinando Simpson 1/3, Simpson 3/8
e la regola dei Trapezi (per n=1).

Autore:      Sicky2005
Corso:       Metodi Numerici per l'Ingegneria
"""

import numpy as np
from .trapezoidal import trapezoidal


def simpson(f, a, b, n):
    """
    Calcola l'integrale definito di f(x) tra a e b.

    Logica:
    - Se n = 1: Usa la regola dei Trapezi.
    - Se n è DISPARI (>1): Usa Simpson 3/8 sugli ultimi 3 intervalli
      e Simpson 1/3 sui rimanenti.
    - Se n è PARI: Usa Simpson 1/3 su tutto l'intervallo.

    Args:
        f (callable): Funzione integranda.
        a (float): Estremo inferiore.
        b (float): Estremo superiore.
        n (int): Numero di intervalli.

    Returns:
        float: Valore approssimato dell'integrale.
    """
    # Controllo input base
    if n < 1:
        raise ValueError("Il numero di intervalli n deve essere almeno 1.")

    h = (b - a) / n

    # CASO n=1:
    if n == 1:
        # Richiamiamo la funzione trapezoidal
        return trapezoidal(f, a, b, n)

    # Generiamo tutti i nodi e le valutazioni della funzione
    x = np.linspace(a, b, n + 1)
    y = f(x)

    sum_val = 0.0
    m = n

    # CASO DISPARI:
    if n % 2 != 0:
        # Applichiamo Simpson 3/8 agli ultimi 3 intervalli (4 punti)
        # Punti coinvolti: n-3, n-2, n-1, n
        # Formula Simp38 (pannello b): 3h * (f0 + 3f1 + 3f2 + f3) / 8
        sum_val += (3 * h / 8) * (y[n - 3] + 3 * y[n - 2] + 3 * y[n - 1] + y[n])

        # Riduciamo m di 3, così trattiamo la parte restante con la 1/3
        m = n - 3

    # CASO PARI (o residuo del dispari):
    if m > 1:
        # Applichiamo Simpson 1/3 Multiplo (Simp13m ) sui primi m intervalli
        # Slice dei dati fino all'indice m
        # Formula vettorializzata di Simp13m:
        # h/3 * (y0 + 4*(dispari) + 2*(pari) + ym)

        # Somma termini indici dispari (1, 3, ..., m-1)
        sum_odds = 4 * np.sum(y[1:m:2])

        # Somma termini indici pari (2, 4, ..., m-2)
        sum_evens = 2 * np.sum(y[2:m:2])

        simpson_13_part = (h / 3) * (y[0] + sum_odds + sum_evens + y[m])

        sum_val += simpson_13_part

    return sum_val