"""
Modulo per l'Interpolazione Polinomiale.

Include metodi per trovare il polinomio che passa per un set di punti (nodi)
e funzioni di utilità per la scelta ottimale dei nodi (Chebyshev).

Autore:      Sicky2005
Corso:       Metodi Numerici per l'Ingegneria
"""

import numpy as np

def lagrange(x_nodes, y_nodes, x_target):
    """
    Calcola il valore interpolato in x_target usando il Polinomio di Lagrange.

    Args:
        x_nodes (list/array): Coordinate x dei nodi.
        y_nodes (list/array): Coordinate y dei nodi.
        x_target (float): Il punto in cui valutare l'interpolazione.

    Returns:
        float: Il valore interpolato y.
    """
    x = np.array(x_nodes, dtype=float)
    y = np.array(y_nodes, dtype=float)
    n = len(x)

    result = 0.0

    for i in range(n):
        L_i = 1.0
        for j in range(n):
            if i != j:
                L_i *= (x_target - x[j]) / (x[i] - x[j])
        result += y[i] * L_i

    return result

def newton(x_nodes, y_nodes, x_target):
    """
    Calcola il valore interpolato usando il metodo di Newton (Differenze Divise).
    Funziona con qualsiasi set di nodi (equidistanti, Chebyshev, casuali).

    Args:
        x_nodes (list/array): Coordinate x dei nodi.
        y_nodes (list/array): Coordinate y dei nodi.
        x_target (float): Il punto in cui valutare l'interpolazione.

    Returns:
        float: Il valore interpolato y.
    """
    x = np.array(x_nodes, dtype=float)
    y = np.array(y_nodes, dtype=float)
    n = len(x)

    # Inizializzazione matrice Differenze Divise
    fdd = np.zeros((n, n))
    fdd[:, 0] = y

    # Calcolo tabella FDD
    for j in range(1, n):
        for i in range(n - j):
            numerator = fdd[i + 1, j - 1] - fdd[i, j - 1]
            denominator = x[i + j] - x[i]
            fdd[i, j] = numerator / denominator

    # Coefficienti del polinomio (prima riga della tabella)
    coeffs = fdd[0, :]

    # Valutazione (Algoritmo di Horner)
    result = coeffs[0]
    x_term = 1.0

    for i in range(1, n):
        x_term *= (x_target - x[i - 1])
        result += coeffs[i] * x_term

    return result

def chebyshev_nodes(a, b, n):
    """
    Genera n nodi di Chebyshev nell'intervallo [a, b].
    Minimizza il fenomeno di Runge rispetto ai nodi equidistanti.

    Args:
        a (float): Estremo sinistro.
        b (float): Estremo destro.
        n (int): Numero di nodi.

    Returns:
        np.array: Array di n punti x.
    """
    # Generazione indici k = 0, 1, ..., n-1
    k = np.arange(n)

    # Nodi standard in [-1, 1]
    # Formula ottimizzata per k che parte da 0: cos((2k+1)/(2n)*pi)
    t_k = np.cos((2 * k + 1) / (2 * n) * np.pi)

    # Mappatura in [a, b]
    x_k = (a + b) / 2 + (b - a) / 2 * t_k

    return x_k