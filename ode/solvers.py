"""
Modulo per la risoluzione di Equazioni Differenziali Ordinarie (ODE).

Include metodi a passo singolo espliciti:
- 1° Ordine: Eulero
- 2° Ordine (RK2): Heun (Semplice e Iterativo), Midpoint, Ralston
- 4° Ordine (RK4): Classico (gestisce anche sistemi di ODE)

Autore:      Sicky2005
Corso:       Metodi Numerici per l'Ingegneria
"""

import numpy as np

def euler(f, x0, y0, x_end, h):
    """
    Risolve una ODE o un sistema di ODE usando il Metodo di Eulero (1° Ordine).
    Formula: y_{i+1} = y_i + f(x_i, y_i) * h

    Args:
        f (callable): Funzione derivata dy/dx = f(x, y). Deve accettare (x, y) e ritornare float o np.ndarray.
        x0 (float): Valore iniziale della variabile indipendente.
        y0 (float | np.ndarray): Condizione iniziale y(x0). Scalare per singola ODE, array per sistemi.
        x_end (float): Valore finale della variabile indipendente.
        h (float): Passo di integrazione.

    Returns:
        tuple[np.ndarray, np.ndarray]: Una tupla (x, y) dove:
            - x è l'array dei passi temporali.
            - y è l'array delle soluzioni (o matrice n_passi x n_variabili per i sistemi).
    """
    x = np.arange(x0, x_end + h/100, h)
    n = len(x)

    y0 = np.array(y0, dtype=float)
    if y0.ndim == 0:
        y = np.zeros(n)
    else:
        y = np.zeros((n, len(y0)))

    y[0] = y0

    for i in range(n - 1):
        y[i + 1] = y[i] + f(x[i], y[i]) * h

    return x, y


def heun(f, x0, y0, x_end, h):
    """
    Risolve una ODE usando il Metodo di Heun Semplice (RK2 - Senza iterazione).
    Usa la media tra la pendenza iniziale e quella stimata alla fine dell'intervallo.

    Args:
        f (callable): Funzione derivata dy/dx = f(x, y).
        x0 (float): Valore iniziale x.
        y0 (float | np.ndarray): Valore iniziale y.
        x_end (float): Valore finale x.
        h (float): Passo di integrazione.

    Returns:
        tuple[np.ndarray, np.ndarray]: Tupla (x, y) con i risultati.
    """
    x = np.arange(x0, x_end + h/100, h)
    n = len(x)

    y0 = np.array(y0, dtype=float)
    if y0.ndim == 0:
        y = np.zeros(n)
    else:
        y = np.zeros((n, len(y0)))

    y[0] = y0

    for i in range(n - 1):
        k1 = f(x[i], y[i])                # Pendenza iniziale
        k2 = f(x[i] + h, y[i] + k1 * h)   # Pendenza finale (stimata con Eulero)
        y[i + 1] = y[i] + (0.5 * k1 + 0.5 * k2) * h

    return x, y


def heun_iterative(f, x0, y0, x_end, h, es=0.01, max_it=20):
    """
    Risolve una ODE usando il Metodo di Heun con Iterazione (Correttore).
    Migliora la stima correggendo y_{i+1} finché l'errore approssimato non scende sotto la soglia.

    Args:
        f (callable): Funzione derivata dy/dx = f(x, y).
        x0 (float): Valore iniziale x.
        y0 (float | np.ndarray): Valore iniziale y.
        x_end (float): Valore finale x.
        h (float): Passo di integrazione.
        es (float, optional): Errore relativo target in percentuale (%). Default 0.01.
        max_it (int, optional): Numero massimo di iterazioni del correttore per ogni passo. Default 20.

    Returns:
        tuple[np.ndarray, np.ndarray]: Tupla (x, y) con i risultati raffinati.
    """
    x = np.arange(x0, x_end + h/100, h)
    n = len(x)

    y0 = np.array(y0, dtype=float)
    if y0.ndim == 0:
        y = np.zeros(n)
    else:
        y = np.zeros((n, len(y0)))

    y[0] = y0

    for i in range(n - 1):
        # 1. Predittore (Eulero standard)
        slope_predictor = f(x[i], y[i])
        y_predict = y[i] + slope_predictor * h

        # 2. Correttore Iterativo
        y_old = y_predict
        for it in range(max_it):
            # Calcoliamo la pendenza nel punto predetto
            slope_end = f(x[i] + h, y_old)

            # Pendenza media
            slope_avg = (slope_predictor + slope_end) / 2.0

            # Nuova stima di y
            y_new = y[i] + slope_avg * h

            # Calcolo errore relativo approssimato (ea)
            # Aggiungiamo epsilon per evitare divisioni per zero
            ea = np.max(np.abs((y_new - y_old) / (y_new + 1e-15))) * 100

            y_old = y_new

            # Uscita anticipata se l'errore è sotto la soglia
            if ea <= es:
                break

        # Salviamo il risultato raffinato
        y[i + 1] = y_old

    return x, y


def midpoint(f, x0, y0, x_end, h):
    """
    Risolve una ODE usando il Metodo del Punto Medio (Midpoint - RK2).
    Valuta la pendenza a metà dell'intervallo h/2.

    Args:
        f (callable): Funzione derivata dy/dx = f(x, y).
        x0 (float): Valore iniziale x.
        y0 (float | np.ndarray): Valore iniziale y.
        x_end (float): Valore finale x.
        h (float): Passo di integrazione.

    Returns:
        tuple[np.ndarray, np.ndarray]: Tupla (x, y) con i risultati.
    """
    x = np.arange(x0, x_end + h/100, h)
    n = len(x)

    y0 = np.array(y0, dtype=float)
    if y0.ndim == 0:
        y = np.zeros(n)
    else:
        y = np.zeros((n, len(y0)))

    y[0] = y0

    for i in range(n - 1):
        k1 = f(x[i], y[i])
        k2 = f(x[i] + 0.5 * h, y[i] + 0.5 * k1 * h)
        y[i + 1] = y[i] + k2 * h

    return x, y


def ralston(f, x0, y0, x_end, h):
    """
    Risolve una ODE usando il Metodo di Ralston (RK2 Ottimizzato).
    Metodo di secondo ordine che minimizza l'errore di troncamento.

    Args:
        f (callable): Funzione derivata dy/dx = f(x, y).
        x0 (float): Valore iniziale x.
        y0 (float | np.ndarray): Valore iniziale y.
        x_end (float): Valore finale x.
        h (float): Passo di integrazione.

    Returns:
        tuple[np.ndarray, np.ndarray]: Tupla (x, y) con i risultati.
    """
    x = np.arange(x0, x_end + h/100, h)
    n = len(x)

    y0 = np.array(y0, dtype=float)
    if y0.ndim == 0:
        y = np.zeros(n)
    else:
        y = np.zeros((n, len(y0)))

    y[0] = y0

    for i in range(n - 1):
        k1 = f(x[i], y[i])
        k2 = f(x[i] + 0.75 * h, y[i] + 0.75 * k1 * h)
        y[i + 1] = y[i] + ((1/3) * k1 + (2/3) * k2) * h

    return x, y


def rk4(f, x0, y0, x_end, h):
    """
    Risolve una ODE o un Sistema di ODE usando Runge-Kutta del 4° Ordine (RK4).
    Metodo standard ad alta precisione (O(h^4)). Gestisce automaticamente input scalari o vettoriali.

    Args:
        f (callable): Funzione derivata dy/dx = f(x, y).
                      Se sistema, deve ritornare un np.array delle derivate.
        x0 (float): Valore iniziale x.
        y0 (float | np.ndarray): Valore iniziale y.
                                 Passare una lista/array [y1, y2...] per i sistemi.
        x_end (float): Valore finale x.
        h (float): Passo di integrazione.

    Returns:
        tuple[np.ndarray, np.ndarray]: Una tupla (x, y) dove:
            - x è l'array dei passi temporali.
            - y è l'array delle soluzioni (o matrice n x m per i sistemi).
    """
    x = np.arange(x0, x_end + h/100, h)
    n = len(x)

    y0 = np.array(y0, dtype=float)

    if y0.ndim == 0:
        y = np.zeros(n)
    else:
        y = np.zeros((n, len(y0)))

    y[0] = y0

    for i in range(n - 1):
        k1 = f(x[i], y[i])
        k2 = f(x[i] + 0.5 * h, y[i] + 0.5 * h * k1)
        k3 = f(x[i] + 0.5 * h, y[i] + 0.5 * h * k2)
        k4 = f(x[i] + h, y[i] + h * k3)

        y[i + 1] = y[i] + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    return x, y