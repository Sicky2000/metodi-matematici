"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Gauss-Seidel con Fattore di Rilassamento
"""

def gseid(a, b, n, x, max_iter, tol, lam):
    """
    a        : Matrice dei coefficienti (lista di liste n x n)
    b        : Vettore dei termini noti (lista di dimensione n)
    n        : Numero di equazioni/incognite
    x        : Vettore per la stima iniziale (viene aggiornato con la soluzione)
    max_iter : Numero massimo di iterazioni consentite
    tol      : Tolleranza per l'errore percentuale (criterio di arresto)
    lam      : Fattore di rilassamento (lambda)
    """
    # Normalizzazione delle righe
    # Divide ogni riga per l'elemento diagonale affinché a[i][i] diventi 1.
    for i in range(n):
        dummy = a[i][i]
        for j in range(n):
            a[i][j] = a[i][j] / dummy
        b[i] = b[i] / dummy

    # Stima iniziale
    # Prima passata per calcolare un valore iniziale migliore per x
    for i in range(n):
        somma = b[i]
        for j in range(n):
            if i != j:
                somma = somma - a[i][j] * x[j]
        x[i] = somma

    # Ciclo iterativo principale
    iter_count = 1

    while True:
        sentinel = 1  # Flag "ottimista": assume convergenza (1). Se un errore > tol, diventa 0 per continuare.

        for i in range(n):
            old = x[i]  # Memorizza il vecchio valore di x[i]
            somma = b[i]

            for j in range(n):
                if i != j:
                    somma = somma - a[i][j] * x[j]

            # Applica la formula SOR (Rilassamento)
            # x_nuovo = lambda * calcolato + (1 - lambda) * vecchio
            x[i] = lam * somma + (1.0 - lam) * old

            # Controllo dell'errore
            # Calcoliamo l'errore solo se sentinel è ancora 1 e x[i] non è zero
            if sentinel == 1 and x[i] != 0:
                ea = abs((x[i] - old) / x[i]) * 100
                if ea > tol:
                    sentinel = 0  # L'errore è troppo alto, continuiamo a iterare

        iter_count += 1

        # Condizione di uscita:
        # Se sentinel è rimasto 1 (errore sotto la soglia per tutti)
        # OPPURE abbiamo raggiunto il massimo delle iterazioni.
        if sentinel == 1 or iter_count >= max_iter:
            break

    return x, iter_count

# --- MAIN DI ESEMPIO

if __name__ == "__main__":
    # SISTEMA
    # 4x - y      = 2
    # -x + 4y - z = 6
    #      -y + 4z = 2
    # Soluzione attesa: x=1, y=2, z=1

    # 1. DEFINIZIONE MATRICE A (COEFFICIENTI)
    # Nota: Usiamo float (es. 4.0) per evitare divisioni intere accidentali
    a = [
        [4.0, -1.0, 0.0],
        [-1.0, 4.0, -1.0],
        [0.0, -1.0, 4.0]
    ]

    # 2. DEFINIZIONE VETTORE b (TERMINI NOTI)
    b = [2.0, 6.0, 2.0]

    # 3. PARAMETRI DI CONFIGURAZIONE
    n = 3  # Numero di equazioni
    x = [0.0, 0.0, 0.0]  # Stima iniziale (tutto a zero)
    max_iter = 100  # Numero massimo di iterazioni
    tol = 0.0001  # Tolleranza errore (%)
    lam = 1.0  # Lambda (1.0 = G-S standard, 1.1 = Sovra-rilassamento)

    # Stampa dei dati di input
    print(f"Numero equazioni: {n}")
    print(f"Lambda (rilassamento): {lam}")
    print(f"Tolleranza errore: {tol}%")
    print("-" * 30)

    # 4. CHIAMATA ALLA FUNZIONE
    soluzione, iter_count = gseid(a, b, n, x, max_iter, tol, lam)

    # 5. OUTPUT RISULTATI
    print("\n--- Risultati finali ---")
    if iter_count >= max_iter:
        print("ATTENZIONE: Il metodo NON ha raggiunto la convergenza nel numero massimo di iterazioni.")
    else:
        print("Convergenza raggiunta!")

    print(f"Iterazioni impiegate: {iter_count}")

    # Formattiamo la stampa del vettore soluzione per renderlo leggibile
    print("Vettore soluzione x:")
    for i, val in enumerate(soluzione):
        print(f"  x[{i}] = {val:.6f}")