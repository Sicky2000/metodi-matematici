"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Interazione di Punto Fisso
"""

def fixed_point(func, x0, tol, max_iter):

    """
        func = funzione di iterazione g(x) (tale che x = g(x))
        x0 = stima iniziale dello zero
        tol = tolleranza
        max_iter = numero massimo di iterazioni
    """

    xr = x0 # xr memorizza la stima iniziale dello zero
    i = 0 # contatore delle iterazioni
    ea = 1.0 # errore

    while i < max_iter: # esegui il ciclo finché non è raggiunto il numero massimo di iterazioni
        xrold = xr # xrold memorizza la stima precedente dello zero
        xr = func(xrold) # calcolo il valore della funzione nel punto di stima dello zero
        i+=1 # incremento il contatore di 1

        if xr != 0: # se xr è diverso da 0
            ea = abs((xr - xrold) / xr) # calcolo l'errore relativo
        else:
            # fallback se la soluzione è 0 (per evitare divisione per zero nell'errore)
            ea = abs(xr - xrold)

        # controllo se l'errore è sceso sotto la tolleranza (SUCCESSO)
        if ea < tol:
            print(f"Radice trovata con successo.")
            break

        # controllo se ho raggiunto il numero massimo di iterazioni (FALLIMENTO)
        elif i >= max_iter:
            print("Numero massimo di iterazioni raggiunto senza successo.")
            break

    return xr, i

# --- MAIN DI ESEMPIO ---
if __name__ == "__main__":
    # 1. DEFINIAMO LA FUNZIONE
    def funzione(x):
        return (x + 2) ** (1 / 3)

    # 2. PARAMETRI
    x0 = 1.5  # Stima iniziale (tra i precedenti estremi 1 e 2)
    tol = 0.001  # Tolleranza
    max_iter = 20  # Interazioni massime

    root, iter_count = fixed_point(funzione, x0, tol, max_iter)
    print(f"Valore: {root:.6f} trovato in {iter_count} iterazioni")