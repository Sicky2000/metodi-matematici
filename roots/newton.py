"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Newton
"""

def newton(func, dfunc, x0, tol, max_iter):

    """
    func = funzione f(x)
    dfunc = derivata della funzione f'(x)
    x0 = stima iniziale
    tol = tolleranza
    max_iter = numero massimo di iterazioni
    """

    xr = x0 # xr memorizza la stima iniziale dello zero
    i = 0  # contatore delle iterazioni
    ea = 1.0  # errore approssimato relativo

    while i < max_iter: # eseguo il ciclo finché non è raggiunto il numero massimo di iterazioni
        xrold = xr # xrold memorizza il valore precedente dello zero stimato
        f_val = func(xrold) # calcolo il valore della funzione nel punto xr
        df_val = dfunc(xrold) # calcolo il valroe della derivata della funzione nel punto xr

        if df_val == 0:
            print("Errore: Derivata nulla. Divisione per zero impossibile.")
            return None, i

        xr = xrold - f_val / df_val # calcolo il nuovo valore dello zero con il metodo di Newton
        i += 1 # incremento il contatore di 1

        if xr != 0:
            # errore relativo
            ea = abs((xr - xrold) / xr)
        else:
            # fallback se la soluzione è 0 (per evitare divisione per zero nell'errore)
            ea = abs(xr - xrold)

        # controllo se l'errore è sceso sotto la tolleranza (SUCCESSO)
        if ea < tol:
            print(f"\n>> Radice trovata con successo in {i} iterazioni.")
            break

        # controllo se ho raggiunto il numero massimo di iterazioni (FALLIMENTO)
        elif i >= max_iter:
            print("\n>> Numero massimo di iterazioni raggiunto senza convergenza.")
            break

    return xr, i

# --- MAIN DI ESEMPIO ---
if __name__ == "__main__":
    # 1. DEFINIAMO LA FUNZIONE: x^2 - 4
    def funzione(x):
        return x ** 2 - 4

    # 2. DEFINIAMO LA DERIVATA DELLA FUNZIONE: 2x
    def derivata(x):
        return 2 * x

    result = newton(funzione, derivata, x0=6.0, tol=1e-6, max_iter=100)

    if result:
        root, iterations = result
        print(f"Valore: {root:.6f}")
    else:
        print("Soluzione non trovata (Derivata uguale a zero)")