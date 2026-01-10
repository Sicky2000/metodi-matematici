"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo delle Secanti
"""

def secant(func, x0, x1, tol, max_iter):

    """
    func = funzione f(x)
    x0 = prima stima iniziale
    x1 = seconda stima iniziale
    tol = tolleranza
    max_iter = numero massimo di iterazioni
    """

    f0 = func(x0) # calcolo il valore della funzione in x0
    f1 = func(x1) # calcolo il valore della funzione in x1

    i = 0  # contatore delle iterazioni
    xr = 0.0  # nuova stima della radice
    ea = 1.0  # errore relativo

    while i < max_iter: # eseguo il ciclo finché non raggiungo il numero massimo di iterazioni
        if abs(f1 - f0) < 1e-12:  # Uso una soglia piccola invece di == 0 per sicurezza numerica
            print("Errore: Il denominatore è nullo (divisione per zero). Secante orizzontale.")
            return None, i

        xr = x1 - f1 * (x1 - x0) / (f1 - f0) # calcolo la stima dello zero con il metodo delle secanti
        i += 1 # incremento il contatore di 1

        if xr != 0:
            # errore relativo
            ea = abs((xr - x1) / xr)
        else:
            # fallback se la soluzione è 0 (per evitare divisione per zero nell'errore)
            ea = abs(xr - x1)  # Fallback errore assoluto

        # controllo se l'errore è sceso sotto la tolleranza (SUCCESSO)
        if ea < tol:
                print(f"\n>> Radice trovata con successo in {i} iterazioni.")
                break

        # controllo se ho raggiunto il numero massimo di iterazioni (FALLIMENTO)
        if i >= max_iter:
            print("\n>> Numero massimo di iterazioni raggiunto senza convergenza.")
            break

        x0 = x1 # x(i-1) prende il valore di x(i)
        f0 = f1 # Aggiorno f(x(i-1))

        x1 = xr # x(i) prende il valore di x(i+1) appena calcolato
        f1 = func(xr) # Calcolo f(x(i)) per il prossimo ciclo

    return xr, i

# --- MAIN DI ESEMPIO ---
if __name__ == "__main__":
    # DEFINIAMO LA FUNZIONE: f(x) = x^2 - 2
    def f_test(x):
        return x**2 - 2

    # Stimiamo valori compresi tra 0 e 10 (abbastanza lontani per testare la stabilità del metodo)
    root, iterations = secant(f_test, 0, 10, 1e-6, 100)

    if root is not None:
        print(f"Valore: {root:.8f}")