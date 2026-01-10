"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Bisezione
"""

def bisect(func, xl, xu, tol, max_iter):

    """
        func = funzione f(x)
        xl = estremo inferiore dell'intervallo
        xu = estremo superiore dell'intervallo
        tol = tolleranza
        max_iter = numero massimo di iterazioni
    """

    if func(xl) * func(xu) >= 0:
        print("Errore: Lo zero non è garantito nell'intervallo. func(xl) e func(xu) devono avere segni opposti.")
        return None, 0

    i = 0 # contatore delle iterazioni
    xr = 0.0 # stima dello zero
    ea = 1.0 # errore

    while i < max_iter: # esegui il ciclo finché non è raggiunto il limite massimo di iterazioni
        xrold = xr # xrold memorizza il valore precedente dello zero stimato
        xr = (xl + xu) / 2 # calcolo la semisomma tra i due estremi e sovrascrivo xr
        i += 1 # incremento il contatore di 1

        if xr != 0 and i > 1:
            # errore relativo (calcolato dalla seconda iterazione in poi)
            ea = abs((xr - xrold) / xr)
        else:
            # fallback se la soluzione è 0 (per evitare divisione per zero nell'errore)
            ea = abs(xu - xl)

        test_val = func(xl) * func(xr) # criterio di selezione del nuovo intervallo

        if test_val < 0: # segno opposto, lo zero si trova nella metà sinistra (tra xl e xr)
            xu = xr # restringo l'intervallo sostituendo l'estremo superiore con il punto medio calcolato
        elif test_val > 0: # stesso segno, lo zero si trova nella metà destra (tra xr e xu)
            xl = xr # restringo l'intervallo sostituendo l'estremo inferiore con il punto medio calcolato
        else: # caso speciale: f(xr) è esattamente 0
            ea = 0.0 # imposto l'errore a 0 per forzare l'uscita dal ciclo

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
        return x ** 3 - x - 2

    # 2. DEFINIZIONE DEI PARAMETRI
    x_lower = 1.0        # Estremo inferiore
    x_upper = 2.0        # Estremo superiore
    tolleranza = 1e-5    # Tolleranza desiderata (0.00001)
    max_iterazioni = 100 # Limite iterazioni

    print(f"Intervallo: [{x_lower}, {x_upper}]")
    print(f"Tolleranza: {tolleranza}")
    print(f"Funzione: x^3 - x - 2\n")

    # 3. CHIAMATA ALLA FUNZIONE BISECT
    radice, iterazioni = bisect(funzione, x_lower, x_upper, tolleranza, max_iterazioni)

    # 4. OUTPUT DEI RISULTATI
    if radice is not None:
        print("-" * 30)
        print(f"Risultato finale:")
        # Formattazione a 6 cifre decimali
        print(f"Radice stimata (xr): {radice:.6f}")
        print(f"Iterazioni totali:   {iterazioni}")
        print(f"Valore di f(xr):     {funzione(radice):.2e}") # Verifica quanto siamo vicini a 0
        print("-" * 30)
    else:
        print("\nImpossibile trovare la radice: intervallo non valido.")