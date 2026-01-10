"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Falsa Posizione (Variante Illinois)
"""

def false_pos(func, xl, xu, tol, max_iter):

    """
        func = funzione f(x)
        xl = estremo inferiore dell'intervallo
        xu = estremo superiore dell'intervallo
        tol = tolleranza
        max_iter = numero massimo di iterazioni
    """

    i = 0 # contatore delle iterazioni
    xr = 0.0 # stima dello zero
    ea = 1.0 # errore
    fl = func(xl) # funzione valutata nell'estremo inferiore
    fu = func(xu) # funzione valutata nell'estremo superiore
    iu = 0 # indice superiore
    il = 0 # indice inferiore

    if fl * fu >= 0:
        print("Errore: Lo zero non è garantito nell'intervallo. func(xl) e func(xu) devono avere segni opposti.")
        return None, 0

    while i < max_iter: # esegui il ciclo finché non è raggiunto il limite massimo di iterazioni
        xrold = xr # xrold memorizza la precedente stima dello zero
        xr = xu - fu * (xl - xu) / (fl - fu) # calcolo il nuovo valore dello zero con la formula di falsa posizione
        fr = func(xr) # calcolo il valore della funzione nel punto xr
        i += 1 # incremento il contatore di 1

        if xr != 0 and i > 1:
            # errore relativo (calcolato dalla seconda iterazione in poi)
            ea = abs((xr - xrold) / xr)
        else:
            # fallback se la soluzione è 0 (per evitare divisione per zero nell'errore)
            ea = abs(xr - xrold)

        test_val = fl * fr # criterio di selezione del nuovo intervallo

        if test_val < 0: # segno opposto, lo zero si trova nella metà sinistra (tra xl e xr)
            xu = xr # restringo l'intervallo sostituendo l'estremo superiore
            fu = fr # aggiorno il valore della funzione nel nuovo estremo
            iu = 0 # reset contatore: xu si è appena mosso, quindi non è stagnante
            il += 1 # incremento contatore: xl non si è mosso, quindi aumenta la sua "anzianità".

            if il >= 2: # se l'estremo inferiore (xl) è rimasto fermo per 2 o più iterazioni consecutive:
                fl = fl / 2 # "penalizzo" il valore f(xl) dimezzandolo.

        elif test_val > 0: # stesso segno, lo zero si trova nella metà destra (tra xr e xu)
            xl = xr # restringo l'intervallo sostituendo l'estremo inferiore
            fl = fr # aggiorno il valore della funzione nel nuovo estremo
            il = 0 # reset contatore: xl si è appena mosso, quindi non è stagnante
            iu += 1 # incremento contatore: xu non si è mosso, quindi aumenta la sua "anzianità".

            if iu >= 2: # se l'estremo superiore (xu) è rimasto fermo per 2 o più iterazioni consecutive:
                fu = fu / 2 # "penalizzo" il valore f(xu) dimezzandolo.

        else: # caso speciale: f(xr) è esattamente 0
            ea = 0 # imposto l'errore a 0 per forzare l'uscita dal ciclo

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
    # DEFINIZAMO LA FUNZIONE
    def funzione(x):
        return x**3 - x - 2

    root, iterations = false_pos(funzione, 1, 2, 0.001, 20)
    print(f"Valore: {root:.6f} trovato in {iterations} iterazioni")