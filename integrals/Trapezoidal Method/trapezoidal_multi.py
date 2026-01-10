"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo dei Trapezi (Multiple-Segment)
"""

def metodo_trapezi_multi(func, a, b, n):
    """
    func: La funzione da integrare
    a: Inizio dell'intervallo
    b: Fine dell'intervallo
    n: Numero di intervalli
    """
    # Calcolo del passo (h)
    h = (b - a) / n

    # Valutazione estremi (f(a) + f(b))
    somma = func(a) + func(b)

    # Sommatoria dei punti interni moltiplicati per 2
    for i in range(1, n):
        x = a + i * h
        somma += 2 * func(x)

    # Calcolo finale
    return (h / 2) * somma

# --- MAIN DI ESEMPIO ---
if __name__ == "__main__":
    # 1. DEFINIZIONE DELLA FUNZIONE
    def funzione(x):
        return x ** 2

    # 2. IMPOSTAZIONE DEI PARAMETRI
    print("--- Calcolo Integrale col Metodo dei Trapezi ---")

    a = 0  # Estremo inferiore
    b = 10  # Estremo superiore
    n = 100  # Numero di intervalli (più è alto, più è preciso)

    # 3. ESECUZIONE DEL CALCOLO
    risultato_stimato = metodo_trapezi_multi(funzione_da_integrare, a, b, n)

    # 4. CALCOLO VALORE ESATTO (Solo per confronto matematico)
    # L'integrale di x^2 è x^3/3. Calcoliamo F(b) - F(a)
    valore_esatto = (b ** 3 / 3) - (a ** 3 / 3)
    errore = abs(valore_esatto - risultato_stimato)

    # 5. STAMPA DEI RISULTATI
    print(f"Intervallo: [{a}, {b}]")
    print(f"Numero intervalli (n): {n}")
    print("-" * 30)
    print(f"Risultato Metodo Trapezi: {risultato_stimato:.5f}")
    print(f"Valore Esatto (Analitico): {valore_esatto:.5f}")
    print(f"Errore assoluto:          {errore:.5f}")