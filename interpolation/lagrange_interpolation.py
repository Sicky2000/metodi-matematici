"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Interpolazione di Lagrange
"""

import numpy as np

def interpolazione_lagrange(x_nodi, y_nodi, x_target):
    """
    x_nodi: Array dei punti x conosciuti (nodi)
    y_nodi: Array dei valori y corrispondenti
    x_target: Il punto x in cui vogliamo trovare il valore y
    """
    # ricavo n dalla lunghezza del vettore
    n = len(x_nodi)

    # inizializzo la somma finale a 0
    risultato_interpolazione = 0.0

    # Ciclo esterno
    for i in range(n):

        # inizializzo il prodotto base Li(x) a 1 (elemento neutro della moltiplicazione)
        polinomio_base_li = 1.0

        # ciclo interno
        for j in range(n):
            if i != j:
                numeratore = x_target - x_nodi[j]
                denominatore = x_nodi[i] - x_nodi[j]
                polinomio_base_li = polinomio_base_li * (numeratore / denominatore)

        # ora moltiplico il polinomio base Li per il valore yi
        termine_corrente = y_nodi[i] * polinomio_base_li

        # aggiungo alla somma totale
        risultato_interpolazione += termine_corrente

        # mostro quanto vale Li e quanto contribuisce questo nodo al totale
        print(f"Nodo {i} (x={x_nodi[i]}): L_{i} = {polinomio_base_li:.5f} -> Contributo (+ y[{i}]*L_{i}) = {termine_corrente:.10f}")

    return risultato_interpolazione

# --- MAIN DI ESEMPIO ---
if __name__ == "__main__":
    # 1. DEFINIZIONE DELLA FUNZIONE
    def funzione(x):
        return np.sin(x)

    # 2. DEFINIAMO I NODI X A PIACERE
    # Usiamo 0, pi/2 e pi per interpolare il seno
    x_dati = [0.0, 1.5, 3.0]

    # 3. CALCOLIAMO I NODI Y AUTOMATICAMENTE USANDO LA FUNZIONE DEFINITA SOPRA
    y_dati = [funzione(x) for x in x_dati]

    # 4. DEFINIAMO IL PUNTO DOVE VOGLIAMO LA STIMA
    punto_da_stimare = 1.0  # Vogliamo sapere quanto vale sin(1.0) interpolando

    print(f"Nodi X scelti: {x_dati}")
    print(f"Nodi Y calcolati (dalla funzione reale): {[round(y, 4) for y in y_dati]}")

    # 5. ESEGUIAMO L'INTERPOLAZIONE
    valore_interpolato = interpolazione_lagrange(x_dati, y_dati, punto_da_stimare)

    # 6. CALCOLIAMO IL VALORE VERO PER CONFRONTO
    valore_vero = funzione(punto_da_stimare)

    # 7. CALCOLO DELL'ERRORE ASSOLUTO
    errore_assoluto = abs(valore_vero - valore_interpolato)

    print("-" * 60)
    print(f"Risultato Lagrange (stima) : {valore_interpolato:.8f}")
    print(f"Risultato Esatto (f(x))    : {valore_vero:.8f}")
    print(f"Errore Assoluto            : {errore_assoluto:.8e}")
    print("-" * 60)