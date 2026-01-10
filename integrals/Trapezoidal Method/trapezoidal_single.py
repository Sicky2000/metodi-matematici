"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo dei Trapezi (Single-Segment)
"""

def metodo_trapezi_single(h, f0, f1):
    """
    h (float): La larghezza dell'intervallo (b - a).
    f0 (float): Il valore della funzione nel punto iniziale, f(a).
    f1 (float): Il valore della funzione nel punto finale, f(b).
    """
    return h * (f0 + f1) / 2

# --- ESEMPIO DI UTILIZZO ---
def funzione(x):
    """
    Definisce la funzione matematica di cui vogliamo calcolare l'area.
    In questo esempio usiamo: f(x) = x^2
    """
    return x ** 2

# --- MAIN ---
if __name__ == "__main__":

    # 1. Definiamo l'intervallo di integrazione [a, b]
    punto_a = 0.0  # Punto iniziale
    punto_b = 2.0  # Punto finale

    # 2. Calcoliamo h (la larghezza del segmento)
    # Questo corrisponde alla base del trapezio
    h = punto_b - punto_a

    # 3. Calcoliamo i valori della funzione agli estremi
    # Queste corrispondono alle due altezze parallele del trapezio
    f0 = funzione(punto_a)
    f1 = funzione(punto_b)

    # 4. Chiamiamo la funzione definita all'inizio per calcolare l'area
    area_calcolata = metodo_trapezi_single(h, f0, f1)

    # 5. Stampiamo i risultati a video
    print(f"Calcolo dell'integrale definito da {punto_a} a {punto_b}")
    print(f"Larghezza intervallo (h): {h}")
    print(f"Valori della funzione -> f({punto_a}): {f0}, f({punto_b}): {f1}")
    print("-" * 40)
    print(f"Risultato Area (Regola Trapezio): {area_calcolata}")