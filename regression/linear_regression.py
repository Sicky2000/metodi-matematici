"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per la Regressione con Generazione Grafico
"""

import matplotlib.pyplot as plt

def lin_reg(x, y):
    """
        x (list): Lista dei valori x
        y (list): Lista dei valori y
    """
    n = len(x)

    # Controllo di sicurezza: servono almeno 3 punti per calcolare syx (n-2)
    if n < 3:
        raise ValueError("Sono necessari almeno 3 punti dati.")

    # Inizializzazione Variabili
    sumx = 0 # Somma di tutte le x (Σx). Serve per calcolare la media xm
    sumy = 0# Somma di tutte le y (Σy). Serve per calcolare la media ym
    sumxy = 0 # Somma dei prodotti x*y (Σxy). Fondamentale per la covarianza
    sumx2 = 0 # Somma dei quadrati di x (Σx²). Fondamentale per la varianza di x
    st = 0 # (Sum Total) Dispersione totale dei dati rispetto alla MEDIA.
    sr = 0 # (Sum Residuals) Dispersione dei dati rispetto alla RETTA (errore).

    # Primo Ciclo (Accumulo Somme)
    for i in range(n):
        sumx = sumx + x[i]
        sumy = sumy + y[i]
        sumxy = sumxy + (x[i] * y[i])
        sumx2 = sumx2 + (x[i] * x[i])

    # Calcolo Medie e Coefficienti
    xm = sumx / n
    ym = sumy / n

    # Calcolo a1 (Pendenza)
    numerator = (n * sumxy) - (sumx * sumy)
    denominator = (n * sumx2) - (sumx * sumx)

    if denominator == 0:
        raise ValueError("Impossibile calcolare: il denominatore è 0 (tutti gli x sono uguali?)")

    a1 = numerator / denominator

    # Calcolo a0 (Intercetta)
    a0 = ym - (a1 * xm)

    # Secondo Ciclo (Calcolo Errori)
    for i in range(n):
        # st: scarto quadratico totale rispetto alla media
        st = st + (y[i] - ym) ** 2
        # sr: scarto quadratico dei residui (errore della regressione)
        sr = sr + (y[i] - a1 * x[i] - a0) ** 2

    # Calcolo Statistiche Finali
    # syx: Errore standard della stima
    syx = (sr / (n - 2)) ** 0.5

    # r2: Coefficiente di determinazione
    # Se st è 0 (tutti i valori y sono uguali), r2 non è definito matematicamente, gestiamo il caso.
    if st == 0:
        r2 = 1.0  # Adattamento perfetto (linea orizzontale)
    else:
        r2 = (st - sr) / st

    return a1, a0, syx, r2

# GRAFICO
def disegna_grafico(x, y, a1, a0, r2):
    # Calcolo i punti della retta per il grafico
    # Creo una lista di y previsti usando l'equazione y = a1*x + a0
    y_pred = []
    for val in x:
        y_pred.append(a1 * val + a0)

    # Imposto le dimensioni del grafico
    plt.figure(figsize=(10, 6))

    # Disegno i punti reali (Scatter plot)
    plt.scatter(x, y, color='blue', s=80, label='Dati Reali (Osservati)', zorder=2)

    # Disegno la retta di regressione (Plot)
    plt.plot(x, y_pred, color='red', linewidth=2, label=f'Modello: y = {a1:.2f}x + {a0:.2f}', zorder=1)

    # Aggiungo le linee di errore (Residui)
    # Disegna una linea tratteggiata verticale tra il punto reale e la retta
    for i in range(len(x)):
        plt.plot([x[i], x[i]], [y[i], y_pred[i]], color='gray', linestyle='--', alpha=0.5)

    # Etichette e Titolo
    plt.title(f'Regressione Lineare Semplice\n$R^2$ = {r2:.4f} (Adattamento: {r2 * 100:.1f}%)', fontsize=14)
    plt.xlabel('Variabile Indipendente (x)', fontsize=12)
    plt.ylabel('Variabile Dipendente (y)', fontsize=12)

    # Mostro la griglia e la legenda
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=11)

    # Visualizza il grafico finale
    plt.show()