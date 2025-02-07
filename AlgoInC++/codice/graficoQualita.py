import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# Verifica che ci sia un argomento e che sia valido (k)
if len(sys.argv) != 2:
    print("Uso: python3 graficoQualita.py [k]")
    sys.exit(1)

# Estrai il valore di k dalla riga di comando
k_value = int(sys.argv[1])  # Valore di k fissato

# Determina il file da leggere
filename = 'quality_results_k=' + str(k_value) + '.csv'  # Modificato per includere k nel nome del file
output_image = f'grafico_qualita_k={k_value}.png'  # Nome del file di output

# Leggi il file CSV
df = pd.read_csv(filename, sep=';')

# Verifica se ci sono dati per il valore di k specificato
if df.empty:
    print(f"Nessun dato trovato per k={k_value}.")
    sys.exit(1)

# Definisce i colori per ogni algoritmo
colors = {'KMH': 'blue',  'FSS': 'green'}

# Crea il grafico 12x8 pollici
plt.figure(figsize=(12, 8))

# Plotta i dati per ogni algoritmo
for algoritmo in df['Algoritmo'].unique():
    algo_data = df[df['Algoritmo'] == algoritmo]

    # Calcola la media degli errori RMSE per ogni valore di Jaccard
    mean_rmse = algo_data.groupby('Jaccard_reale').apply(
        lambda x: np.sqrt(np.mean((x['Jaccard_reale'] - x['Jaccard_stimata'])**2))
    ).reset_index()
    mean_rmse.columns = ['Jaccard_reale', 'RMSE']

    variance_rmse = algo_data.groupby('Jaccard_reale').apply(
        lambda x: np.var(np.sqrt((x['Jaccard_reale'] - x['Jaccard_stimata'])**2))
    ).reset_index()

    mean_rmse.columns = ['Jaccard_reale', 'RMSE']
    variance_rmse.columns = ['Jaccard_reale', 'Variance']
    # Plotta la curva per l'algoritmo corrente
    plt.plot(mean_rmse['Jaccard_reale'], mean_rmse['RMSE'], label=algoritmo, color=colors[algoritmo])

    # Aggiungi l'area sfumata per la varianza
    plt.fill_between(mean_rmse['Jaccard_reale'],
                     mean_rmse['RMSE'] - np.sqrt(variance_rmse['Variance']),
                     mean_rmse['RMSE'] + np.sqrt(variance_rmse['Variance']),
                     color=colors[algoritmo], alpha=0.2)

# Personalizza il grafico
plt.xlabel('Jaccard Similarity Reale')
plt.ylabel('RMSE')
plt.title(f'Errore RMSE rispetto alla Jaccard Similarity per k={k_value}')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

# Salva il grafico
plt.savefig(output_image)
plt.close()

print(f"Grafico salvato come {output_image}.")