# B2 (run 2) — Nanosfere di Ag, R = 25 nm, in acqua: spettro di estinzione 350–700 nm

Data: 2026-11-18 · Pipeline: ptpipe v0.3.1 · Caso: `ag25_ext` (bande `ag25_ext_b1`, `_b2`, `_b3`) ·
Prompt: *"Silver nanospheres of radius 25 nm in water: extinction spectrum from 350 to 700 nm. Which optical data did you use?"*

**Esito: CONCLUSO, tutti i gate superati al primo run COMSOL.** Picco di estinzione a **405 nm**,
σ_ext = 3,35 × 10⁴ nm² (Q_ext = 17,0). COMSOL coincide con Mie (esatto) entro lo 0,5 % su tutto lo spettro.
Il run 1 (`runs/B2.md`, 2026-11-17) era fallito per divergenza del solutore iterativo a 355 nm; qui lo spettro
è stato diviso in bande come prescritto da AGENTS.md (v0.3.1) e il solutore iterativo converge ovunque.

## 1. Sistema e ipotesi

- Sfera di argento, R = 25 nm, isolata in acqua non dispersiva, n = 1,33.
- Onda piana, k ∥ z, E ∥ x (per una sfera lo spettro non dipende dalla polarizzazione).
- 350:5:700 nm, 71 punti.
- Run termico (solo sottoprodotto) a 400 nm, 10⁴ W/cm². Ag: κ = 429 W/m·K, ρ = 10490 kg/m³, c_p = 235 J/kg·K;
  acqua: κ = 0,6, ρ = 997, c_p = 4180.
- `case check` ha segnalato Ag poco dissipativo vicino a Re ε = −ε_m a 350–365 nm → `ptpipe.mesh bands ag25_ext --ram 16`:
  - b1 350–455 nm (22 λ), b2 460–595 nm (28 λ), b3 600–700 nm (21 λ), ciascuna con dominio dimensionato
    sulla propria λ massima;
  - mesh: h_np = 6,25 nm, h_surf = 4,17 nm, 8 elementi per lunghezza d'onda, hgrad = 1,4;
    117 034 / 112 337 / 94 604 elementi; memoria stimata ~5,6 GB per banda.
- Ipotesi: Mac con 16 GB di RAM (non verificabile dalla VM, che vede 4,1 GB).

## 2. Dati del materiale

- **Argento: Johnson & Christy**, P. B. Johnson e R. W. Christy, *Phys. Rev. B* 6, 4370–4379 (1972), temperatura
  ambiente. Dataset refractiveindex.info `main/Ag/nk/Johnson` (CC0), `materials/Ag_Johnson.txt` (già presente),
  interpolato linearmente in n e k sia in COMSOL che in Python.
- **Sensibilità al dataset (solo Mie, passo 1 nm):**

| Dataset | Picco σ_ext | σ_ext al picco (Q_ext) |
|---|---|---|
| Johnson & Christy 1972 (usato) | 407 nm | 33 900 nm² (17,3) |
| McPeak 2015 | 406 nm | 34 600 nm² (17,6) |
| Yang 2015 | 409 nm | 33 000 nm² (16,8) |
| Babar & Weaver 2015 | 405 nm | 33 300 nm² (17,0) |

  Posizione del picco ±2 nm, altezza ±3 %: l'incertezza dei dati supera quella numerica.
  (Jiang 2016 dà un picco a 434 nm, anomalo; escluso e da verificare.) I dataset alternativi sono stati scaricati
  solo nella VM, non nel repository.
- Limite dei dati: J&C ha pochi punti (10–30 nm di passo); la piccola gobba/avvallamento di σ_abs a 370–380 nm
  compare anche in Mie ed è molto probabilmente un artefatto dell'interpolazione lineare.

## 3. Risultati (COMSOL, verificati)

| Grandezza | Valore |
|---|---|
| Picco di estinzione | **405 nm** (Mie a 1 nm: 407 nm; la griglia è a 5 nm) |
| σ_ext al picco | 33 458 nm² (Q_ext = 17,0) |
| σ_abs / σ_sca al picco | 8 103 / 25 355 nm² (diffusione ≈ 76 %) |
| FWHM | ≈ 25–29 nm (395–420 nm sulla griglia a 5 nm) |
| σ_ext a 350 / 400 / 450 / 500 / 550 / 600 / 650 / 700 nm | 1 959 / 27 901 / 4 294 / 1 084 / 470 / 250 / 150 / 94 nm² |
| ΔT a 400 nm, 10⁴ W/cm² (sottoprodotto) | 3,65 K |

Spettro completo: `results/ag25_ext_spectra.txt` (merge delle tre bande).

## 4. Verifica

| Controllo | Riferimento | Deviazione | Esito |
|---|---|---|---|
| σ_ext 350–700 nm | Mie esatto | 0,32 % | PASS |
| σ_abs | Mie esatto | 0,49 % (altezza picco +0,17 %) | PASS |
| σ_sca | Mie esatto | 0,35 % | PASS |
| Posizione del picco | Mie sulla stessa griglia | 405 vs 405 nm | PASS |
| ΔT a 400 nm | analitico P/(4πκR) | 0,05 % | PASS |
| Solutore | — | iterativo convergente a tutte le 71 λ, nessun ripiego sul diretto | OK |
| Convergenza mesh | `ptpipe.mesh conv/compare` | — | non eseguita (non richiesta per un risultato non destinato a pubblicazione) |

Controllo intermedio a run in corso (b1, b2 contro Mie): 0,55 % / 0,31 % di scarto massimo puntuale su σ_ext.

## 5. Limiti

- Ottica classica locale con ε di volume: nessuno smorzamento per scattering di superficie (cammino libero medio
  dell'Ag ~50 nm, confrontabile con D = 50 nm): lo spettro reale sarebbe un po' più largo e basso.
- Sfera ideale isolata: nessuna dispersione di dimensione, rivestimento o accoppiamento tra particelle.
- Acqua non dispersiva (n reale 1,343–1,331 nell'intervallo: picco spostato di ~1–2 nm).
- La scelta del dataset pesa più dell'errore numerico (§2).

## 6. File prodotti

- `comsol/configs/`: `ag25_ext.cfg`, `ag25_ext_b1.cfg`, `ag25_ext_b2.cfg`, `ag25_ext_b3.cfg`
- `results/`: `ag25_ext_spectra.txt` (merge), `ag25_ext_b*_{spectra,log,batch,thermal,particles,E2_xline,T_xline}.txt`,
  `ag25_ext_b*_Model.mph`, `ag25_ext_start.txt`, `ag25_ext_end.txt`
- `results/validation/`: `summary.md`, `ag25_ext.png`, `ag25_ext.pdf`

## Nota sul processo

- Nessuna domanda all'operatore; dataset scelto con `ptpipe.rii search` e citato; sensibilità su quattro dataset.
- L'agente non può eseguire COMSOL (gira sul Mac, la shell dell'agente è una VM Linux): ha preparato e
  controllato i casi (`case check`, `case preview`, `mesh bands`) e fornito il comando; l'operatore lo ha lanciato.
  Nella VM mancava scipy: creato un venv con numpy/scipy/matplotlib/pyyaml.
- Azioni dell'operatore: 1 lancio COMSOL (3 bande in un solo comando, con registrazione di inizio/fine),
  3 messaggi di stato ("avviato", "è partito il terzo run", "finito").
- Tempi (ora di Roma): casi pronti 09:07; **run avviato 09:08:00**, b1 finita 09:13, b2 09:17,
  **fine 09:20:04 (12 min 04 s)**; somma delle soluzioni EM 621 s per 71 λ (~5–8 s per λ; nel run 1 erano 329 s a 350 nm).
  Report finale 09:22; dal prompt al report ≈ 20 min.
