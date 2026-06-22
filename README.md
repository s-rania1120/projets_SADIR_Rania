# Projets SADIR Rania

Biostatisticienne - Master 2 B2 AI-EPI, Université Claude Bernard Lyon 1  
Stage de recherche : IUCPQ, Université Laval, Québec (Équipe Arsenault)  
Contact : s.rania1120@gmail.com

---

## Vue d'ensemble

Ce dossier regroupe l'ensemble des travaux académiques et de recherche réalisés entre 2025 et 2026, couvrant la modélisation prédictive, l'analyse de survie, l'inférence statistique et l'analyse de données longitudinales.

| Projet | Contexte | Méthodes clés | Langage |
|--------|----------|---------------|---------|
| [PRED-CAD - Risque coronarien (All of Us)](#1-pred-cad--prédiction-du-risque-coronarien-all-of-us) | Stage M2, IUCPQ | DSM, MICE, OMOP, PRS | Python, R, SQL |
| [MOBILS - Activité physique hospitalière](#2-mobils--activité-physique-et-mobilité-active-hôpital-lyon-sud) | Stage M1, CRNH | Longitudinal, PAEE, régression | R, SAS |
| [Analyse de survie - Cancer](#3-analyse-de-survie--mortalité-par-cancer) | M2 UE Survie | KM, Cox, survPen, splines | R |
| [Censure et inférence - Ebola](#4-censure-et-inférence--charge-virale-résiduelle-ebola) | M2 UE Inférence | MV censurée, simulation MC | R |

---

## 1. PRED-CAD - Prédiction du risque coronarien (All of Us)

**Contexte :** Stage M2, IUCPQ - Université Laval, Québec (Fév. - Juil. 2026). Superviseurs : Louis-Jacques Ruel (doctorant), Pr Benoît Arsenault.  
**Repo principal :** [github.com/s-rania1120/DSM-Pred-AoU](https://github.com/s-rania1120/DSM-Pred-AoU)

### Objectif

Reproduire et étendre le modèle Deep Survival Machines (DSM, Chen et al. 2025, *Nature Medicine*) pour prédire le risque incident de maladie coronarienne (CAD) sur la cohorte All of Us (AoU), en intégrant variables cliniques et score polygénique de risque (PRS).

### Cohorte analytique

- **N = 85 276** participants (30 640 hommes, 54 636 femmes)
- **2 887 événements CAD** (3,4%) définis via 57 codes SNOMED
- Suivi moyen : 3,7 ans | Date de censure : 2023-10-01 (CDR v8)
- Âge cible : 18-79 ans

### Pipeline

```
Extraction OMOP (BigQuery, CDR C2024Q3R8)
    └── Feature engineering (HDL, LDL, PAS, tabac, statines, antihypertenseurs...)
        └── Gestion des valeurs manquantes par MICE (N = 241 600)
            └── Entraînement DSM stratifié par sexe (80/20)
                └── Optimisation Optuna (k, layer_size, lr)
                    └── Validation interne (AoU 20%)
                        └── Validation externe (CARTaGENE)
```

### Résultats

| Population | C-index DSM | FRS | PCE | SCORE2 | PREVENT |
|------------|-------------|-----|-----|--------|---------|
| Hommes | **0,719** | 0,62 | 0,61 | 0,63 | 0,64 |
| Femmes | **0,698** | 0,62 | 0,61 | 0,63 | 0,64 |

Validation externe CARTaGENE : C-index **0,68-0,70** (hommes), **0,76-0,78** (femmes).

**Hyperparamètres optimaux (Optuna) :**

| Sexe | k | layer_size | lr | Distribution |
|------|---|------------|----|--------------|
| Hommes | 8 | 184 | 2,01×10⁻⁵ | Weibull |
| Femmes | 4 | 80 | 4,81×10⁻⁵ | Weibull |

### Phase 2 (en cours)

Intégration d'un PRS CAD (PRSmix, 2,8M variants, PLINK2) pour tester l'apport incrémental de l'information génomique au-delà des variables cliniques (comparaison C-index, LRT, NRI/IDI).

### Technologies

`Python` `R` `SQL` `GCP BigQuery` `Hail` `PLINK2` `Optuna` `MICE` `Verily Workbench` `Git`

---

## 2. MOBILS - Activité physique et mobilité active, Hôpital Lyon Sud

**Contexte :** Stage M1, CRNH Rhône-Alpes - Hôpital Lyon Sud (Avr. - Juin 2025).  
Encadrantes : Louise Seconda, Muriel Rabilloud.

### Objectif

Evaluer l'impact de l'extension de la ligne B du métro lyonnais (mise en service oct. 2023) sur les comportements de mobilité et la dépense énergétique liée à l'activité physique (PAEE) du personnel hospitalier du site Lyon Sud des HCL.

### Design

Cohorte longitudinale prospective à **3 temps de mesure** :

| Temps | Période | Description |
|-------|---------|-------------|
| T0 | Avant ouverture | Mesure de référence |
| T6 | 6 mois après | Premier suivi |
| T12 | 12 mois après | Suivi principal |

- **502 agents** inclus (volontariat, toutes catégories professionnelles)
- **Analyse principale : N = 129** répondants communs T0-T12
- Questionnaires validés : STAQ (activité physique et sédentarité), EQ-5D-5L (qualité de vie)

### Méthodes

**Calcul de la PAEE**

Estimation individuelle de la dépense énergétique liée à l'activité physique par la méthode MET (Metabolic Equivalent of Task), adaptée au contexte hospitalier :

```
PAEE_activité = [durée × MET × poids × (3,5 × 20,3 × 60 / 1000)]
              − [RMR_predicted × (durée / 24)]
              + RT(h/j)
```

**Variable synthétique d'évolution des transports en commun (3 groupes)**

- Groupe 1 : adoption des TC après ouverture (non-utilisateurs → utilisateurs)
- Groupe 2 : maintien d'un usage déjà présent à T0
- Groupe 3 : maintien d'une absence d'usage à T0 et T12

**Tests statistiques**

- Normalité : Shapiro-Wilk
- Evolution T0→T12 : Student apparié (test principal) + Wilcoxon apparié (robustesse)
- Comparaison inter-groupes : Kruskal-Wallis (3 groupes) ; Wilcoxon non apparié (groupes 1 vs 3)

**Régression linéaire multiple progressive (4 modèles emboîtés)**

| Modèle | Variables ajoutées | Critère de comparaison |
|--------|--------------------|------------------------|
| M0 | Evolution TC (groupes 1/2/3) | - |
| M1 | + Variables de mobilité | AIC + ANOVA séquentielle |
| M2 | + Variables sociodémographiques | AIC + ANOVA séquentielle |
| M3 | + Contexte de vie et de travail | AIC + ANOVA séquentielle |

Contrôle de la multicolinéarité : VIF (seuil < 5). Diagnostics sur modèle final : Cook, normalité et homoscédasticité des résidus.

### Résultats principaux

- Adoption des transports en commun associée à une **augmentation significative de la PAEE transport**
- Pas d'effet global significatif sur la PAEE totale (résultat cohérent avec la littérature sur les interventions de mobilité à court terme)
- Score EQ-5D-5L stable entre T0 et T12 (moyenne 0,923 vs 0,920 ; médiane 0,929 aux deux temps)
- Population majoritairement féminine (78,3%), âge moyen centré sur 35-50 ans

### Limite principale

Attrition importante : 502 participants à T0 → 129 analysés (25,7% de rétention). Risque de biais de sélection entre répondants et non-répondants à T12 non quantifiable dans le cadre du stage.

### Fichiers

- `MemoireStage.docx` — Mémoire complet (40 pages, incluant plan d'analyse statistique, résultats et discussion)

### Technologies

`R` `SAS` `Excel`

---

## 3. Analyse de survie — Mortalité par cancer

**Contexte :** Projet M2, UE M2_5 Analyse de survie — Université Lyon 1 (2025-2026).  
Enseignants : Pr Delphine Maucort-Boulch, Dr Mathieu Fauvernier.  
**Auteurs :** SADIR Rania, FAKHAR Salma

### Objectif

Modéliser la relation entre l'âge au diagnostic, le sexe, le temps depuis le diagnostic et la mortalité sur une cohorte simulée de patients atteints de cancer (N = 20 000 ; 86,1% de décès ; suivi max 2 ans).

### Méthodes

**Estimation non paramétrique**
- Kaplan-Meier global, stratifié par sexe et par 6 classes d'âge
- Log-rank : sexe (χ² = 931, p < 2×10⁻¹⁶), classes d'âge (χ² = 394, p < 2×10⁻¹⁶)

**Modèles de Cox (stratégie progressive)**

| Modèle | Variables | HR sexe (H vs F) | LRT vs modèle sans |
|--------|-----------|------------------|-------------------|
| M1 | Sexe seul | 1,59 [1,54-1,64] | — |
| M2 | Age (classes) | — | — |
| M3 | Sexe + âge (additif) | 1,60 [1,55-1,65] | χ² = 940, p < 2,2×10⁻¹⁶ |
| M4 | Sexe × âge (interaction) | — | χ² = 8,96, p = 0,111 (NS) |

Modèle retenu : **M3** (additif)

**Validation du modèle de Cox**
- Schoenfeld : pas de violation pour les classes d'âge (p = 0,39) ; violation pour le sexe (p < 2×10⁻¹⁶) → effet du sexe interprété comme effet moyen sur la période
- Martingale : distributions centrées sur 0, pas de défaut d'ajustement majeur par groupe
- Comparaison graphique KM vs Cox prédit : concordance satisfaisante dans toutes les classes

**Modélisation paramétrique du taux de mortalité**

| Modèle | Caractéristique | Ajustement vs KM |
|--------|----------------|------------------|
| Taux constant | h(t) = exp(β₀) = 0,837/an | Insuffisant |
| Taux constant par intervalles (Δ = 0,25 an) | Taux de 0,57 → 1,52 sur 2 ans | Meilleur |
| Spline cubique temps (nœud 1 an) + âge linéaire + sexe | Termes spline tous significatifs | Bon |
| Spline temps + spline âge (nœud médian) | LRT χ² = 3,47, p = 0,33 vs modèle précédent | Equivalent, non retenu |

### Fichiers

- `codeSurvie_M2_SADIR_Rania_et_FAKHAR_Salma.R` — Code R (532 lignes)
- `ProjetSurvie_M2_SADIR_Rania_et_FAKHAR_Salma.pdf` — Rapport complet

### Technologies

`R` `survival` `survPen` `splines` `ggplot2`

---

## 4. Censure et inférence — Charge virale résiduelle Ebola

**Contexte :** Projet M2, UE Inférence fréquentiste et bayésienne — Université Lyon 1 (2025-2026).  
**Auteurs :** SADIR Rania, FAKHAR Salma

### Objectif

Modéliser la relation entre le nombre de cycles PCR (Ct) et le délai depuis la sortie du traitement chez 409 survivants d'Ebola, en prenant correctement en compte la censure à droite imposée par la limite de détection PCR (seuil C = 40 cycles).

### Méthodes

**Modèle 1 — Régression linéaire sans censure (référence)**

```
ct_i = β₀ + β₁ × time_i + ε_i,   ε_i ~ N(0, σ²)
```

Résultats : β₀ = 32,00 ; β₁ = 0,0187 ; σ = 2,17 ; R² = 0,598  
Biais : surestimation β₀, sous-estimation β₁ et σ ; hétéroscédasticité sur les résidus.

**Modèle 2 — Vraisemblance censurée**

Pour les valeurs censurées (ct = 40) : contribution = P(ct* ≥ C) = 1 − Φ((C − μᵢ)/σ).  
Maximisation numérique via `nlm()`.

Résultats : β₀ = 30,93 ; β₁ = 0,0268 (+43% vs modèle naïf) ; σ = 2,674

**Simulation Monte Carlo (1 000 réplications × 3 seuils)**

| Seuil C | Méthode | β₁ estimé moyen | Biais relatif |
|---------|---------|-----------------|---------------|
| 40 | Naïf | 0,01882 | −30% |
| 40 | Censuré | 0,02683 | ~0% |
| 38 | Naïf | 0,01399 | −48% |
| 38 | Censuré | 0,02685 | ~0% |
| 35 | Naïf | 0,00676 | −75% |
| 35 | Censuré | 0,02700 | ~0% |

Le biais du modèle naïf augmente avec la proportion de censure ; l'estimateur censuré reste non biaisé à tous les seuils, mais sa variance croît (perte d'information).

### Fichiers

- `CodeR_IFB_Final.R` — Code R (144 lignes)
- `Projet_IFB.docx` — Rapport d'analyse

### Technologies

`R` `nlm()` `ggplot2`

---

## Compétences transversales

| Domaine | Détail |
|---------|--------|
| Langages | R (avancé), Python (intermédiaire), SAS (intermédiaire), SQL (intermédiaire), Bash |
| Survie et modélisation | Cox, Kaplan-Meier, survPen, splines cubiques, DSM, MV censurée |
| Données longitudinales | Plan d'analyse 3 temps, tests appariés, régression progressive, AIC, VIF |
| Données génomiques | PLINK2, Hail, PRS computation, OMOP CDM, SNOMED |
| Infrastructure | GCP BigQuery, Dataproc, Verily Workbench, Git/GitHub |
| Données manquantes | MICE (241 600 observations), imputation multiple |
| Optimisation | Optuna, simulation Monte Carlo (1 000 réplications) |

---

## Organisation du dossier

```
projets_SADIR_Rania/
├── README.md                               (ce fichier)
├── PRED-CAD/                               (voir repo dédié : DSM-Pred-AoU)
├── MOBILS-M1/
│   └── MemoireStage.docx
├── survie-cancer-M2/
│   ├── codeSurvie_M2_SADIR_Rania_et_FAKHAR_Salma.R
│   └── ProjetSurvie_M2_SADIR_Rania_et_FAKHAR_Salma.pdf
└── censure-ebola-IFB/
    ├── CodeR_IFB_Final.R
    └── Projet_IFB.docx
```

---

*Dernière mise à jour : juin 2026*
