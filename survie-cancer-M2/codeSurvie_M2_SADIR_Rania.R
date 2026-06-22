#__________________________________________________
# Projet Survie - SADIR Rania
#__________________________________________________

# Packages
library(survival)
library(splines)
library(survPen)

# Chargement du jeu de données 
Dir <- "/Users/rania/Downloads/UE Survie 2025-2026-20260116/Projet_2025_2026/"
load(paste0(Dir, "don.RData"))

#__________________________________________________
# Statistiques descriptives
#__________________________________________________

# Description globale de l’âge et du temps de suivi
summary(don$age)
summary(don$fu)

# Répartition du sexe
table(don$sex)
prop.table(table(don$sex)) * 100

# Répartition du statut vital (0 = censure, 1 = décès)
table(don$event)                   
prop.table(table(don$event)) * 100

# Répartition des classes d’âge
table(don$class.age)
prop.table(table(don$class.age)) * 100

# Quartiles de l’âge et du temps de suivi
quantile(don$age, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
quantile(don$fu,  probs = c(0.25, 0.50, 0.75), na.rm = TRUE)

# Comparaison de l’âge selon le sexe
quantile(don$age[don$sex == "F"], probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
quantile(don$age[don$sex == "H"], probs = c(0.25, 0.50, 0.75), na.rm = TRUE)

# Comparaison du temps de suivi selon le sexe
quantile(don$fu[don$sex == "F"], probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
quantile(don$fu[don$sex == "H"], probs = c(0.25, 0.50, 0.75), na.rm = TRUE)

# Statut vital selon le sexe
table(don$event[don$sex == "F"])
prop.table(table(don$event[don$sex == "F"])) * 100

table(don$event[don$sex == "H"])
prop.table(table(don$event[don$sex == "H"])) * 100

# Classes d’âge par sexe (pourcentages par sexe)
table(don$class.age, don$sex)
prop.table(table(don$class.age, don$sex), margin = 2) * 100

# Figures 1 et 2 (une seule figure par fenetre)
par(mfrow = c(1, 1))
hist(don$age, breaks = 30, main = "Distribution de l'age au diagnostic", xlab = "Age (ans)")

par(mfrow = c(1, 1))
hist(don$fu, breaks = 30, main = "Distribution du temps de suivi", xlab = "Temps de suivi (annees)")

#__________________________________________________
# Estimation non parametrique (Kaplan-Meier)
#__________________________________________________

# Estimation de la survie globale par la méthode de Kaplan-Meier
KM_global <- survfit(Surv(fu, event) ~ 1, data = don)

# Résumé de la courbe de survie : médiane et intervalles de confiance
summary(KM_global)
summary(KM_global)$table
summary(KM_global)$table["median"]
summary(KM_global)$table[c("0.95LCL", "0.95UCL")]

# Probabilité de survie à 2 ans et intervalle de confiance
summary(KM_global, times = 2)$surv
summary(KM_global, times = 2)$lower
summary(KM_global, times = 2)$upper

# Courbe de survie globale de Kaplan-Meier avec intervalle de confiance
par(mfrow = c(1, 1))
plot(KM_global, conf.int = TRUE,
     main = "Survie globale (Kaplan-Meier)",
     xlab = "Temps (annees)", ylab = "Probabilite de survie")

# KM selon sexe (meme couleurs que dans la validation Cox)
KM_sex <- survfit(Surv(fu, event) ~ sex, data = don)

# Représentation graphique des courbes de survie selon le sexe
par(mfrow = c(1, 1))
plot(KM_sex, conf.int = TRUE, col = c("forestgreen", "blue"),
     main = "Survie selon le sexe (Kaplan-Meier)",
     xlab = "Temps (annees)", ylab = "Probabilite de survie")
legend("topright", legend = c("Femmes", "Hommes"),
       col = c("forestgreen", "blue"), lty = 1, bty = "n")

# Test du log-rank pour comparer les courbes de survie entre femmes et hommes
survdiff(Surv(fu, event) ~ sex, data = don)

# Estimation des courbes de survie de Kaplan-Meier selon les classes d’âge
KM_ageclass <- survfit(Surv(fu, event) ~ class.age, data = don)

# Représentation graphique des courbes de survie selon les classes d’âge
par(mfrow = c(1, 1))
plot(KM_ageclass, conf.int = FALSE, col = 1:nlevels(don$class.age),
     main = "Survie selon les classes d'age (Kaplan-Meier)",
     xlab = "Temps (annees)", ylab = "Probabilite de survie")
legend("topright", legend = levels(don$class.age),
       col = 1:nlevels(don$class.age), lty = 1, cex = 0.8, bty = "n")


# Test du log-rank pour comparer les courbes de survie entre les classes d’âge
survdiff(Surv(fu, event) ~ class.age, data = don)

#__________________________________________________
# Modelisation semi-parametrique  
#__________________________________________________

###  Ajustement de modèles de Cox avec complexité croissante

# Modèle de Cox 1 : le sexe comme seule covariable
mod1 <- coxph(Surv(fu, event) ~ sex, data = don)
summary(mod1)

# Modèle 2 cox : l’âge (en classes) comme seule covariable
mod2 <- coxph(Surv(fu, event) ~ class.age, data = don)
summary(mod2)

# Modèle 3 de Cox additif : sexe et âge
mod3 <- coxph(Surv(fu, event) ~ sex + class.age, data = don)
summary(mod3)

# Modèle 4 de Cox avec interaction entre le sexe et l’âge
mod4 <- coxph(Surv(fu, event) ~ sex * class.age, data = don)  # équivaut à sex + age + sex:age
summary(mod4)

### # Comparaison des modèles emboîtés par test du rapport de vraisemblance (LRT)

# Apport de l'âge au-delà du sexe : mod1 vs mod3
anova(mod1, mod3, test = "LRT")
# Apport du sexe au-delà de l'âge : mod2 vs mod3
anova(mod2, mod3, test = "LRT")
# Apport de l'interaction sexe × âge : mod3 vs mod4
anova(mod3, mod4, test = "LRT")


# Hypothese des risques proportionnels (Schoenfeld)  
test_PH <- cox.zph(mod3)
test_PH

# Représentation graphique des résidus de Schoenfeld pour l’âge et le sexe
par(mfrow = c(1, 2))
plot(test_PH, var = "class.age", main = "Residus de Schoenfeld (age)",
     xlab = "Temps", ylab = "Coefficient estime")
plot(test_PH, var = "sex", main = "Residus de Schoenfeld (sexe)",
     xlab = "Temps", ylab = "Coefficient estime")

# Analyse des résidus de Martingale pour l’évaluation de la forme fonctionnelle
mart <- residuals(mod3, type = "martingale")


# Résidus de Martingale selon l’âge et le sexe
par(mfrow = c(1, 2))
boxplot(mart ~ class.age, data = don,
        xlab = "Classe d'âge",
        ylab = "Résidus de Martingale",
        main = "Martingale vs classes d'âge")

boxplot(mart ~ sex, data = don,
        xlab = "Sexe",
        ylab = "Résidus de Martingale",
        main = "Martingale vs sexe")

abline(h = 0, lty = 2, col = "red")
par(mfrow = c(1, 1))


# Comparaison des courbes de survie Kaplan-Meier et des prédictions du modèle de Cox selon le sexe
age_med <- median(don$age, na.rm = TRUE)

# on fixe la classe  
class_fix <- levels(don$class.age)[ findInterval(age_med, don$breaks.age, rightmost.closed = TRUE) ]
class_fix <- names(which.max(table(don$class.age)))

pred_F <- survfit(mod3, newdata = data.frame(sex = "F", class.age = class_fix))
pred_H <- survfit(mod3, newdata = data.frame(sex = "H", class.age = class_fix))

par(mfrow = c(1, 1))
plot(KM_sex, col = c("forestgreen", "blue"), conf.int = TRUE,
     main = "Kaplan-Meier vs Cox selon le sexe (classe d'age fixee)",
     xlab = "Temps (annees)", ylab = "Probabilite de survie")
lines(pred_F$time, pred_F$surv, lwd = 2, lty = 2)
lines(pred_H$time, pred_H$surv, lwd = 2, lty = 2)
legend("bottomleft",
       legend = c("KM Femmes", "KM Hommes", "Cox Femmes", "Cox Hommes"),
       col = c("forestgreen", "blue", "black", "black"),
       lty = c(1, 1, 2, 2), bty = "n")

# Comparaison KM vs Cox par classes d'age a sexe fixe (H)
sexe <- "H"

par(mfrow = c(2, 3), mar = c(3, 3, 1.5, 0.5), mgp = c(3, 1, 0))

for (i in 1:nlevels(don$class.age)) {
  
  data.cst <- don[don$sex == sexe & don$class.age == levels(don$class.age)[i], ]
  
  KM <- survfit(Surv(fu, event) ~ 1, data = data.cst)
  
  Cox_pred <- survfit(mod3,
                      newdata = data.frame(sex = sexe,
                                           class.age = levels(don$class.age)[i]))
  
  plot(KM, col = "forestgreen", conf.int = TRUE,
       main = levels(don$class.age)[i],
       xlab = "Temps (annees)", ylab = "Survie")
  
  lines(Cox_pred$time, Cox_pred$surv, lwd = 2)
  
  if (i == nlevels(don$class.age)) {
    legend("bottomleft",
           legend = c("Kaplan-Meier", "Cox"),
           col = c("forestgreen", "black"),
           lty = 1, bty = "n")
  }
}

par(mfrow = c(1, 1))

# Comparaison KM vs Cox par classes d'age a sexe fixe (F)
sexe <- "F"

par(mfrow = c(2, 3), mar = c(3, 3, 1.5, 0.5), mgp = c(3, 1, 0))

for (i in 1:nlevels(don$class.age)) {
  
  data.cst <- don[don$sex == sexe & don$class.age == levels(don$class.age)[i], ]
  
  KM <- survfit(Surv(fu, event) ~ 1, data = data.cst)
  
  Cox_pred <- survfit(mod3,
                      newdata = data.frame(sex = sexe,
                                           class.age = levels(don$class.age)[i]))
  
  plot(KM, col = "forestgreen", conf.int = TRUE,
       main = levels(don$class.age)[i],
       xlab = "Temps (annees)", ylab = "Survie")
  
  lines(Cox_pred$time, Cox_pred$surv, lwd = 2)
  
  if (i == nlevels(don$class.age)) {
    legend("bottomleft",
           legend = c("Kaplan-Meier", "Cox"),
           col = c("forestgreen", "black"),
           lty = 1, bty = "n")
  }
}

par(mfrow = c(1, 1))

#__________________________________________________
#### 4) Modelisation parametrique du taux de mortalite (survPen)
#__________________________________________________

horizon <- 2
ntime <- seq(0, horizon, length.out = 200)

# Modele a taux constant  
mod_const <- survPen(~ 1, t1 = fu, event = event, data = don)
summary(mod_const)

# Modele a taux constant par intervalles (piecewise constant)
temps.pwcst <- seq(0, 2, by = 0.25)   # meme choix que dans les TD (modifiable)
mod_pwcst_global <- survPen(~ pwcst(breaks = temps.pwcst), t1 = fu, event = event, data = don)
summary(mod_pwcst_global)

# Comparaison (taux constant) vs (piecewise constant) VS Kaplan-Meier
# Choix du sexe fixe
sexe <- "H"  

# Bornes pour le modele piecewise constant  
temps.pwcst <- seq(0, 2, by = 0.25)

par(mfrow = c(2, 3), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))

for(i in 1:6){
  
  # Sous-groupe : sexe fixe + classe d'age i
  data.cst <- don[don$sex == sexe & don$class.age == levels(don$class.age)[i], ]
  
  # Kaplan-Meier dans le sous-groupe
  KM <- survfit(Surv(fu, event) ~ 1, data = data.cst)
  
  # Modèle à taux constant (dans le sous-groupe)
  mod_const_i <- survPen(~ 1, t1 = fu, event = event, data = data.cst)
  
  # Modèle à taux constant par intervalles (dans le sous-groupe)
  mod_pwcst_i <- survPen(~ pwcst(breaks = temps.pwcst), t1 = fu, event = event, data = data.cst)
  
  # Survie prédite par les deux modèles
  pred_const <- predict(mod_const_i, data.frame(fu = ntime))$surv
  pred_pwcst <- predict(mod_pwcst_i, data.frame(fu = ntime))$surv
  
  # Graphe
  plot(KM, conf.int = TRUE, col = "forestgreen",
       main = levels(don$class.age)[i],
       xlab = "Temps (annees)", ylab = "Probabilite de survie")
  
  lines(ntime, pred_const, lwd = 2, lty = 2, col = "black")
  lines(ntime, pred_pwcst, lwd = 2, lty = 1, col = "black")
  
  if(i == 6){
    par(xpd = TRUE)  # autorise le dessin en dehors du cadre
legend("topright",
       inset = c(-0.25, 0),
       legend = c("Kaplan-Meier", "Taux constant", "Taux intervalles"),
       col = c("forestgreen", "black", "black"),
       lty = c(1, 2, 1),
       lwd = c(1, 2, 2),
       cex = 0.75,
       bty = "n")
par(xpd = FALSE)
  }
}

par(mfrow = c(1, 1))

# Comparaison (taux constant) vs (piecewise constant) VS Kaplan-Meier
# Choix du sexe fixe femmes
sexe <- "F"  

# Bornes pour le modele piecewise constant  
temps.pwcst <- seq(0, 2, by = 0.25)

par(mfrow = c(2, 3), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))

for(i in 1:6){
  
  # Sous-groupe : sexe fixe + classe d'age i
  data.cst <- don[don$sex == sexe & don$class.age == levels(don$class.age)[i], ]
  
  # Kaplan-Meier dans le sous-groupe
  KM <- survfit(Surv(fu, event) ~ 1, data = data.cst)
  
  # Modèle à taux constant (dans le sous-groupe)
  mod_const_i <- survPen(~ 1, t1 = fu, event = event, data = data.cst)
  
  # Modèle à taux constant par intervalles (dans le sous-groupe)
  mod_pwcst_i <- survPen(~ pwcst(breaks = temps.pwcst), t1 = fu, event = event, data = data.cst)
  
  # Survie prédite par les deux modèles
  pred_const <- predict(mod_const_i, data.frame(fu = ntime))$surv
  pred_pwcst <- predict(mod_pwcst_i, data.frame(fu = ntime))$surv
  
  # Graphe
  plot(KM, conf.int = TRUE, col = "forestgreen",
       main = levels(don$class.age)[i],
       xlab = "Temps (annees)", ylab = "Probabilite de survie")
  
  lines(ntime, pred_const, lwd = 2, lty = 2, col = "black")
  lines(ntime, pred_pwcst, lwd = 2, lty = 1, col = "black")
  
  if(i == 6){
    par(xpd = TRUE)  # autorise le dessin en dehors du cadre
    legend("topright",
           inset = c(-0.25, 0),
           legend = c("Kaplan-Meier", "Taux constant", "Taux intervalles"),
           col = c("forestgreen", "black", "black"),
           lty = c(1, 2, 1),
           lwd = c(1, 2, 2),
           cex = 0.75,
           bty = "n")
    par(xpd = FALSE)
  }
}

par(mfrow = c(1, 1))

#_________________________________________________________

# Modele spline cubique sur le temps (noeud a 1 an) + sexe + age (lineaire)
mod_spline_t <- survPen(
  ~ bs(fu, knots = 1, Boundary.knots = c(0, 2)) + sex + age,
  t1 = fu, event = event, data = don
)
summary(mod_spline_t)

# Modele spline temps (noeud 1 an) + spline age (noeud a l'age median) + sexe
age_med <- median(don$age, na.rm = TRUE)

mod_spline_t_age <- survPen(
  ~ bs(fu,  knots = 1,       Boundary.knots = c(0, 2)) +
    bs(age, knots = age_med, Boundary.knots = c(40, 80)) +
    sex,
  t1 = fu, event = event, data = don
)
summary(mod_spline_t_age)

# Comparaison spline temps vs spline temps+age
## Test du rapport de vraisemblance (LRT)  
LRT_survPen <- abs(2 * (mod_spline_t_age$ll.unpen - mod_spline_t$ll.unpen))
df_test <- abs(mod_spline_t_age$p - mod_spline_t$p)
pval_LRT <- pchisq(LRT_survPen, df = df_test, lower.tail = FALSE)

LRT_survPen
df_test
pval_LRT

## Comparaison des survies predites aux KM dans chaque classe d'age (sexe fixe)
sexe_fix <- "H"   

par(mfrow = c(2, 3), mar = c(3, 3, 1.5, 0.5), mgp = c(3, 1, 0))

for (i in 1:nlevels(don$class.age)) {
  
  data.cst <- don[don$sex == sexe_fix & don$class.age == levels(don$class.age)[i], ]
  age_med_i <- median(data.cst$age, na.rm = TRUE)
  
  KM_i <- survfit(Surv(fu, event) ~ 1, data = data.cst)
  
  pred_surv_t <- predict(mod_spline_t,
                         data.frame(fu = ntime, sex = sexe_fix, age = age_med_i))$surv
  
  pred_surv_t_age <- predict(mod_spline_t_age,
                             data.frame(fu = ntime, sex = sexe_fix, age = age_med_i))$surv
  
  plot(KM_i, conf.int = TRUE, col = "forestgreen",
       main = levels(don$class.age)[i],
       xlab = "Temps (annees)", ylab = "Survie")
  
  lines(ntime, pred_surv_t, lwd = 2)                 # spline temps + age lineaire
  lines(ntime, pred_surv_t_age, lwd = 2, col = "red")# spline temps + spline age
  
  if (i == nlevels(don$class.age)) {
    legend("bottomleft",
           legend = c("Kaplan-Meier", "Spline temps", "Spline temps + spline age"),
           col = c("forestgreen", "black", "red"),
           lty = 1, bty = "n", cex = 0.9)
  }
}

par(mfrow = c(1, 1))

## Comparaison des survies predites aux KM dans chaque classe d'age (sexe fixe)
sexe_fix <- "F"   

par(mfrow = c(2, 3), mar = c(3, 3, 1.5, 0.5), mgp = c(3, 1, 0))

for (i in 1:nlevels(don$class.age)) {
  
  data.cst <- don[don$sex == sexe_fix & don$class.age == levels(don$class.age)[i], ]
  age_med_i <- median(data.cst$age, na.rm = TRUE)
  
  KM_i <- survfit(Surv(fu, event) ~ 1, data = data.cst)
  
  pred_surv_t <- predict(mod_spline_t,
                         data.frame(fu = ntime, sex = sexe_fix, age = age_med_i))$surv
  
  pred_surv_t_age <- predict(mod_spline_t_age,
                             data.frame(fu = ntime, sex = sexe_fix, age = age_med_i))$surv
  
  plot(KM_i, conf.int = TRUE, col = "forestgreen",
       main = levels(don$class.age)[i],
       xlab = "Temps (annees)", ylab = "Survie")
  
  lines(ntime, pred_surv_t, lwd = 2)                 # spline temps + age lineaire
  lines(ntime, pred_surv_t_age, lwd = 2, col = "red")# spline temps + spline age
  
  if (i == nlevels(don$class.age)) {
    legend("bottomleft",
           legend = c("Kaplan-Meier", "Spline temps", "Spline temps + spline age"),
           col = c("forestgreen", "black", "red"),
           lty = 1, bty = "n", cex = 0.9)
  }
}

par(mfrow = c(1, 1))


## Comparaison des hazards (splines vs piecewise constant) par classes d'age (sexe fixe)

sexe_fix <- "H"   # changer en "F" si besoin
par(mfrow = c(2, 3), mar = c(3, 3, 1.5, 0.5), mgp = c(3, 1, 0))

for (i in 1:nlevels(don$class.age)) {
  
  data.cst <- don[don$sex == sexe_fix & don$class.age == levels(don$class.age)[i], ]
  age_med_i <- median(data.cst$age, na.rm = TRUE)
  
  # Modele piecewise constant (dans le sous-groupe)
  mod_pwcst_i <- survPen(~ pwcst(breaks = temps.pwcst),
                         t1 = fu, event = event, data = data.cst)
  
  # Hazards predits par les 2 modeles spline (sexe fixe, age fixe)
  haz_t <- predict(mod_spline_t,
                   data.frame(fu = ntime, sex = sexe_fix, age = age_med_i))$haz
  
  haz_t_age <- predict(mod_spline_t_age,
                       data.frame(fu = ntime, sex = sexe_fix, age = age_med_i))$haz
  
  # Hazard piecewise constant (par intervalle)
  log_taux <- summary(mod_pwcst_i)$coefficients[, "Estimate"]
  taux <- exp(log_taux)
  
  # Ylim adapte (evite de couper des segments)
  ylim_max <- max(c(haz_t, haz_t_age, taux), na.rm = TRUE)
  ylim_max <- 1.05 * ylim_max  # petite marge
  
  # Tracés
  plot(ntime, haz_t, type = "l", lwd = 2, ylim = c(0, ylim_max),
       main = levels(don$class.age)[i],
       xlab = "Temps (annees)", ylab = "Taux instantane (hazard)")
  
  lines(ntime, haz_t_age, lwd = 2, col = "red")
  
  # Segments piecewise constant
  for (j in 1:(length(temps.pwcst) - 1)) {
    segments(x0 = temps.pwcst[j],   y0 = taux[j],
             x1 = temps.pwcst[j+1], y1 = taux[j],
             col = "forestgreen", lwd = 2)
  }
  
  if (i == nlevels(don$class.age)) {
    legend("topleft",
           legend = c("Spline temps", "Spline temps + spline age", "Taux constant par intervalles"),
           col = c("black", "red", "forestgreen"),
           lty = 1, lwd = 2, bty = "n", cex = 0.8)
  }
}

par(mfrow = c(1, 1))
