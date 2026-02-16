# ================================================================================
# STANDARDISATION DE LA RICHESSE SPÉCIFIQUE PAR COVERAGE (iNEXT)
# ================================================================================
# Ce script standardise la richesse entre échantillons en utilisant le 
# "sample coverage" comme base de comparaison avec le package iNEXT
# ================================================================================

# 1. CHARGEMENT DES PACKAGES ====================================================
library(iNEXT)      # Pour les analyses de raréfaction/extrapolation
library(tidyverse)

# 2. PRÉPARATION DES DONNÉES ====================================================

# Option A: Données sous forme de matrice (espèces × sites)
# ---------------------------------------------------------
# Créer ou importer vos données d'abondance
donnees_abondance <- data.frame(
  Site1 = c(45, 12, 89, 3, 23, 7, 15, 2, 34, 8),
  Site2 = c(23, 34, 45, 7, 11, 2, 8, 19, 21, 5),
  Site3 = c(67, 8, 12, 45, 34, 19, 25, 11, 42, 15)
)
rownames(donnees_abondance) <- paste0("Espece_", 1:10)

# Convertir en liste (format requis pour iNEXT)
donnees_list <- list(
  Site1 = donnees_abondance$Site1,
  Site2 = donnees_abondance$Site2,
  Site3 = donnees_abondance$Site3
)

# Option B: Importer depuis un fichier CSV
# -----------------------------------------
# donnees_csv <- read.csv("mes_donnees.csv", row.names = 1)
# donnees_list <- as.list(donnees_csv)


# 3. VÉRIFICATION DES DONNÉES ===================================================

# Résumé des échantillons
cat("\n=== RÉSUMÉ DES ÉCHANTILLONS ===\n")
for(i in 1:length(donnees_list)) {
  cat(paste0("\n", names(donnees_list)[i], ":\n"))
  cat(paste0("  - Nombre d'individus: ", sum(donnees_list[[i]]), "\n"))
  cat(paste0("  - Richesse observée: ", sum(donnees_list[[i]] > 0), "\n"))
}


# 4. ANALYSE iNEXT AVEC COVERAGE ================================================

# Calculer la couverture d'échantillonnage et les estimations de diversité
# datatype = "abundance" car nous avons des données d'abondance
# endpoint spécifie jusqu'où extrapoler (par défaut: 2× la taille de l'échantillon)

cat("\n=== ANALYSE iNEXT ===\n")
resultats_inext <- iNEXT(
  x = donnees_list,
  q = 0,                    # q = 0 pour la richesse spécifique
  datatype = "abundance",   # Type de données
  endpoint = NULL,          # NULL = extrapolation jusqu'à 2× la taille
  knots = 40,              # Nombre de points pour les courbes
  conf = 0.95              # Intervalle de confiance à 95%
)

# Afficher le résumé
print(resultats_inext)


# 5. CALCUL DE LA COUVERTURE DE CHAQUE ÉCHANTILLON ==============================

cat("\n=== COUVERTURE D'ÉCHANTILLONNAGE ===\n")

# Calculer la couverture observée pour chaque site
coverages <- sapply(donnees_list, function(x) {
  n <- sum(x)
  f1 <- sum(x == 1)  # Nombre de singletons
  coverage <- 1 - (f1 / n)
  return(coverage)
})

print(data.frame(
  Site = names(coverages),
  Coverage = round(coverages, 4),
  Pourcentage = paste0(round(coverages * 100, 2), "%")
))


# 6. STANDARDISATION PAR COVERAGE COMMUN ========================================

# Déterminer la couverture minimale (base de comparaison)
coverage_min <- min(coverages)
cat(paste0("\nCouverture minimale: ", round(coverage_min, 4), 
           " (", round(coverage_min * 100, 2), "%)\n"))

# Estimer la richesse standardisée pour chaque site à la couverture minimale
cat("\n=== STANDARDISATION PAR COVERAGE ===\n")

richesse_standardisee <- estimateD(
  x = donnees_list,
  datatype = "abundance",
  base = "coverage",           # Standardisation par coverage
  level = coverage_min,        # Niveau de coverage cible
  conf = 0.95
)

# Afficher les résultats
print(richesse_standardisee)


# 7. COMPARAISON AVANT/APRÈS STANDARDISATION ====================================

# Richesse observée (non standardisée)
richesse_observee <- sapply(donnees_list, function(x) sum(x > 0))

# Créer un tableau comparatif
comparaison <- data.frame(
  Site = names(donnees_list),
  N_individus = sapply(donnees_list, sum),
  Richesse_observee = richesse_observee,
  Coverage_observe = round(coverages, 4),
  Richesse_standardisee = round(richesse_standardisee$qD, 2),
  IC_inf = round(richesse_standardisee$qD.LCL, 2),
  IC_sup = round(richesse_standardisee$qD.UCL, 2)
)

cat("\n=== TABLEAU COMPARATIF ===\n")
print(comparaison)

# Exporter les résultats
write.csv(comparaison, "resultats_standardisation_coverage.csv", row.names = FALSE)


# 8. VISUALISATIONS =============================================================

# 8.1 Courbes de raréfaction/extrapolation
# -----------------------------------------
plot_courbes <- ggiNEXT(resultats_inext, type = 1, color.var = "Assemblage") +
  theme_bw() +
  labs(
    title = "Courbes de raréfaction/extrapolation",
    x = "Nombre d'individus",
    y = "Richesse spécifique (q = 0)"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(plot_courbes)
ggsave("courbes_rarefaction.png", width = 10, height = 6, dpi = 300)


# 8.2 Couverture d'échantillonnage
# ---------------------------------
plot_coverage <- ggiNEXT(resultats_inext, type = 2, color.var = "Assemblage") +
  theme_bw() +
  labs(
    title = "Couverture d'échantillonnage",
    x = "Nombre d'individus",
    y = "Couverture (%)"
  ) +
  geom_hline(yintercept = coverage_min, linetype = "dashed", color = "red") +
  annotate("text", x = Inf, y = coverage_min, 
           label = paste0("Coverage min = ", round(coverage_min, 3)),
           hjust = 1.1, vjust = -0.5, color = "red") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(plot_coverage)
ggsave("courbes_coverage.png", width = 10, height = 6, dpi = 300)


# 8.3 Richesse en fonction de la couverture
# ------------------------------------------
plot_richesse_coverage <- ggiNEXT(resultats_inext, type = 3, 
                                   color.var = "Assemblage") +
  theme_bw() +
  labs(
    title = "Richesse standardisée par coverage",
    x = "Couverture d'échantillonnage",
    y = "Richesse spécifique (q = 0)"
  ) +
  geom_vline(xintercept = coverage_min, linetype = "dashed", color = "red") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(plot_richesse_coverage)
ggsave("richesse_vs_coverage.png", width = 10, height = 6, dpi = 300)


# 8.4 Graphique comparatif personnalisé
# --------------------------------------
ggplot(comparaison, aes(x = Site)) +
  geom_point(aes(y = Richesse_observee, color = "Observée"), size = 4) +
  geom_point(aes(y = Richesse_standardisee, color = "Standardisée"), size = 4) +
  geom_errorbar(aes(ymin = IC_inf, ymax = IC_sup, color = "Standardisée"), 
                width = 0.2, linewidth = 1) +
  scale_color_manual(values = c("Observée" = "#E74C3C", 
                                 "Standardisée" = "#3498DB")) +
  theme_bw() +
  labs(
    title = "Comparaison richesse observée vs standardisée (coverage)",
    subtitle = paste0("Coverage de référence: ", round(coverage_min, 4)),
    x = "Site",
    y = "Richesse spécifique",
    color = "Type"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave("comparaison_richesse.png", width = 10, height = 6, dpi = 300)


# 9. ANALYSE STATISTIQUE (OPTIONNEL) ===========================================

# Test de chevauchement des intervalles de confiance
cat("\n=== COMPARAISONS ENTRE SITES ===\n")
for(i in 1:(nrow(comparaison)-1)) {
  for(j in (i+1):nrow(comparaison)) {
    site1 <- comparaison$Site[i]
    site2 <- comparaison$Site[j]
    
    ic1 <- c(comparaison$IC_inf[i], comparaison$IC_sup[i])
    ic2 <- c(comparaison$IC_inf[j], comparaison$IC_sup[j])
    
    overlap <- (ic1[2] >= ic2[1]) & (ic2[2] >= ic1[1])
    
    cat(paste0(site1, " vs ", site2, ": "))
    if(overlap) {
      cat("Pas de différence significative (IC se chevauchent)\n")
    } else {
      cat("Différence potentiellement significative (IC ne se chevauchent pas)\n")
    }
  }
}


# 10. NOTES ET INTERPRÉTATION ==================================================

cat("\n")
cat("================================================================================\n")
cat("NOTES IMPORTANTES:\n")
cat("================================================================================\n")
cat("1. Sample Coverage:\n")
cat("   - Représente la proportion d'espèces détectées dans l'échantillon\n")
cat("   - Varie de 0 à 1 (0% à 100%)\n")
cat("   - Plus la couverture est élevée, plus l'inventaire est complet\n\n")
cat("2. Standardisation:\n")
cat("   - La richesse est estimée au niveau du coverage minimum commun\n")
cat("   - Permet de comparer équitablement des échantillons de tailles différentes\n")
cat("   - Les intervalles de confiance indiquent l'incertitude des estimations\n\n")
cat("3. Interprétation:\n")
cat("   - Si les IC se chevauchent: pas de différence significative\n")
cat("   - Si les IC ne se chevauchent pas: différence potentiellement significative\n")
cat("   - Regarder aussi les graphiques pour une vision d'ensemble\n")
cat("================================================================================\n")


# 11. FONCTION RÉUTILISABLE ====================================================

standardiser_richesse_coverage <- function(data_list, confidence = 0.95) {
  #' Standardise la richesse spécifique entre échantillons par coverage
  #'
  #' @param data_list Liste nommée contenant les abondances par site
  #' @param confidence Niveau de confiance pour les IC (défaut: 0.95)
  #' @return Liste contenant les résultats et les graphiques
  
  # Analyse iNEXT
  inext_results <- iNEXT(data_list, q = 0, datatype = "abundance", conf = confidence)
  
  # Calcul des coverages
  coverages <- sapply(data_list, function(x) {
    1 - (sum(x == 1) / sum(x))
  })
  
  # Standardisation
  coverage_min <- min(coverages)
  richness_std <- estimateD(data_list, datatype = "abundance", 
                            base = "coverage", level = coverage_min, 
                            conf = confidence)
  
  # Résultats
  results <- list(
    inext = inext_results,
    coverages = coverages,
    coverage_min = coverage_min,
    richness_standardized = richness_std,
    comparison = data.frame(
      Site = names(data_list),
      N = sapply(data_list, sum),
      Obs_Richness = sapply(data_list, function(x) sum(x > 0)),
      Coverage = coverages,
      Std_Richness = richness_std$qD,
      CI_lower = richness_std$qD.LCL,
      CI_upper = richness_std$qD.UCL
    )
  )
  
  return(results)
}

# Exemple d'utilisation de la fonction:
# resultats <- standardiser_richesse_coverage(donnees_list)
# print(resultats$comparison)

cat("\n✓ Script terminé avec succès!\n")
cat("Fichiers générés:\n")
cat("  - resultats_standardisation_coverage.csv\n")
cat("  - courbes_rarefaction.png\n")
cat("  - courbes_coverage.png\n")
cat("  - richesse_vs_coverage.png\n")
cat("  - comparaison_richesse.png\n")
