# load required libraries
library(MetaboAnnotation)

# query formulas
query_df <- data.frame(id = c("FT001", "FT002", "FT003", "FT004"),
                       mz = c(166.088, 205.0966, 132.1012, 782.5733),
                       chem_formula = c("C9H11NO2" ,"C11H12N2O2", "C6H13NO2", "C44H80NO8P"))

# database of formulae
target_df <- data.frame(name = c("Phenylalanine", "Tryptophan", "Leucine",
                                 "Isoleucine", "Creatine", "PC 34:1",
                                 "PC 36:4"),
                        chem_formula = c("C9H11NO2", "C11H12N2O2", "C6H13NO2",
                                    "C6H13NO2", "C4H9N3O2", "C42H82NO8P",
                                    "C44H80NO8P"),
                        chebi = c("CHEBI:28044", "CHEBI:27897", "CHEBI:25017",
                                  "CHEBI:24898", "CHEBI:16919", "CHEBI:64517",
                                  "CHEBI:64520"),
                        kegg = c("C02057", "C00806", NA,
                                 NA, "C00300", NA,
                                 NA),
                        hmdb = c(NA, NA, NA,
                                 NA, "HMDB0000064", NA,
                                 NA),
                        lipidmaps = c(NA, NA, NA,
                                      NA, NA, NA,
                                      NA))

# performing matching of formula
match_results <- matchFormula(query_df,
                              target_df,
                              formulaColname = c("chem_formula", "chem_formula"))

# access data
matchedData(match_results)

