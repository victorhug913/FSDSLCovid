library(readr)
library(dplyr)
library(rstatix)
library(corrplot)

# Read the data
importdb <- read.csv2(file = "Dados_Brutos.csv", stringsAsFactors = TRUE, fileEncoding = "WINDOWS-1252")

# Filter the confirmed cases with laboratory confirmation and closed notification status
confirmados <- importdb %>%
  filter(Classificacao == "Confirmados",
         CriterioConfirmacao == "Laboratorial",
         StatusNotificacao == "Encerrado") %>%
  select(DataDiagnostico, DataEncerramento, DataObito, Evolucao, CriterioConfirmacao, Municipio, 
         FaixaEtaria, Sexo, RacaCor, Escolaridade, Gestante, Febre, DificuldadeRespiratoria, Tosse, 
         Coriza, DorGarganta, Diarreia, Cefaleia, ComorbidadePulmao, ComorbidadeCardio, ComorbidadeRenal, 
         ComorbidadeDiabetes, ComorbidadeTabagismo, ComorbidadeObesidade, FicouInternado, ViagemBrasil, 
         ViagemInternacional, ProfissionalSaude, PossuiDeficiencia, MoradorDeRua, ResultadoSorologia_IGG, 
         ResultadoSorologia, ResultadoRT_PCR, ResultadoTesteRapido, TipoTesteRapido)

# Filter the closed cases with outcome information
ndesfecho <- confirmados %>%
  filter(!Evolucao %in% c("-", "Ignorado", "Óbito por outras causas"))

# Recode the outcome variable
confirmados$Obito <- ifelse(confirmados$Evolucao == "Óbito pelo COVID-19", "Sim", "Não")

# Remove cases with missing values in selected variables
vars_to_clean <- c("Escolaridade", "RacaCor", "ComorbidadeCardio", "ComorbidadeDiabetes", "ComorbidadePulmao", 
                   "ComorbidadeRenal", "ComorbidadeTabagismo", "ComorbidadeObesidade")
for (var in vars_to_clean) {
  confirmados <- confirmados %>%
    filter(!var %in% c("-", "Ignorado", "Não se aplica"))
}

# Show the number of cases
cat(paste("Number of cases:", nrow(confirmados), "\n"))

# Load the required packages
library(mfx)
library(factoextra)
library(ggplot2)
library(psych)
library(car)

# Tabela de variáveis categóricas

IcProporcao <- function(a, b){
  
  db <<- db %>% mutate_all(as.character)
  db <<- db %>% mutate_all(as.factor)
  
  tabela22 <- table(db$desfecho, {{a}})
  testeq <- chisq.test(tabela22, correct = FALSE)
  print(testeq)
  tabela <- db[complete.cases({{a}}),] %>% group_by({{a}}[complete.cases({{a}})]) %>% 
    summarise("n" = length(db$Obito),
              "Totaln" = length(desfecho),
              "Total%" = length(desfecho)/length(db$desfecho[complete.cases({{a}})])*100,
    )
  names(tabela)[1] <- "Valor"
  DescCat <<- rbind(DescCat, cbind(Variável=paste({{b}}, ", n (%)"),tabela))
  
  tabela2 <- db[complete.cases({{a}}),] %>% group_by({{a}}[complete.cases({{a}})]) %>% 
    summarise("Óbiton" = length(desfecho[desfecho=="Sim"]),
              "Óbito%" = length(desfecho[desfecho=="Sim"])/(length(db$desfecho[(db$desfecho=="Sim") & (complete.cases({{a}}))]))*100,
              "Curan" = length(desfecho[desfecho=="Não"]),
              "Cura%" = length(desfecho[desfecho=="Não"])/(length(db$desfecho[(db$desfecho=="Não") & complete.cases({{a}})]))*100,
              "Valor-p"=testeq["p.value"],
              "V de cramer" = cramer_v(db$desfecho, {{a}})
    )
  names(tabela2)[1] <- "Valor"
  tabela2 <- cbind(tabela2, resajus)
  AssocCat <<- rbind(AssocCat, cbind(Variável=paste({{b}}, ", n (%)"),tabela2))

  nome <- b
  
  regressao <- glm(db$desfecho~a, family = binomial(link="logit"))
  tabela <- cbind("Desfecho"={{b}}, exp(cbind(OR=coef(regressao), confint(regressao, level = 0.99))))
  print(tabela[-1,])
  RegBru <<- rbind(RegBru, tabela)
}

RegBru <- data.frame()
confirmados$desfecho <- confirmados$Obito;
confirmados$desfecho <- as.factor(confirmados$desfecho)
db <- confirmados
DescCat <- data.frame()
PAnalise <- data.frame()
AssocCat <- data.frame()

db$FaixaEtaria <- as.character(db$FaixaEtaria)
db$FaixaEtaria[db$FaixaEtaria=="0 a 4 anos"] <- "0 a 19 anos"
db$FaixaEtaria[db$FaixaEtaria=="05 a 9 anos"] <- "0 a 19 anos"
db$FaixaEtaria[db$FaixaEtaria=="10 a 19 anos"] <- "0 a 19 anos"
db$FaixaEtaria[db$FaixaEtaria=="20 a 29 anos"] <- "20 a 39 anos"
db$FaixaEtaria[db$FaixaEtaria=="30 a 39 anos"] <- "20 a 39 anos"
db$FaixaEtaria[db$FaixaEtaria=="40 a 49 anos"] <- "40 a 59 anos"
db$FaixaEtaria[db$FaixaEtaria=="50 a 59 anos"] <- "40 a 59 anos"
db$FaixaEtaria[db$FaixaEtaria=="60 a 69 anos"] <- "60 anos ou mais"
db$FaixaEtaria[db$FaixaEtaria=="70 a 79 anos"] <- "60 anos ou mais"
db$FaixaEtaria[db$FaixaEtaria=="80 a 89 anos"] <- "80 anos ou mais"
db$FaixaEtaria[db$FaixaEtaria=="90 anos ou mais"] <- "80 anos ou mais"
db$FaixaEtaria[db$FaixaEtaria=="-"] <- NA
db$FaixaEtaria <- as.factor(db$FaixaEtaria)
#IcProporcao(db$FaixaEtaria, "Faixa Etária")
summary(db$FaixaEtaria)

db1 <- db[db$FaixaEtaria=="60 anos ou mais"|db$FaixaEtaria=="80 anos ou mais",]
db <- db1

db$Sexo <- as.character(db$Sexo)
db$Sexo[db$Sexo=="I"] <- NA
db$Sexo <- as.factor(db$Sexo)
IcProporcao(db$Sexo, "Sexo")

db$RacaCor <- as.character(db$RacaCor)
db$RacaCor[db$RacaCor=="Branca"] <- "ABranca"
db$RacaCor[db$RacaCor=="Ignorado"] <- NA
db$RacaCor <- as.factor(db$RacaCor)
summary(db$RacaCor)
IcProporcao(db$RacaCor, "Raça")

summary(as.factor(db$Escolaridade))
db$Escolaridade <- as.character(db$Escolaridade)
db$Escolaridade[db$Escolaridade=="Ignorado"] <- NA
db$Escolaridade[db$Escolaridade=="Não se aplica"] <- NA
db$Escolaridade[db$Escolaridade=="Analfabeto"] <- "5) Analfabeto"
db$Escolaridade[db$Escolaridade=="1ª a 4ª série incompleta do EF (antigo primário ou 1º grau)"] <- "4) Ensino Fundamental Incompleto"
db$Escolaridade[db$Escolaridade=="4ª série completa do EF (antigo primário ou 1º grau)"] <- "4) Ensino Fundamental Incompleto"
db$Escolaridade[db$Escolaridade=="5ª à 8ª série incompleta do EF (antigo ginásio ou 1º grau)"] <- "4) Ensino Fundamental Incompleto"
db$Escolaridade[db$Escolaridade=="Ensino fundamental completo (antigo ginásio ou 1º grau) "] <- "3) Ensino Fundamental Completo"
db$Escolaridade[db$Escolaridade=="Ensino médio incompleto (antigo colegial ou 2º grau )"] <- "3) Ensino Fundamental Completo"
db$Escolaridade[db$Escolaridade=="Ensino médio completo (antigo colegial ou 2º grau ) "] <- "2) Ensino Médio Completo"
db$Escolaridade[db$Escolaridade=="Educação superior incompleta "] <- "2) Ensino Médio Completo"
db$Escolaridade[db$Escolaridade=="Educação superior completa"] <- "1) Ensino Superior Completo"
summary(as.factor(db$Escolaridade))
db$Escolaridade <- as.factor(db$Escolaridade)
IcProporcao(db$Escolaridade, "Escolaridade")

db$ComorbidadeCardio[db$ComorbidadeCardio=="-"] <- NA
IcProporcao(db$ComorbidadeCardio, "Cardíaca")

db$ComorbidadePulmao[db$ComorbidadePulmao=="-"] <- NA
IcProporcao(db$ComorbidadePulmao, "Pulmonar")

db$ComorbidadeRenal[db$ComorbidadeRenal=="-"] <- NA
IcProporcao(db$ComorbidadeRenal, "Renal")

db$ComorbidadeDiabetes[db$ComorbidadeDiabetes=="-"] <- NA
IcProporcao(db$ComorbidadeDiabetes, "Diabetes")

db$ComorbidadeTabagismo[db$ComorbidadeTabagismo=="-"] <- NA
IcProporcao(db$ComorbidadeTabagismo, "Tabagismo")

db$ComorbidadeObesidade[db$ComorbidadeObesidade=="-"] <- NA
IcProporcao(db$ComorbidadeObesidade, "Obesidade")

# Ajustes
ajuste1 <- glm(db$desfecho~db$Sexo+db$RacaCor+db$Escolaridade, family = binomial(link="logit"))
ajuste1 <- exp(cbind(OR=coef(ajuste1), confint(ajuste1, level = 0.99)))

ajuste2 <- glm(db$desfecho~db$Sexo+db$RacaCor+db$Escolaridade+db$ComorbidadeCardio+db$ComorbidadeDiabetes+db$ComorbidadeObesidade+db$ComorbidadePulmao+db$ComorbidadeRenal+db$ComorbidadeTabagismo, family = binomial(link="logit"))
ajuste2 <- exp(cbind(OR=coef(ajuste2), confint(ajuste2, level = 0.99)))


require('data.table')
fwrite(DescCat, file ="Saida/DescCat_Maiores_60.csv.csv")
RegBru$Variavel <- rownames(RegBru)
fwrite(RegBru, file ="Saida/RegBru_Maiores_60.csv.csv")
fwrite(AssocCat, file ="Saida/AssocCat_Maiores_60.csv.csv")
write.csv(ajuste1, file="Saida/ajuste1_Maiores_60.csv.csv")
write.csv(ajuste2, file="Saida/ajuste2_Maiores_60.csv.csv")