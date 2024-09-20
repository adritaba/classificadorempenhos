suppressWarnings(library('tidyverse'))
suppressWarnings(library('tidytext'))
suppressWarnings(library('glmnet'))
suppressWarnings(library('broom'))
suppressWarnings(library('vroom'))
suppressWarnings(library('wordcloud'))
suppressWarnings(library('wordcloud2'))
suppressWarnings(library('yardstick'))
suppressWarnings(library('doParallel'))
suppressWarnings(library('hms'))
suppressWarnings(library('beepr'))
registerDoParallel(4)

#Funções auxiliares
totidy <- function(ptext){
  return(
    tidytext <- data.frame(text=strsplit(tolower(ptext), "\\s+")[[1]]) %>% mutate(document=1)
  )
}

classify <- function(pintercept, pcoefs, ptidytext){
  return(
    prediction <- ptidytext %>%
      inner_join(coefs, by = c("text" = "term")) %>%
      group_by(document) %>%
      summarize(score = sum(estimate)) %>%
      mutate(probability = plogis(pintercept + score))
  )
}

# Guarda diretório base do projeto/script
diretorio <- base::getwd()
subdiretorios <- c(
  '2022' = paste0(diretorio,'/empenho.2022/'), 
  '2023' = paste0(diretorio,'/empenho.2023/'))

# Objeto principal que armazenará os .CSVs
empenhos2022 <- NULL
empenhos2023 <- NULL

# Utiliza os N primeiros subdiretórios para análise (critério pode ser diferente)
arquivos <- list.files(path = paste0(diretorio,'/empenho.2022/'),
                       pattern = '*empenho.empenho.csv$',
                       full.names = TRUE,
                       recursive = TRUE)

columns <- c("tdt","id","pseed","psample","wfreq","wmin","auc","cm.precision","cm.recall","cm.accuracy","cm.f1_score","texec","fam")
stats <- data.frame(matrix(ncol=length(columns)))
colnames(stats) <- columns
stats <- stats  %>% drop_na() %>% mutate(fam=as.character(fam),texec=as.character(texec))
hist <- stats 

# stats <- stats %>% drop_na() %>% add_row(id=1,pseed=575,psample=0.01,wfreq=1000,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=2,pseed=575,psample=0.02,wfreq=1000,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=3,pseed=575,psample=0.03,wfreq=1000,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=4,pseed=575,psample=0.04,wfreq=1000,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=5,pseed=575,psample=0.05,wfreq=1000,wmin=2,fam="logit") #
# 
# stats <- stats %>% add_row(id=6,pseed=575,psample=0.01,wfreq=25,wmin=2,fam="logit") #
# stats <- stats %>% add_row(id=7,pseed=575,psample=0.01,wfreq=50,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=8,pseed=575,psample=0.01,wfreq=100,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=9,pseed=575,psample=0.01,wfreq=200,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=10,pseed=575,psample=0.01,wfreq=400,wmin=2,fam="logit")
# 
# stats <- stats %>% add_row(id=11,pseed=575,psample=0.01,wfreq=25,wmin=2,fam="logit") #
# stats <- stats %>% add_row(id=12,pseed=575,psample=0.01,wfreq=25,wmin=3,fam="logit")
# stats <- stats %>% add_row(id=13,pseed=575,psample=0.01,wfreq=25,wmin=4,fam="logit")
# stats <- stats %>% add_row(id=14,pseed=575,psample=0.01,wfreq=25,wmin=5,fam="logit")
# stats <- stats %>% add_row(id=15,pseed=575,psample=0.01,wfreq=25,wmin=6,fam="logit")
# 
# sp <- sample(1:100000,5, replace=F)
# stats <- stats %>% add_row(id=16,pseed=575,psample=0.01,wfreq=25,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=17,pseed=75086,psample=0.01,wfreq=25,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=18,pseed=27156,psample=0.01,wfreq=25,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=19,pseed=94747,psample=0.01,wfreq=25,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=20,pseed=69619,psample=0.01,wfreq=25,wmin=2,fam="logit")
# stats <- stats %>% add_row(id=21,pseed=sp[5],psample=0.01,wfreq=25,wmin=2, fam="logit") #pseed=13020, auc=0.9801486, cm=0.9428005

#stats <- stats %>% drop_na() %>% add_row(id=99,pseed=69619,psample=0.05,wfreq=25,wmin=2,fam="logit")
stats <- stats %>% drop_na() %>% add_row(id=99,pseed=69619,psample=0.03,wfreq=25,wmin=2,fam="logit")

for (row in 1:nrow(stats)) {
  
  print("Iniciando análise...") 
  start.time <- Sys.time()
  stats[row, 'tdt'] <- format(start.time)
  
  #Parâmetros
  print(stats[row,])
  id <- stats[row, "id"]
  pseed <- stats[row, "pseed"]
  psample <- stats[row, "psample"]
  wfreq <- stats[row, "wfreq"]
  wmin <- stats[row, "wmin"]
  fam <- stats[row, "fam"]
  
  set.seed(pseed)
  empenhos2022 <- vroom(sample(arquivos,size=psample*length(arquivos)))
  
  # Seleciona as colunas de interesse para a análise
  empenhos2022 <- empenhos2022 %>% dplyr::select('seq_empenho','cod_municipio','seq_orgao','num_anoexercicio'
                                                 ,'dsc_tipo_empenho','num_empenho','dsc_empenho'
                                                 ,'dsc_funcao','dsc_subfuncao')
  
  # Realiza a quebra das colunas funcao/subfuncao separadas por traço
  empenhos2022 <- empenhos2022 %>% 
    tidyr::separate(col=dsc_funcao,into=c('cod_funcao','dsc_funcao'),sep=' - ') %>%
    dplyr::mutate(cod_funcao=as.numeric(cod_funcao))
  
  empenhos2022 <- empenhos2022 %>% 
    tidyr::separate(col=dsc_subfuncao,into=c('cod_subfuncao','dsc_subfuncao'),sep=' - ') %>%
    dplyr::mutate(cod_funcao=as.numeric(cod_funcao)) %>%
    dplyr::mutate(dsc_empenho=tolower(dsc_empenho))
  
  empenhos2022 <- empenhos2022 %>%
    dplyr::mutate(cod_funcao=replace(cod_funcao, cod_funcao != 10, 99), dsc_funcao=replace(dsc_funcao, dsc_funcao != 'SAÚDE', 'Z-OUTRA')) %>%
    dplyr::select('seq_empenho', 'dsc_empenho', 'dsc_funcao')
  
  # Classifica cada empenho em um documento específico (poderia usar o seq_empenho)
  empenhos2022 <- tibble(empenhos2022) %>% 
    mutate(document = row_number())
  
  # Muda para formato tidy e elimina palavras (arbitrário, parâmetro = 10 é original) - mudar para ser uma função da amostra
  tidy_empenhos2022 <- empenhos2022 %>%
    tidytext::unnest_tokens(word, dsc_empenho) %>%
    dplyr::group_by(word) %>%
    dplyr::filter(n() > wfreq & nchar(word) > wmin) %>%
    dplyr::ungroup()
  
  # Verifica a frequência de palavras
  frequencia <- tidy_empenhos2022 %>%
    dplyr::count(dsc_funcao, word, sort = TRUE) %>%
    dplyr::anti_join(tidytext::get_stopwords(language="pt")) %>%
    dplyr::filter(is.na(as.numeric(word)) & nchar(word) > wmin) %>%
    dplyr::group_by(dsc_funcao)
  
  frequencia %>%
    dplyr::top_n(25) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(aes(reorder_within(word, n, dsc_funcao), n, fill=dsc_funcao)) +
    geom_col(alpha = 0.8, show.legend = FALSE) +
    geom_bar(stat="identity", fill="grey40") +
    theme(panel.background = element_rect(fill = "grey70",
                                          colour = "black",
                                          size = 0.5, 
                                          linetype = "solid")) +
    scale_x_reordered() + 
    coord_flip() +
    facet_wrap(~dsc_funcao, scales="free") +
    scale_y_continuous(expand = c(0,0)) +
    labs (
      x = NULL, y = "Contagem de palavras",
      title = "Palavras mais frequentes após remoção de stop words",
      subtitle = "Algumas palavras tem ranking parecido entre as diferentes funções"
    )
  
  # Nuvem de palavras
  wordcloud2(data=frequencia[frequencia$dsc_funcao=="SAÚDE",][c('word','n')], 
             size=1.0,
             color="grey",
             shape='circle')
  
  # Divisão em base de treino e teste
  empenhos2022_split <- tidy_empenhos2022 %>%
    select(document) %>%
    rsample::initial_split()
  empenhos2022_train <- rsample::training(empenhos2022_split)
  empenhos2022_test <- rsample::testing(empenhos2022_split)
  
  # Gera matriz esparsa de distribuição (treino)
  empenhos2022_sparse <- tidy_empenhos2022 %>%
    count(document, word) %>%
    inner_join(empenhos2022_train) %>%
    cast_sparse(document, word, n)
  
  word_rownames <- as.integer(rownames(empenhos2022_sparse))
  empenhos2022_join <- tibble(document = word_rownames) %>%
    left_join(empenhos2022 %>% select(document,dsc_funcao))
  
  # Geração do modelo de regressão linear, regularização LASSO (alpha=1), RIDGE (alpha=0)
  print("Gerando modelo...")
  is_saude <- empenhos2022_join$dsc_funcao == "SAÚDE"
  
  model <- cv.glmnet(empenhos2022_sparse, is_saude, family="binomial", parallel=TRUE, trace.it=1, nfolds = 8, alpha = 1)
  plot(model)
  plot(model$glmnet.fit)
  print("Modelo gerado!")
  
  # Apresentação dos coeficientes/variáveis com mais relevância
  coefs <- model$glmnet.fit %>%
    tidy() %>%
    filter(lambda == model$lambda.1se)
  
  coefs %>%
    group_by(estimate > 0) %>%
    top_n(15, abs(estimate)) %>%
    ungroup() %>%
    ggplot(aes(fct_reorder(term,estimate), estimate, fill = estimate > 0)) +
    geom_col(alpha = 0.8, show.legend = FALSE) +
    geom_bar(stat="identity", fill="grey40") +
    theme(panel.background = element_rect(fill = "grey70",
                                          colour = "black",
                                          size = 0.5, 
                                          linetype = "solid")) +
    coord_flip() +
    labs(
      x = NULL,
      title = "Coeficientes que aumentam/diminuem a probabilidade de ser da classe Saúde"
    ) +
    ylab("estimativa")
  
  #Extrai o intercept/constante
  intercept <- coefs %>%
    filter(term == "(Intercept)") %>%
    pull(estimate)
  
  #Gera as classificações a partir da junção dos empenhos (formato tidy) com coeficentes e seus pesos
  classif <- tidy_empenhos2022 %>%
    inner_join(empenhos2022_test) %>%
    inner_join(coefs, by = c("word" = "term")) %>%
    group_by(document) %>%
    summarize(score = sum(estimate)) %>%
    mutate(probability = plogis(intercept + score))
  
  #Faz a junção para comparar a base de teste predita com os dados reais
  classif_real <- classif %>%
    left_join(empenhos2022 %>%
                select(dsc_funcao, document), by = "document") %>%
    mutate(dsc_funcao = as.factor(dsc_funcao))
  
  #Gera curva ROC - indica o quão bom é o classificador binário
  curva_roc <- classif_real %>%
    roc_curve(dsc_funcao, probability) %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(size = 1.0, color = "grey40") +
    geom_abline(lty = 2, alpha = 0.5, size = 1.0, color = "slategray4") +
    theme(panel.background = element_rect(fill = "grey70",
                                          colour = "black",
                                          size = 0.5, 
                                          linetype = "solid")) +
    labs(
      title = "Curva ROC",
      subtitle = "Predição de classificação de empenhos em Saúde",
    ) +
    xlab("1-especificidade") + 
    ylab("sensibilidade")
  
  curva_roc
  
  rauc <- classif_real %>%
    roc_auc(dsc_funcao, probability)
  
  #Acurácia do ROC - área sob a curva (quanto mais próxima de 1, melhor)
  stats[row, 'auc'] <- rauc$.estimate[1]
  
  #Matriz confusão
  matriz_confusao <- classif_real %>%
    mutate(
      prediction = case_when(
        probability > 0.4 ~ "SAÚDE",
        TRUE ~ "Z-OUTRA"
      ),
      prediction = as.factor(prediction)
    ) %>%
    conf_mat(dsc_funcao, prediction)
  
  autoplot(matriz_confusao,type="heatmap")
  
  #Matricas da matriz confusão
  tp <- matriz_confusao$table[1,1]
  fp <- matriz_confusao$table[1,2]
  fn <- matriz_confusao$table[2,1]
  tn <- matriz_confusao$table[2,2]
  
  precisao <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  mc_acuracia <- (tp +tn) / sum(matriz_confusao$table)
  mc_f1_score <- 2* (precisao * recall) / (precisao + recall)
  
  #Precisão de matriz confusão
  stats[row, 'cm.precision'] <- precisao
  stats[row, 'cm.recall'] <- recall
  stats[row, 'cm.accuracy'] <- mc_acuracia
  stats[row, 'cm.f1_score'] <- mc_f1_score
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  stats[row, 'texec'] <- as.character(as.hms(time.taken))
  
  print("Fim da análise atual!")
  
}
beep(5)

hist <- hist %>% rbind(stats)
#hist <- hist %>% drop_na()
#hist <- hist[-c(21),]
#hist[27,'id'] <- 27

#Teste de uma instância/empenho
print(classify(pintercept=intercept
               ,pcoefs=coefs
               ,ptidytext=totidy("REFERENTE AO PAGAMENTO DE 13 (TREZE) DIAS TRABALHADOS DO MES DE JANEIRO2022.")
))

#Teste de uma instância/empenho
print(classify(pintercept=intercept
               ,pcoefs=coefs
               ,ptidytext=totidy("VALOR QUE SE EMPENHA REFERENTE A AQUISIÇÃO DE LEITOR DE CARTÃO EXTERNO E CRAZY BOLSA/CÂMERA PARA SECRETARIA DE SAÚDE")
))

#hist_saude <- hist