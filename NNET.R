

reticulate::use_condaenv("selene", conda = "/home/patrycja/minocona3/condabin/conda",required = TRUE)
reticulate::py_config()
#tensorflow::install_tensorflow(envname = "python310")
#tf$constant("Good")
#reticulate::py_list_packages()
#reticulate::py_install(c("pandas", "scikit-learn", "matplotlib"), envname = "python310")
#K0=keras::backend()
####

library("DeepGenomeScan")
library("caret")### for ML calling functions and performance estimation
library("keras") ### for DL
library("tensorflow")
library("RSNNS")
library("NeuralNetTools")
library(dplyr)

#library("caretEnsemble")
#library("kerasR")
  
  
#### nnet
gen = read.table('merged_freq_sample.txt', header = T)
env = read.csv('land_rda.csv')

gen$Pop = sub("pops/", "", gen$Pop)
gen = filter(gen, Pop != "NOM16") %>% arrange(Pop) %>% select(-Pop)
env = filter(env, Landscape != "NOM16") %>% arrange(Landscape) %>% select(-Landscape)

para = colnames(env)

econtrol <- caret::trainControl( ## 5-fold CV, repeat 5 times
    method = "adaptive_cv",
    number = 5,
    ## repeated ten times
    repeats = 5,
    adaptive = list(min = 5, alpha = 0.05,method = "gls", complete = TRUE),search = "random")

model_nnet_crop = DeepGenomeScan(genotype=gen, env=env[,1],
                               method="nnet",
                               metric = "MAE",## "Accuracy", "RMSE"
                               preProcess=c("scale"),
                               tuneLength = 10, ### search 10 combinations of parameters
                               # verbose=1 is reporting the progress,o is sclience
                               trControl = econtrol,maxit=100,MaxNWts=94581)

nnet_imp_crop = varImp(model_nnet_crop$finalModel)
nnet_oldenimp_crop = NeuralNetTools::olden(model_nnet_crop$finalModel,bar_plot=FALSE)
