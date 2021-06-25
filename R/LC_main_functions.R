#LC main analyses functions

#Define the functions





#this is for normalizing some dataset
function_DESeq2_null_model=function(raw_counts_table,main_clinical_df,treshold=5){
  colnames(raw_counts_table)[1] = "ID"
  raw_counts_table$ID = as.character(raw_counts_table$ID)
  #select design variables and remove NAs
  dataset_df = main_clinical_df
  
  training_df = dataset_df
  
  rownames(training_df) = training_df$sample #assign row names
  
  training_ids=training_df$sample
  
  #in the null model function, i made filtering less relaxed
  #21/03/2019 I removed filtering and added BDg as a factor
  filtered_table_training=raw_counts_table %>% group_by(ID) %>% dplyr::select(one_of(training_ids)) %>% gather(sample,read.counts,one_of(training_ids)) %>% 
    mutate(quant=quantile(read.counts,0.20)) %>% 
    filter(quant >=treshold)  %>% dplyr::select(ID,sample,read.counts) %>% spread(sample,read.counts) 
  
  
  filtered_counts_training_df=data.frame(filtered_table_training %>% ungroup() %>% dplyr::select(-ID),row.names = filtered_table_training$ID)
  
  filtered_counts_training_df=filtered_counts_training_df[,training_ids]
  
  formula_=as.formula(rlang::parse_expr(paste("~",paste("BDg",collapse = " + "))))
  
  dds_training=DESeq(DESeqDataSetFromMatrix(countData = filtered_counts_training_df, colData = training_df, design =~1))
  
  #vst_training=varianceStabilizingTransformation(dds_training,blind = T)
  
  #rpm_training=function_rpm_normalization(filtered_counts_training_df)
  
  #return(list(vst_training=assay(vst_training)))
  
  return(list(dds_training=dds_training))
}


function_pathway_analyses=function(tables_df,p_value_treshold=0.05){
  #tables_df = tables_list %>% map(~.$results_training) %>% bind_rows(.id="Type")
  
  gene_symbols=tables_df %>% filter(Type == "mRNA",padj <= p_value_treshold) %>% .$ID
  mirna_ids=tables_df %>% dplyr::mutate(ID=str_replace_all(ID,"_","-"))  %>% filter(Type == "miRNA",padj <= p_value_treshold) %>% .$ID
  isomir_ids=tables_df %>% filter(Type == "isomiR",padj <= p_value_treshold) %>% .$ID
  
  converted_isomir_ids=isomir.reference.table.v22 %>% filter(sequence %in% isomir_ids) %>% .$MiRBase.ID
  
  mirna_targets=mirdb_predictions %>% filter(ID %in% c(mirna_ids,converted_isomir_ids)) %>% .$target
  
  mrna_gene_ids=entrez_gene_symbols_table %>% filter(alias_symbol %in% gene_symbols) %>% .$gene_id
  
  mirna_target_ids=entrez_gene_ids_table %>% filter(accession %in% mirna_targets) %>% .$gene_id
  
  print(mirna_target_ids)
  kegga(unique(c(mrna_gene_ids,mirna_target_ids)),species="Hs") %>% mutate(p.adjust=p.adjust(P.DE,"fdr")) %>% arrange(p.adjust)
  #goana(unique(c(mrna_gene_ids,mirna_target_ids)),species="Hs") %>% mutate(p.adjust=p.adjust(P.DE,"fdr")) %>% arrange(p.adjust)
  
  #enrichDO(unique(c(mrna_gene_ids,mirna_target_ids)),maxGSSize = 10000) %>% as.data.frame()
}





#######################
###Modelling functions#
#######################
#######################

function_roc_plot = function(lc_match_with_vst,design_variables,time_=c("(7,8.4]"),stage_ = c("Localized"),cancertype_ = c("NSCLC","Control"),bdg_=c("Grp4","Grp5"),frac=0.6) {
  
  model_data_df = lc_match_with_vst %>% bind_rows() %>% 
    dplyr::select(-groups) %>% distinct() %>% 
    filter(time_group %in% time_,ID %in% design_variables,stage_group %in% stage_) %>% 
    dplyr::select(sample,cancertype,condition,vst,ID,BDg) %>% 
    distinct() %>% filter(cancertype %in% cancertype_) %>%
    filter(!sample %in% problematic_samples) %>%
    filter(!BDg %in% bdg_) %>% dplyr::select(-BDg) %>% spread(ID,vst)
  
  set.seed(255)
  model_training_cases_df = model_data_df %>% filter(condition == "LC") %>% sample_frac(frac)
  set.seed(255)
  model_training_df = model_data_df %>% filter(condition == "C") %>% sample_n((model_training_cases_df %>% nrow()) + 3) %>% bind_rows(model_training_cases_df)
  model_test_cases_df = model_data_df %>% filter(condition == "LC") %>% anti_join(model_training_cases_df)
  set.seed(255)
  model_test_df = model_data_df %>% filter(condition == "C") %>% anti_join(model_training_df %>% filter(condition == "C")) %>% sample_n((model_test_cases_df %>% nrow())) %>% bind_rows(model_test_cases_df) %>% distinct()
  
  
  
  formula_=as.formula(rlang::parse_expr(paste("condition ~ ",paste(glue::backtick(design_variables),collapse = " + "))))
  
  glm_model_=glm(condition ~  .  ,family=binomial, data = model_training_df %>% dplyr::select(-sample,-cancertype))
  
  
  
  predict_=predict(glm_model_,model_test_df %>% dplyr::select(-sample,-cancertype) ,type=c("response"))
  
  prediction_df = model_test_df %>% dplyr::select(-sample,-cancertype) %>% mutate(prob=predict_,color="prediction")
  model_df = model_training_df %>% dplyr::select(-sample,-cancertype) %>% mutate(prob=glm_model_$fitted.values,color="model")
  
  rocfit_=roc(condition ~  prob ,data =prediction_df)
  auc_value=round(as.numeric(auc(rocfit_)),3)
  
  
  
  
  ggthemr("fresh")
  plot_=ggplot() + 
    geom_roc(data=prediction_df,aes(d=condition,m=prob,color=color),labels = F,n.cuts = 0,linejoin = "round",pointsize = 0) +
    geom_roc(data=model_df,aes(d=condition,m=prob,color = color),labels = F,n.cuts = 0,linejoin = "round",pointsize = 0) +  
    scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1)) + ylab("Sensitivity")  +  geom_label(aes(x=Inf,y=-Inf,label=paste0("AUC = ",auc_value),hjust=1,vjust=0),size=6,inherit.aes = FALSE) + ggtitle(label = paste0(design_variables)) + theme(legend.title = element_blank())
  
  return(list(plot=plot_,training_df=model_training_df,test_df=model_test_df,model=glm_model_))
  
}


function_auc = function(model_test_df,glm_model_) {
  
  
  predict_=predict(glm_model_,model_test_df %>% dplyr::select(-sample,-cancertype) ,type=c("response"))
  
  prediction_df = model_test_df %>% dplyr::select(-sample,-cancertype) %>% mutate(prob=predict_,color="prediction")
  
  rocfit_=roc(condition ~  prob ,data =prediction_df)
  auc_value=round(as.numeric(auc(rocfit_)),3)
  
  return(auc_value)
  
  
  
}




#######################
#######################
###Obsolete functions##
#######################


















function_data_split=function(x,ratio=c(.6,.2,.2),seed=255) {
  
  x_=unique(x)
  set.seed(seed)
  training_=x_[caTools::sample.split(x_,ratio[1]/(sum(ratio)))]
  y_=x_[!x_ %in% training_]
  
  if(ratio[2] != 0) {
    
    set.seed(seed)
    validate_=y_[caTools::sample.split(y_,ratio[2]/(sum(ratio[2:3])))]
    testing_=y_[!y_ %in% validate_]
    
    x=as.character(x)
    
    x[x %in% training_] <- "training"
    x[x %in% validate_] <- "validate"
    x[x %in% testing_] <- "testing"
    
    
    
  } else  {
    
    testing_=y_
    
    x=as.character(x)
    
    x[x %in% training_] <- "training"
    x[x %in% testing_] <- "testing"
    
    
  }
  
  
  return(x)
  
}



function_roc_plot_for_application=function(grid_result=grid_all_runs_sglfast_parsed %>% group_by(t1,w,cancertype) %>% mutate(auc_median=median(auc)) %>% ungroup() %>% arrange(desc(auc_median)),slice=2) {
  
  grid_result %>% slice(slice) -> df
  
  auc=as.vector(df$result[[1]]$rocfit$auc)
  df <- tibble(rocfit=list(df$result[[1]]$rocfit)) %>% hoist(rocfit,Specificity="specificities",Sensitivity="sensitivities") %>% unnest(c(Specificity,Sensitivity)) 
  
  
  ggplot(df,aes(x=Specificity,y=Sensitivity)) + geom_path() + scale_x_reverse() + annotate("text",x=0.15,y=0.05,label=paste("AUC:",round(auc,digits = 2)),size=4)
  
}


function_evaluate_model=function(grid_result=grid_all_runs_sglfast_parsed %>% group_by(t1,t2,cancertype,type) %>% mutate(auc_median=median(auc)) %>% ungroup() %>% arrange(desc(auc_median)),slice=2,model="sglfast") {
  
  grid_result %>% slice(slice) -> df
  
  if(model== "sglfast") {
    
    df$result[[1]]$test %>% select(sample:age_continuous) %>% bind_cols(df$result[[1]]$predict) %>% dplyr::rename(probability=LC) -> df_prob
    features=NULL
  } else {
    
    probs <- df$result[[1]]$predict$response[[sglOptim::best_model(df$result[[1]]$model)]]["LC",]
    df$result[[1]]$test %>% select(sample:age_continuous) %>% mutate(probability=probs) -> df_prob
    
    features=coef(df$result[[1]]$model_fit,sglOptim::best_model(df$result[[1]]$model))["LC",]

  }


  plot1 <- df_prob %>% mutate(timetodiagnose=factor(timetodiagnose,levels=c("[0,1]", "(1,2]", "(2,3]", "(3,4]","(4,5]", "(5,6]",
"(6,7]", "(7,8]", "(8,9]" ,"(9,10]","C"))) %>% ggplot(aes(y=probability,shape=metastase,x=timetodiagnose,color=BDg)) +
    geom_point(size=4) + theme_light() + xlab("") + facet_wrap(BDg~.)
  
  
  return(list(plot=plot1,df=df_prob,slice=df,features=features))
}


function_confidence_interval=function(x) {
  
  n=length(x)
  error=qt(0.975,df=n-1)*sd(x)/sqrt(n)

  
  return(list(l=mean(x) - error,r=mean(x) + error,m=mean(x)))
  
  
}

function_standart_error=function(x) {
  
  n=length(x)
  se=sqrt(var(x)/n)
  
  
  return(list(l=mean(x) - se,r=mean(x) + se,m=mean(x)))
  
  
}


