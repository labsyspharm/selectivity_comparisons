library(dplyr)
library(ggplot2)
library(reshape2)
`%nin%`<-Negate(`%in%`)

#####################################################################################################################################################################T
# set directories  ------------
#####################################################################################################################################################################T

home_dir<-"MyDir"

#####################################################################################################################################################################T
# make score_mat  ------------
#####################################################################################################################################################################T

#diff_dists = list(1,150,300,450,600,750,900,1050,1200,1350,1500,1650,1800,1950,1999)
diff_dists = c(seq(1,11,1))
num_runs = length(diff_dists)
score_mat = matrix(0, nrow = 4, ncol = length(diff_dists)) #Comparing all four scoring methodologies
score_mat = as.data.frame(score_mat)
row.names(score_mat) = list("Entropy","Gini","High Affinity On-Target","Low Affinity On-Target")
colnames(score_mat) = diff_dists

#####################################################################################################################################################################T
# specificity calculations  ------------
#####################################################################################################################################################################T
#run_num<-1
for (run_num in 1:num_runs){
	#Establish Kd profile for drug 
	num_prot_genes = 2000
	binding_prof_Kd = matrix(0, nrow = 1, ncol = num_prot_genes)%>%
	  as.data.frame()
	#binding_prof_Kd = as.data.frame(binding_prof_Kd)
	selectivity_prof = list(10*10^-9, 10*10^-6)
	num_high_aff = diff_dists[run_num]
	binding_prof_Kd[1,1:num_high_aff[[1]]] = selectivity_prof[1]
	binding_prof_Kd[1,(num_high_aff[[1]]+1):num_prot_genes] = selectivity_prof[2]
	
	#Method 1: Calculation of Entropy Value
	binding_prof_Ka = 1/binding_prof_Kd
	Ka_total = rowSums(binding_prof_Ka)
	binding_prof_Ka = (binding_prof_Ka/Ka_total)*log(binding_prof_Ka/Ka_total)
	Ka_adj_total = rowSums(binding_prof_Ka)
	selectivity_entropy = Ka_adj_total*-1
	score_mat[1,run_num] = round(selectivity_entropy,2)
	
	#Method 2: Calculation of Gini Coefficient
	#How To Take Activity Values and Convert To Kd (Or is it even necessary?)
	binding_prof_Ka = 1/binding_prof_Kd
	binding_prof_Ka = binding_prof_Ka[order(binding_prof_Ka)]
	Ka_total = rowSums(binding_prof_Ka)
	binding_prof_Ka = binding_prof_Ka/Ka_total
	binding_prof_Ka_cum = binding_prof_Ka
	AUC_total = 0
	for (i in 1:num_prot_genes){binding_prof_Ka_cum[1,i] = sum(binding_prof_Ka[1,1:i])}
	for (j in 1:(num_prot_genes-1)){
		#binding_prof_Ka[1,j] = sum(binding_prof_Ka[1,1:j])
		#binding_prof_Ka[1,j+1] = sum(binding_prof_Ka[1,1:(j+1)])
		AUC = (binding_prof_Ka_cum[1,j])*(1/num_prot_genes) + (1/2)*(binding_prof_Ka_cum[1,(j+1)]-binding_prof_Ka_cum[1,j])*(1/num_prot_genes)
		AUC_total = AUC_total + AUC
	}
	#Account for area from 0 to 1 
	AUC_total = AUC_total + (1/2)*(binding_prof_Ka_cum[1,1])*(1/num_prot_genes)
	#Calculate Gini coefficient using area under curve
	gini_coeff = 1 - (2*AUC_total)
	score_mat[2,run_num] = round(gini_coeff,3)
	
	#Method 3: Selectivity Score as Rank Sum Test
	
	#Calculate Wilcox test p-value for high-affinity target
	binding_prof_Kd_asc = binding_prof_Kd[order(binding_prof_Kd)]
	output_list = wilcox.test(binding_prof_Kd_asc[1,1], unlist(binding_prof_Kd_asc[1,2:num_prot_genes]), alternative = "less")
	p_value = output_list[[3]]
	if (p_value < 10^-15){p_value = 10^-15}
	adj_p_value = -log10(p_value)/15
	on_target_quart = quantile(binding_prof_Kd_asc[1,1], 0.25)
	off_target_quart = quantile(binding_prof_Kd_asc[1,2:num_prot_genes], 0.25)
	diff_value = (log10(off_target_quart) - log10(on_target_quart))/3
	research_bias = 1 - (1/num_prot_genes) #Only one measurement for on target, rest are for off target
	selectivity_score = (diff_value[[1]] + adj_p_value + research_bias)/3
	score_mat[3,run_num] = round(selectivity_score,2)
	
	#Calculate Wilcox test p-value for low-affinity target
	binding_prof_Kd_asc = binding_prof_Kd[order(binding_prof_Kd)]
	output_list <- wilcox.test(binding_prof_Kd_asc[1,(num_high_aff[[1]]+1)]%>%as.numeric(), unlist(binding_prof_Kd_asc[1,1:(num_prot_genes-1)]), alternative = "less")
	p_value = output_list[[3]]
	if (p_value < 10^-15){p_value = 10^-15}
	adj_p_value = -log10(p_value)/15
	on_target_quart = quantile(binding_prof_Kd_asc[1,(num_high_aff[[1]]+1)], 0.25)
	off_target_quart = quantile(binding_prof_Kd_asc[1,1:(num_prot_genes-1)], 0.25)
	diff_value = (log10(off_target_quart) - log10(on_target_quart))/3
	research_bias = 1 - (1/num_prot_genes) #Only one measurement for on target, rest are for off target
	selectivity_score = (diff_value[[1]] + adj_p_value + research_bias)/3
	score_mat[4,run_num] = round(selectivity_score,2)
	print(run_num)
}

#####################################################################################################################################################################T
# normalization  ------------
#####################################################################################################################################################################T

#for (i in 1:nrow(score_mat)){score_mat[i,] = score_mat[i,]/min(score_mat[i,])}

#####################################################################################################################################################################T
# plot table  ------------
#####################################################################################################################################################################Tß
score_mat_rev = melt(t(score_mat))
colnames(score_mat_rev) = list("Number_of_High_Affinity_Targets","Method","Values_raw")#"Values_Normalized_to_Min")
Values_raw_plot <- 
  ggplot(score_mat_rev,aes(x=Number_of_High_Affinity_Targets, y = Values_raw, color = Method)) + 
  geom_line() + 
  geom_point() + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5 ),
        axis.text.y = element_text(angle = 0, hjust = 0.5 ),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black") )+
  ggtitle("Comparisons Drug Selectivity Methods")
Values_raw_plot

#setwd(home_dir)
#ggsave("MethodComparisonValue.png",first_plot)


score_mat_dup = score_mat
for (k in 1:(length(diff_dists)-1))
{
	score_mat_dup[,k+1] = (((score_mat[,k+1] - score_mat[,k])/(score_mat[,k]))*100)
}
score_mat_dup[,1] = 0
score_mat_dup_rev = melt(t(score_mat_dup))
colnames(score_mat_dup_rev) = list("Number_of_High_Affinity_Targets","Method","Percent_Change")
second_plot <- 
  ggplot(score_mat_dup_rev,aes(x=Number_of_High_Affinity_Targets, y = Percent_Change, color = Method)) + 
  geom_line() + geom_point() + 
  ggtitle("Percent Change of Different Drug Selectivity Metrics Across A Varying Number of High Affinity Targets")

#setwd(home_dir)
#ggsave("MethodComparisonPercent.png", second_plot)