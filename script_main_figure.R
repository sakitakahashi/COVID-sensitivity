## Read in source code
source("functions_sensitivity.R")

### Set the sensitivity in each category
Se.mild <- 0.6
Se.severe <- 0.95

## Set the ratio between Se.asx and Se.mild
Se_asx_to_Se_mild_ratio <- 2/3
Se.asx <- Se_asx_to_Se_mild_ratio * Se.mild

# Set the specificity
Sp <- 1.0

## Set the proportions of each category in the population
prop.asx_true <- 0.43
prop.mild_true <- 0.52
prop.severe_true <- 0.05

## Set the color scheme
colors_for_scenarios <- c("#7C14CC","#CC8A00","#0ACCC5", "lightgrey")

################################################################################

## [1] Generate scenarios
data.frame(which_scenario="1", prop.asx=0.43, prop.mild=0.52, prop.severe=0.05) %>%
	add_row(which_scenario="2", prop.asx=0.05, prop.mild=0.45, prop.severe=0.50) %>%
	add_row(which_scenario="3", prop.asx=0, prop.mild=0, prop.severe=1) -> scenarios

## Plot them
scenarios %>%
	## Add data on true clinical symptoms
	add_row(which_scenario="True clinical\noutcomes", prop.asx=prop.asx_true, prop.mild=prop.mild_true, prop.severe=prop.severe_true, .before=1) %>%
	melt() %>%
	mutate(variable=factor(variable, levels=c("prop.asx", "prop.mild", "prop.severe"))) %>%
	mutate(which_scenario=factor(which_scenario, levels=c("True clinical\noutcomes", "1", "2", "3"))) %>%
	ggplot(aes(fill=variable, x=which_scenario, y=value)) +
	geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
	scale_fill_manual(name="", breaks=c("prop.asx", "prop.mild", "prop.severe"), labels=c("Asymptomatic", "Mild", "Severe"), values=c("#97BAFC", "#EDFCCA", "#FDBCB1")) +
	xlab("Assay validation controls") +
	ylab("Proportion") +
	theme_classic() +
	theme(
		legend.position="top",
		axis.text=element_text(size=12),
		axis.title=element_text(size=12),
		legend.text=element_text(size=12)) -> plot_scenarios

## [2] Simulations over a range of true prevalence
true.prev <- seq(0, 0.5, by=0.025)
output.mat1 <- matrix(ncol=nrow(scenarios)+2, nrow=length(true.prev))
output.mat1[,1] <- true.prev
colnames(output.mat1) <- c("true.prev", "Uncorrected", "1", "2", "3")

for(i in 1:nrow(scenarios)) {
	
	for(j in 1:nrow(output.mat1)) {
		
		## Simulate observed prevalence (uncorrected)
		prev.sim <- sim.prev_symptoms(
			Se.asx=Se.asx,
			Se.mild=Se.mild,
			Se.severe=Se.severe,
			Sp=1,
			prop.asx=prop.asx_true,
			prop.mild=prop.mild_true,
			prop.severe=prop.severe_true,
			true.prev=true.prev[j])
		output.mat1[j,2] <- prev.sim
		
		## Calculate the weighted Se based on the distribution of positive controls
		Se.wtd <- Se.weighted_symptoms(
			Se.asx=Se.asx,
			Se.mild=Se.mild,
			Se.severe=Se.severe,
			v.asx=scenarios$prop.asx[i],
			v.mild=scenarios$prop.mild[i],
			v.severe=scenarios$prop.severe[i])
		
		## Calculate corrected prevalence
		adjust.prev <- adj.func(prev.obs=prev.sim, Se=Se.wtd, Sp=1)
		output.mat1[j,i+2] <- adjust.prev
		
	}
	
}

output.mat1 %>%
	as.data.frame() %>%
	gather(which_scenario, measured.prev, Unadjusted:`3`) %>%
	mutate(which_scenario = factor(which_scenario, levels=c("1", "2", "3", "Uncorrected"))) -> data_plot_sims

data_plot_sims %>%
	filter(true.prev>0) %>%
	mutate(ratio = measured.prev/true.prev) %>%
	group_by(which_scenario) %>%
	slice(1) %>%
	select(which_scenario, ratio) %>%
	as.data.frame() -> data_plot_sims_points

data_plot_sims %>%
	ggplot() +
		geom_segment(aes(x=0, y=0, xend=max(true.prev), yend=max(true.prev)*data_plot_sims_points$ratio[1]), colour=colors_for_scenarios[1], size=0.7) +
		geom_segment(aes(x=0, y=0, xend=max(true.prev), yend=max(true.prev)*data_plot_sims_points$ratio[2]), colour=colors_for_scenarios[2], size=0.7) +
		geom_segment(aes(x=0, y=0, xend=max(true.prev), yend=max(true.prev)*data_plot_sims_points$ratio[3]), colour=colors_for_scenarios[3], size=0.7) +
		geom_segment(aes(x=0, y=0, xend=max(true.prev), yend=max(true.prev)*data_plot_sims_points$ratio[4]), colour=colors_for_scenarios[4], size=0.7) +
		geom_point(aes(x=true.prev, y=measured.prev, group=which_scenario, shape=which_scenario), size=1.5) +
		xlab("True prevalence") +
		ylab("Measured prevalence") +
		theme_bw() +
		guides(shape=guide_legend(override.aes=list(size=2))) +
		labs(shape="Validation\nset used for\ncorrection") +
		annotate(geom="text", x=0.0, y=0.5, label=paste0("Se in severe: ", sprintf(Se.severe, fmt='%#.2f')), size=4, hjust=0) +
		annotate(geom="text", x=0.0, y=0.475, label=paste0("Se in mild: ", sprintf(Se.mild, fmt='%#.2f')), size=4, hjust=0) +
		annotate(geom="text", x=0.0, y=0.45, label=paste0("Se in asymptomatic: ", sprintf(Se.asx, fmt='%#.2f')), size=4, hjust=0) +
		scale_shape_manual(values=c(16,15,17,3)) +
		theme(
			aspect.ratio=1,
			legend.position="top",
			legend.title.align=0.5,
			panel.grid.minor=element_blank(),
			axis.text=element_text(size=12),
			axis.title=element_text(size=12),
			legend.text=element_text(size=12),
			legend.title=element_text(size=12)) -> plot_sims

## [3] Sweep over sensitivity in asx/mild for each scenario (fixing Se in severe to Se.severe)
true.prev.fixed <- 0.5 ## Doesn't matter what this is for calculating ratios
Se.mild_to_sweep <- seq(0.25, 0.95, by=0.05)
scenarios_to_sweep <- 1:nrow(scenarios)

output.mat2 <- matrix(ncol=4, nrow=length(Se.mild_to_sweep)*length(scenarios_to_sweep))
output.mat2[,1] <- rep(Se.mild_to_sweep, times=length(scenarios_to_sweep))
output.mat2[,2] <- rep(scenarios_to_sweep, each=length(Se.mild_to_sweep))
output.mat2[,3] <- true.prev.fixed
colnames(output.mat2) <- c("Se.mild_to_sweep", "scenarios_to_sweep", "true.prev.fixed", "adjust.prev")

for(j in 1:nrow(output.mat2)) {
	
	## Simulate observed prevalence (uncorrected)
	prev.sim <- sim.prev_symptoms(
		Se.asx=output.mat2[j,"Se.mild_to_sweep"]*Se_asx_to_Se_mild_ratio,
		Se.mild=output.mat2[j,"Se.mild_to_sweep"],
		Se.severe=Se.severe,
		Sp=1,
		prop.asx=prop.asx_true,
		prop.mild=prop.mild_true,
		prop.severe=prop.severe_true,
		true.prev=true.prev.fixed)
	
	## Calculate the weighted Se based on the distribution of positive controls
	Se.wtd <- Se.weighted_symptoms(
		Se.asx=as.numeric(output.mat2[j,"Se.mild_to_sweep"])*Se_asx_to_Se_mild_ratio,
		Se.mild=as.numeric(output.mat2[j,"Se.mild_to_sweep"]),
		Se.severe=Se.severe,
		v.asx=scenarios$prop.asx[output.mat2[j,"scenarios_to_sweep"]],
		v.mild=scenarios$prop.mild[output.mat2[j,"scenarios_to_sweep"]],
		v.severe=scenarios$prop.severe[output.mat2[j,"scenarios_to_sweep"]])
	
	## Calculated corrected prevalence
	adjust.prev <- adj.func(prev.obs=prev.sim, Se=Se.wtd, Sp=1)
	output.mat2[j,"adjust.prev"] <- adjust.prev
	
}

output.mat2 %>%
	as.data.frame() %>%
	mutate(ratio = adjust.prev/true.prev.fixed) %>%
	mutate(which_scenario = factor(scenarios_to_sweep)) %>%
	select(-scenarios_to_sweep) -> data_plot_sweep

## Single points to highlight on the figure (from panel B)
data_plot_sweep %>%
	mutate(Se.mild_to_sweep=factor(Se.mild_to_sweep)) %>%
	dplyr::filter(Se.mild_to_sweep==Se.mild) %>%
	mutate(color=colors_for_scenarios[1:3]) %>%
	mutate(Se.mild_to_sweep=as.numeric(as.character(Se.mild_to_sweep))) %>%
	as.data.frame() -> data_plot_sweep_single

data_plot_sweep %>%
	ggplot() +
		geom_point(aes(x=as.numeric(Se.mild_to_sweep), y=as.numeric(ratio), group=which_scenario, shape=which_scenario), size=1.5) +
		geom_point(data=data_plot_sweep_single, aes(x=Se.mild_to_sweep, y=ratio, colour=color), shape=1, size=4, stroke=1.2) +
		xlab("Sensitivity in mild") +
		ylab("Ratio, measured to true prevalence") +
		theme_bw() +
		scale_x_continuous(breaks=seq(0.2,1,by=0.1), limits=c(0.2,1.001)) +
		scale_y_continuous(breaks=seq(0.2,1,by=0.1), limits=c(0.2,1.001)) +
		scale_shape_manual(values=c(16,15,17)) +
		scale_colour_identity(breaks=data_plot_sweep_single$color) +
		annotate(geom="text", x=0.7, y=0.35, label=paste0("Se in severe: ", Se.severe), size=4) +
		annotate(geom="text", x=0.7, y=0.30, label=paste0("Se in asymptomatic: Se in mild x ", MASS::fractions(Se_asx_to_Se_mild_ratio)), size=4) +
		guides(shape=guide_legend(ncol=2, override.aes=list(size=2))) +
		theme(
			aspect.ratio=1,
			legend.position="none",
			panel.grid.minor=element_blank(),
			axis.text=element_text(size=12),
			axis.title=element_text(size=12)) -> plot_sweep

## Plot them together
plot_all_main <- plot_scenarios + plot_sims + plot_sweep + plot_layout(ncol=3) + plot_annotation(tag_levels='A')

# ggsave("main_figure.pdf", plot_all_main, width=35, height=20, units="cm")
