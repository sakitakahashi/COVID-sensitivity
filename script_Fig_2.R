## Read in source code
source("functions_sensitivity.R")

### Set the sensitivity in each category (for last curve of panel C)
Se.asx <- 0.4
Se.mild <- 0.6
Se.severe <- 0.95

# Set the specificity
Sp <- 1.0

## Set the proportions of each category in the population
prop.asx_true <- 0.43
prop.mild_true <- 0.52
prop.severe_true <- 0.05

## Decide cutoffs for times since infection (in days)
## (1) Very recent: 0 to VERY_RECENT
## (2) Recent: VERY_RECENT to RECENT
## (3) Not recent: RECENT onwards
VERY_RECENT <- 60
RECENT <- 180

## Pick values of the reduction in Se at the different times since infection
mult.Se.2 <- 0.8
mult.Se.3 <- 0.6

################################################################################

## [S1] Draw an epi curve

## Parameter values
beta_fix <- 0.000003
gamma_fix <- 1/7
sigma_fix <- 1/3.5
N_fix <- 100000

R0 <- beta_fix*N_fix/(gamma_fix)
R0

## Run the ODE solver
lsoda(y=c(S=N_fix-1, E=0, I=1, R=0, C=0), times=seq(from=0, to=365, by=1), func=dx.dt.SEIR, parms=c(beta=beta_fix, gamma=gamma_fix, sigma=sigma_fix)) %>%
	as.data.frame() -> SEIR.output

## Pick times to highlight, systematically
SEIR.output %>%
	mutate(new_cases=C-lag(C, default=0)) %>%
	mutate(new_cases_round = ceiling(new_cases)) -> SEIR.output_2

TSI <- rep(SEIR.output_2$time, times=SEIR.output_2$new_cases_round)

## Get TSI for Scenario 1 (day 60)
N1 <- 60
N1-TSI[TSI<=N1] %>% quantile(probs=seq(0,1,by=0.01)) -> TSI_1
ecdf_TSI_1 <- ecdf(TSI_1)

## Get TSI for Scenario 2 (day 180)
N2 <- 180
N2-TSI[TSI<=N2] %>% quantile(probs=seq(0,1,by=0.01)) -> TSI_2
ecdf_TSI_2 <- ecdf(TSI_2)

## Get TSI for Scenario 3 (day 300)
N3 <- 300
N3-TSI[TSI<=N3] %>% quantile(probs=seq(0,1,by=0.01)) -> TSI_3
ecdf_TSI_3 <- ecdf(TSI_3)

SEIR.output %>%
	filter(time %in% c(N1, N2, N3)) %>%
	mutate(which_scenario=c(paste0("Day ", N1), paste0("Day ", N2), paste0("Day ", N3))) -> SEIR.output.highlight

SEIR.output %>%
	ggplot(aes(x=time, y=E/N_fix)) +
		geom_line(size=0.8) +
		geom_point(data=SEIR.output.highlight, aes(x=time, y=E/N_fix), size=3) +
		geom_text_repel(data=SEIR.output.highlight, aes(x=time, y=E/N_fix, label=which_scenario), size=4.5, segment.colour=NA, box.padding=1) +
		xlab("Time (days)") +
		ylab("Proportion infected") +
		theme_classic(base_size=16) +
		theme(
			legend.position="none",
			panel.grid.minor=element_blank()) -> plot_SEIR

## [S2] Generate scenarios for the epidemic timing (TSI = time since infection)
data.frame(which_scenario=paste0("Day ", N1), prop.1=1, prop.2=0, prop.3=0) %>%
	add_row(which_scenario=paste0("Day ", N2), prop.1=ecdf_TSI_2(VERY_RECENT), prop.2=ecdf_TSI_2(RECENT)-ecdf_TSI_2(VERY_RECENT), prop.3=(1-ecdf_TSI_2(RECENT))) %>%
	add_row(which_scenario=paste0("Day ", N3), prop.1=ecdf_TSI_3(VERY_RECENT), prop.2=ecdf_TSI_3(RECENT)-ecdf_TSI_3(VERY_RECENT), prop.3=(1-ecdf_TSI_3(RECENT))) -> scenarios_TSI

## Which assay validation controls?
scenarios_TSI %>%
	add_row(which_scenario="Assay\nvalidation\ncontrols", prop.1=1, prop.2=0, prop.3=0, .after=4) -> scenarios_TSI_controls

## Plot them
scenarios_TSI_controls %>%
	melt() %>%
	mutate(variable=factor(variable, levels=c("prop.1", "prop.2", "prop.3"))) %>%
	mutate(which_scenario=factor(which_scenario, levels=c(paste0("Day ", N1), paste0("Day ", N2), paste0("Day ", N3), "Assay\nvalidation\ncontrols"))) %>%
	ggplot(aes(fill=variable, x=which_scenario, y=value)) +
	geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
	scale_fill_manual(name="Time since\ninfection", breaks=c("prop.1", "prop.2", "prop.3"), values=c("grey10", "grey45", "grey80"), labels=c(paste0("0-", VERY_RECENT, " days"), paste0(VERY_RECENT, "-", RECENT, " days"), paste0(RECENT, "+ days"))) +
	xlab("Time point in the epidemic") +
	ylab("Proportion") +
	guides(fill=guide_legend(override.aes=list(size=2), nrow=3)) +
	theme_classic(base_size=16) +
	theme(legend.position="top") -> plot_scenarios_TSI

## [S3] Simulations over a range of true prevalence
true.prev <- seq(0, 0.5, by=0.025)
output.mat3 <- matrix(ncol=nrow(scenarios_TSI)+2, nrow=length(true.prev))
output.mat3[,1] <- true.prev
colnames(output.mat3) <- c("true.prev", paste0("Day ", N1), paste0("Day ", N2), paste0("Day ", N3), paste0("Day ", N2, " & clinical spectrum"))

for(i in 1:nrow(scenarios_TSI)) {
	
	for(j in 1:length(true.prev)) {
		
		## Simulate observed prevalence (unadjusted), assuming that earlier infection reduces Se by mult.Se.2 & mult.Se.3
		prev.sim <- sim.prev_TSI(
			Se.1=Se.severe,
			Se.2=Se.severe*mult.Se.2,
			Se.3=Se.severe*mult.Se.3,
			Sp=1,
			prop.1=scenarios_TSI[i,"prop.1"],
			prop.2=scenarios_TSI[i,"prop.2"],
			prop.3=scenarios_TSI[i,"prop.3"],
			true.prev=true.prev[j])
		
		## Calculate the weighted Se based on the distribution of positive controls
		Se.wtd <- Se.weighted_TSI(
			Se.1=Se.severe,
			Se.2=Se.severe*mult.Se.2,
			Se.3=Se.severe*mult.Se.3,
			v.1=scenarios_TSI_controls[4,"prop.1"],
			v.2=scenarios_TSI_controls[4,"prop.2"],
			v.3=scenarios_TSI_controls[4,"prop.3"])
		
		## Calculated adjusted prevalence
		adjust.prev <- adj.func(prev.obs=prev.sim, Se=Se.wtd, Sp=1)
		output.mat3[j,i+1] <- adjust.prev
		
	}
	
}

## Also accounting for clinical spectrum
data.frame(
	which_scenario=paste0("Day ", N2, " & clinical spectrum"),
	prop.asx.1=ecdf_TSI_2(VERY_RECENT)*prop.asx_true,
	prop.asx.2=(ecdf_TSI_2(RECENT)-ecdf_TSI_2(VERY_RECENT))*prop.asx_true,
	prop.asx.3=0,
	prop.mild.1=ecdf_TSI_2(VERY_RECENT)*prop.mild_true,
	prop.mild.2=(ecdf_TSI_2(RECENT)-ecdf_TSI_2(VERY_RECENT))*prop.mild_true,
	prop.mild.3=0,
	prop.severe.1=ecdf_TSI_2(VERY_RECENT)*prop.severe_true,
	prop.severe.2=(ecdf_TSI_2(RECENT)-ecdf_TSI_2(VERY_RECENT))*prop.severe_true,
	prop.severe.3=0) -> scenarios_TSI_symptoms

## Which assay validation controls?
scenarios_TSI_symptoms %>%
	add_row(
		which_scenario="Assay\nvalidation\ncontrols", 
		prop.asx.1=0,
		prop.asx.2=0,
		prop.asx.3=0,
		prop.mild.1=0,
		prop.mild.2=0,
		prop.mild.3=0,
		prop.severe.1=1,
		prop.severe.2=0,
		prop.severe.3=0,
		.after=1) -> scenarios_TSI_symptoms_controls

for(j in 1:length(true.prev)) {
	
	## Simulate observed prevalence (unadjusted), assuming that earlier infection reduces Se by mult.Se.2 & mult.Se.3
	prev.sim <- sim.prev_TSI_symptoms(
		Se.asx.1=Se.asx,
		Se.asx.2=Se.asx*mult.Se.2,
		Se.asx.3=Se.asx*mult.Se.3,
		Se.mild.1=Se.mild,
		Se.mild.2=Se.mild*mult.Se.2,
		Se.mild.3=Se.mild*mult.Se.3,
		Se.severe.1=Se.severe,
		Se.severe.2=Se.severe*mult.Se.2,
		Se.severe.3=Se.severe*mult.Se.3,
		Sp=1,
		prop.asx.1=scenarios_TSI_symptoms[1,"prop.asx.1"],
		prop.asx.2=scenarios_TSI_symptoms[1,"prop.asx.2"],
		prop.asx.3=scenarios_TSI_symptoms[1,"prop.asx.3"],
		prop.mild.1=scenarios_TSI_symptoms[1,"prop.mild.1"],
		prop.mild.2=scenarios_TSI_symptoms[1,"prop.mild.2"],
		prop.mild.3=scenarios_TSI_symptoms[1,"prop.mild.3"],
		prop.severe.1=scenarios_TSI_symptoms[1,"prop.severe.1"],
		prop.severe.2=scenarios_TSI_symptoms[1,"prop.severe.2"],
		prop.severe.3=scenarios_TSI_symptoms[1,"prop.severe.3"],
		true.prev=true.prev[j])
	
	## Calculate the weighted Se based on the distribution of positive controls
	Se.wtd <- Se.weighted_symptoms_TSI(
		Se.asx.1=Se.asx,
		Se.asx.2=Se.asx*mult.Se.2,
		Se.asx.3=Se.asx*mult.Se.3,
		Se.mild.1=Se.mild,
		Se.mild.2=Se.mild*mult.Se.2,
		Se.mild.3=Se.mild*mult.Se.3,
		Se.severe.1=Se.severe,
		Se.severe.2=Se.severe*mult.Se.2,
		Se.severe.3=Se.severe*mult.Se.3,
		v.asx.1=scenarios_TSI_symptoms_controls[2,"prop.asx.1"],
		v.asx.2=scenarios_TSI_symptoms_controls[2,"prop.asx.2"],
		v.asx.3=scenarios_TSI_symptoms_controls[2,"prop.asx.3"],
		v.mild.1=scenarios_TSI_symptoms_controls[2,"prop.mild.1"],
		v.mild.2=scenarios_TSI_symptoms_controls[2,"prop.mild.2"],
		v.mild.3=scenarios_TSI_symptoms_controls[2,"prop.mild.3"],
		v.severe.1=scenarios_TSI_symptoms_controls[2,"prop.severe.1"],
		v.severe.2=scenarios_TSI_symptoms_controls[2,"prop.severe.2"],
		v.severe.3=scenarios_TSI_symptoms_controls[2,"prop.severe.3"])
	
	## Calculate adjusted prevalence
	adjust.prev <- adj.func(prev.obs=prev.sim, Se=Se.wtd, Sp=1)
	output.mat3[j,5] <- adjust.prev
	
}

output.mat3 %>%
	as.data.frame() %>%
	gather(which_scenario, measured.prev, paste0("Day ", N1):paste0("Day ", N2, " & clinical spectrum")) %>%
	mutate(which_scenario = factor(which_scenario, levels=c(paste0("Day ", N1), paste0("Day ", N2), paste0("Day ", N3), paste0("Day ", N2, " & clinical spectrum")))) -> data_plot_sims_TSI

data_plot_sims_TSI %>%
	ggplot() +
		geom_point(aes(x=true.prev, y=measured.prev, group=which_scenario, shape=which_scenario), size=2) +
		xlab("True prevalence") +
		ylab("Estimated prevalence") +
		theme_bw(base_size=16) +
		scale_shape_manual(name="Time point in\nthe epidemic", breaks=c(paste0("Day ", N1), paste0("Day ", N2), paste0("Day ", N3), paste0("Day ", N2, " & clinical spectrum")), labels=c(paste0("Day ", N1), paste0("Day ", N2), paste0("Day ", N3), paste0("Day ", N2, " &\nclinical spectrum")), values=c(1,0,2,7)) +
		guides(shape=guide_legend(override.aes=list(size=2), nrow=2, byrow=TRUE)) +
		annotate(geom="text", x=0.5, y=0.062, label=paste0("Se, ", "0-", VERY_RECENT, " days: x ", sprintf(1, fmt='%#.2f')), size=4.5, hjust=1, vjust=0) +
		annotate(geom="text", x=0.5, y=0.031, label=paste0("Se, ", VERY_RECENT, "-", RECENT, " days: x ", sprintf(mult.Se.2, fmt='%#.2f')), size=4.5, hjust=1, vjust=0) +
		annotate(geom="text", x=0.5, y=0, label=paste0("Se, ", RECENT, "+ days: x ", sprintf(mult.Se.3, fmt='%#.2f')), size=4.5, hjust=1, vjust=0) +
		theme(
			aspect.ratio=1,
			legend.position="top",
			legend.title.align=0.5,
			panel.grid.minor=element_blank()) -> plot_sims_TSI

## Plot them together
plot_all_fig_2 <- plot_SEIR + plot_scenarios_TSI + plot_sims_TSI + plot_layout(ncol=3) + plot_annotation(tag_levels='A', theme=theme(plot.title=element_text(size=20)))

# ggsave(paste0("fig_2.pdf"), plot_all_fig_2, width=35, height=20, units="cm")
