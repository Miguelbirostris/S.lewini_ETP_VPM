
#These scripts are intended as supplementary code for the main file "VPA.Slewini.anotated.R". Contains code for several data visualization options and plots that were ultimately discarted for the publication of:

#Gomez-Garcia et al (in press) "Virtual Population Analysis of the critically endangered Scalloped Hammerhead (Sphyrna lewini) in the Eastern Tropical Pacific" 

#This code has not been optimized for speed

#Created by Miguel de Jesus Gomez Garcia
#Created: 05-April-2026
#Last edited: 07-Apr-2026
#First uploaded to GitHub on 07-Apr-2026





#Additional plots


# Calculate the mean mortality across stages for each draw
mean_Z <- rowMeans(Z_all)
mean_M <- rowMeans(M_all)
mean_Fs <- rowMeans(Fs_all)

Mortality_df<-data.frame(Mortality =c(
  rep("Total", length(mean_Z)),
  rep("Natural", length(mean_M)),
  rep("Fishing", length(mean_Fs))),
  Estimate=c(mean_Z,mean_M,mean_Fs)
)

# Mean mortalities for ploting
mortality_means <- Mortality_df %>%
  group_by(Mortality) %>%
  summarise(mean_est = mean(Estimate))


#Plot mortality densities

ggplot(data=Mortality_df, aes(x = Estimate)) +
  geom_density(fill = "gold", alpha = 0.6) +
  geom_vline(data=mortality_means, linetype = "dashed",aes(xintercept = mean_est)) +
  labs(title = "Sample Distribution of Mortalities",
       x = "Mortality estimates", y = "Density") +
  theme_minimal()+
  facet_wrap(~Mortality,scales="free_y",ncol=1)


#Densities per age class

Z_long <- as.data.frame(Z_all) %>%
  setNames(paste0("Stage_", 1:8)) %>%  # Rename columns for clarity
  pivot_longer(cols = everything(), names_to = "Stage", values_to = "Estimate") %>%
  mutate(Mortality = "Total")

M_long <- as.data.frame(M_all) %>%
  setNames(paste0("Stage_", 1:8)) %>%
  pivot_longer(cols = everything(), names_to = "Stage", values_to = "Estimate") %>%
  mutate(Mortality = "Natural")

Fs_long <- as.data.frame(Fs_all) %>%
  setNames(paste0("Stage_", 1:8)) %>%
  pivot_longer(cols = everything(), names_to = "Stage", values_to = "Estimate") %>%
  mutate(Mortality = "Fishing")

Mortality_long <- bind_rows(Z_long, M_long, Fs_long)

Mortality_short <-dplyr::filter(Mortality_long, Stage== c("Stage_1","Stage_2","Stage_3"))%>%
  mutate(Stage = recode(Stage,
                        "Stage_1" = "Neonates",
                        "Stage_2" = "Juveniles",
                        "Stage_3" = "Adults"))


ggplot(Z_long, aes(x = Estimate, fill = Stage)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  facet_wrap(~Stage, scales = "free_y", ncol = 1) +
  labs(title = "Stage-specific Mortality Density by Type",
       x = "Mortality Estimate", y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(Mortality_long, aes(x = Estimate, fill = Stage)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  facet_wrap(~Mortality, scales = "free_y", ncol = 1) +
  labs(title = "Stage-specific Mortality Density by Type",
       x = "Mortality Estimate", y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(Mortality_short, aes(x = Estimate, fill = Stage)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  facet_wrap(~Mortality, scales = "free_y", ncol = 1) +
  labs(title = "Stage-specific Mortality Density by Type",
       x = "Mortality Estimate", y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")


# Overall fecundity
f_all <- out_matrix$BUGSoutput$sims.list$f  

# Vectorize fecundity. single value for each draw
mean_f <- rowMeans(f_all)



#Itterative population

# Plot with ribbon for 95% credible interval

ggplot(popb_summary, aes(x = Year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = mean), color = "darkblue", size = 1) +
  labs(title = "Population Projection with Uncertainty",
       x = "Year", y = "Population") +
  theme_minimal()



### Individual Violin Plot ------------------------------------------------------------

base_violin_df<-popb_proj_df|>
  group_by(Iteration)|>
  summarise(Growth_Rate=first(lambda))

#shades or ablines

v_plot_full <- ggplot(base_violin_df, aes(x = "", y = Growth_Rate)) +
  geom_violin(fill = viridis(12, option = "mako")[10]) +
  geom_boxplot(width = 0.15,
               color = viridis(12, option = "mako")[1],  
               fill  = viridis(12, option = "mako")[8],
               linewidth = 1.1) +
  geom_point(
    stat = "summary",
    fun = "mean",
    color = viridis(12, option = "mako")[1],  # dark cyan
    fill  = viridis(12, option = "mako")[6],
    shape = 21,
    size = 10
  ) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "red",
             linewidth = 1)+
  theme_classic(base_size = 18) +
  theme() +
  labs(y = "Lambda", x="Hammerheads")


print (v_plot_full)

ggsave("ViolinPlot_Full.jpg", plot = v_plot_full,
       width = 12, height = 10, dpi = 300)



# Additional critical value plots -----------------------------------------
#Extract critical densities

pop_proj_density_critical<-popb_proj_dftrim2%>%group_by(Iteration)%>%
  summarise(Natural = mean(Natural),
            Fishing= mean(Fishing),
            Total = mean(Total),
            Fecundity = mean(Fecundity),
            lambda = mean(lambda))


ggplot(data=pop_proj_density_critical, aes(x = lambda)) +
  geom_density(fill = "orange", alpha = 0.6) +
  geom_vline(xintercept = mean(pop_proj_density_critical$lambda), linetype = "dashed") +
  labs(title = "Sample Distribution of lambda_critical",
       x = "Lambda", y = "Density") +
  theme_minimal()

ggplot(data=pop_proj_density_critical, aes(x = Fecundity)) +
  geom_density(fill = "orange", alpha = 0.6) +
  geom_vline(xintercept = mean(pop_proj_density_critical$Fecundity), linetype = "dashed") +
  labs(title = "Sample Distribution of Fencundity_critical",
       x = "Critical fecundity", y = "Density") +
  theme_minimal()



Mortality_df_critical <- pop_proj_density_critical %>%
  pivot_longer(cols = c(Total, Natural, Fishing),
               names_to = "Mortality",
               values_to = "Estimate")

# Mean mortalities for ploting
mortality_means_critical <- Mortality_df_critical %>%
  group_by(Mortality) %>%
  summarise(mean_est = mean(Estimate))


#Plot mortality densities

mortality_crit<-ggplot(data=Mortality_df_critical, aes(x = Estimate)) +
  geom_density(fill = "gold", alpha = 0.6) +
  geom_vline(data=mortality_means_critical, linetype = "dashed",aes(xintercept = mean_est)) +
  labs(title = "",
       x = "Mortality estimates", y = "Density") +
  theme_minimal(base_size = 30)+
  facet_wrap(~Mortality,scales="free_y",ncol=1)+
  theme(strip.text = element_text(size = 32, face = "bold"))

print(mortality_crit)

ggsave("Total_critical_mortality_density.jpg", plot = mortality_crit,
       width = 12, height = 10, dpi = 300)


#Plot mortality density by stage


stage_pattern <- "(Neonates|Juvenile_Males|Subadult_Males|Adult_Males|Juvenile_Females|Subadult_Females|Adult_Females|Resting_Females)"


pop_proj_density_critical_stage <- popb_proj_df %>%
  group_by(Iteration) %>%
  summarise(across(c(starts_with("Natural"),
                     starts_with("Fishing"),
                     starts_with("Total"),
                     Fecundity, lambda),
                   mean, na.rm = TRUE))


Mortality_stage_df_critical <- pop_proj_density_critical_stage %>%
  pivot_longer(
    cols = matches(stage_pattern),
    names_to = "StageFull",
    values_to = "Estimate"
  ) %>%
  mutate(
    Mortality = str_extract(StageFull, "^(Natural|Fishing|Total)"),
    Stage = str_extract(StageFull, stage_pattern)
  ) %>%
  select(Iteration, Mortality, Stage, Estimate)


# Plot
mortality_plot <- ggplot(Mortality_stage_df_critical, aes(x = Estimate, fill = Stage)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  facet_wrap(~Mortality, scales = "free_y", ncol = 1) +
  labs(
    title = "",
    x = "Mortality Estimate",
    y = "Density"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 32, face = "bold")
  )

print(mortality_plot)



ggsave("stage_specific_mortality_density.jpg", plot = mortality_plot,
       width = 12, height = 10, dpi = 300)





Mortality_split <- Mortality_stage_df_critical %>%
  mutate(
    Sex = case_when(
      str_detect(Stage, "Males")   ~ "Male",
      str_detect(Stage, "Females") ~ "Female",
      TRUE                         ~ "Neonates"
    ),
    Stage_clean = str_remove(Stage, "_Males|_Females")
  ) %>%
  rename(Stage_original = Stage,
         Stage = Stage_clean)


mortality_plot_merged <- ggplot(Mortality_split, aes(x = Estimate, fill = Stage)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  facet_wrap(~Mortality, scales = "free_y", ncol = 1) +
  labs(
    title = "",
    x = "Mortality Estimate",
    y = "Density"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 32, face = "bold")
  )

print(mortality_plot_merged)


ggsave("stage_specific_mortality_density_merged.jpg", plot = mortality_plot_merged,
       width = 12, height = 10, dpi = 300)

##Simmulated vs Critical VPs mortality ----------------------------------------------

stage_pattern <- "(Neonates|Juvenile_Males|Subadult_Males|Adult_Males|Juvenile_Females|Subadult_Females|Adult_Females|Resting_Females)"


sim_df <- sim_pop %>%
  dplyr::select(-Population, -SubPop, -Year) %>%
  pivot_longer(
    cols = matches(stage_pattern),
    names_to = "StageFull",
    values_to = "Estimate"
  ) %>%
  mutate(
    Mortality = str_extract(StageFull, "^(Natural|Fishing|Total)"),
    Stage     = str_extract(StageFull, stage_pattern)
  ) %>%
  select(Iteration, Mortality, Stage, Estimate, lambda, Fecundity)


PrevDF <- popb_proj_df %>%
  group_by(Iteration) %>%
  summarise(across(c(starts_with("Natural"),
                     starts_with("Fishing"),
                     starts_with("Total"),
                     Fecundity, lambda), mean, na.rm = TRUE))



Full_df <- PrevDF %>%
  pivot_longer(
    cols = matches(stage_pattern),
    names_to = "StageFull",
    values_to = "Estimate"
  ) %>%
  mutate(
    Mortality = str_extract(StageFull, "^(Natural|Fishing|Total)"),
    Stage     = str_extract(StageFull, stage_pattern)
  ) %>%
  select(Iteration, Mortality, Stage, Estimate, lambda, Fecundity)


sim_Mortality_plot <- ggplot(sim_df, aes(x = Estimate, fill = Stage)) +
  geom_density(alpha = 0.8, linewidth = 0.8) +
  facet_wrap(~Mortality, scales = "free_y", ncol = 1) +
  labs(title = "",
       x = "Mortality Estimate", y = "Density") +
  theme_minimal(base_size = 20) +
  scale_fill_viridis_d(option = "mako",end = 0.3) +   # mako palette for stages
  theme(legend.position = "bottom",
        strip.text = element_text(size = 32, face = "bold") )

print(sim_Mortality_plot)


ggsave("stage_specific_mortality_density_Current.jpg", plot = sim_Mortality_plot,
       width = 12, height = 10, dpi = 300)

#Critical mortality
print(mortality_plot)

#mixed

sim_df$Model<-"Current"
Full_df$Model <-"Full"
Critical_df<-Full_df[Full_df$lambda>0.975,]
Critical_df$Model <-"Critical"


mixed_df<-rbind(sim_df,Full_df,Critical_df)


Mixed_Mortality_plot <- ggplot(mixed_df, aes(x = Estimate, fill = Model)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  facet_wrap(~Mortality, scales = "free_y", ncol = 1) +
  labs(title = "",
       x = "Mortality Estimate", y = "Density") +
  theme_minimal(base_size = 20) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 32, face = "bold") )

print(Mixed_Mortality_plot)

#ggsave("stage_specific_mortality_density.jpg", plot = mortality_plot,
#       width = 12, height = 10, dpi = 300)


#Pie chart by Stage

Pie_stage_elas<-ie_stage_elas<-Av_stage_elas |>
  filter(sex != "Neonates") |> 
  ggplot( aes(x = "", y = elasticity, fill = Stage)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # Add white borders for separation
  coord_polar("y", start = 0) +  # Convert to pie chart
  scale_fill_viridis_d(option = "mako") +   # Custom colors
  labs(
    title = "Stage Elasticity",
    fill = "Stages"
  ) +
  theme_void() +  # Clean layout by removing background elements
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and style the title
    legend.position = "bottom",  # Position the legend
    legend.title = element_text(face = "bold")  # Bold legend title
  ) +
  facet_wrap(~ sex, scales = "free_y", ncol = 2)  # Facet by Divers

print(Pie_stage_elas)

#Pie chart by transition

Pie_Trans_elas<-Av_Trans_elas |>
  filter(sex != "Neonates") |> 
  ggplot(aes(x = "", y = elasticity, fill = Transition_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # Add white borders for separation
  coord_polar("y", start = 0) +  # Convert to pie chart
  scale_fill_viridis_d(option = "mako") +   # Custom colors
  labs(
    title = "Ttransition Elasticity",
    fill = "Stages"
  ) +
  theme_void(base_size = 30) +  # Set base text size
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and style the title
    legend.position = "bottom",  # Position the legend
    legend.title = element_text(face = "bold")  # Bold legend title
  ) +
  facet_wrap(~ sex, scales = "free_y", ncol = 2)  # Facet by Divers




print(Pie_Trans_elas)

