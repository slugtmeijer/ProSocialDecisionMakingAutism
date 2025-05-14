library(ggplot2)
library(reshape2)
library(dplyr)

bic_data <- read.csv("bic_data.csv")

# Get the model names (excluding the 'Group' column)
model_names <- colnames(bic_data)[-1] 

# Find the index of the model with the lowest BIC for each subject
bic_data$winning_model <- apply(bic_data[, model_names], 1, which.min)

# Count the frequency of each winning model
model_counts <- table(bic_data$winning_model)

# Convert to percentage
model_percentages <- (model_counts / sum(model_counts)) * 100

# Convert to a data frame for plotting
model_percent_df <- data.frame(Model = names(model_percentages),
                               Percentage = as.numeric(model_percentages))

# Pie chart
ggplot(model_percent_df, aes(x = "", y = Percentage, fill = Model)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +  # Turns bar chart into a pie chart
  theme_minimal() +
  labs(title = "Winning Model Distribution", fill = "Model") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# Winning model based on model that wins most often
# overall_winning_model <- model_names[as.numeric(names(sort(model_counts, decreasing = TRUE))[1])]

# Winning model based on sum score
# Sum the BIC scores for each model across all subjects
total_bic_per_model <- colSums(bic_data[, model_names])

# Group by 'Group' and sum the BIC scores for each model
bic_sums_by_group <- bic_data %>%
  group_by(Group) %>%
  summarise(across(all_of(model_names), sum, na.rm = TRUE))

print(bic_sums_by_group)

# Find the model with the lowest total BIC score
overall_winning_model <- names(which.min(total_bic_per_model))

# Compute relative BIC scores per subject
for (model in model_names) {
  bic_data[[paste0(model, "_relative")]] <- bic_data[[model]] - bic_data[[overall_winning_model]]
}

# Melt the data for ggplot
bic_relative_melted <- melt(
  bic_data, 
  id.vars = "Group", 
  measure.vars = paste0(model_names, "_relative"),
  variable.name = "Model", 
  value.name = "Relative_BIC"
)

# Remove "_relative" from model names for cleaner plot labels
bic_relative_melted$Model <- gsub("_relative", "", bic_relative_melted$Model)

# Convert 'group' to a factor with correct labels
bic_relative_melted$group <- factor(bic_relative_melted$Group, 
                                    levels = c(1, 2), 
                                    labels = c("ASD", "Neurotypical"))

# Preserve the original model order
bic_relative_melted$Model <- factor(bic_relative_melted$Model, levels = model_names)

# **Summing the Relative BIC per Model and Group**
bic_summed <- bic_relative_melted %>%
  group_by(Model, group) %>%
  summarise(Sum_Relative_BIC = sum(Relative_BIC), .groups = "drop")

# **Plot the summed Relative BIC values**
ggplot(bic_summed, aes(x = Model, y = Sum_Relative_BIC, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Side-by-side bars
  theme_minimal() +
  labs(
    title = "Summed Relative BIC Scores per Model and Group",
    x = "Model",
    y = "Summed Relative BIC",
    fill = "Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability