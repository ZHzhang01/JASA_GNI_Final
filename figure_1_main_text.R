library(latex2exp)

# install.packages("gridExtra")

# Load the gridExtra package
library(gridExtra)

methods <- c("HT", "Haj", "F", ("L"), "ND-F", "ND-G_1", "ND-L",  "ND-G_2")
tauhat <- c(0.0082, 0.0146, 0.0170, 0.0168, 0.0167 ,0.0161, 0.0164 ,0.0175)
se <- c(0.0272, 0.0225, 0.0218, 0.0197, 0.0211, 0.0205, 0.0195, 0.0193)

lower_bounds <- tauhat - 1.96 * se
upper_bounds <- tauhat + 1.96 * se
midpoints <- (lower_bounds + upper_bounds) / 2


library(ggplot2)

data <- data.frame(direct = c("HT", "Haj", "F", ("L"), "ND-F", "ND-G_1", "ND-L",  "ND-G_2"),
                   Interval = midpoints,
                   lower_ci = lower_bounds,
                   upper_ci = upper_bounds)

data$direct <- factor(data$direct, levels = c("HT", "Haj", "F", ("L"), "ND-F", "ND-G_1", "ND-L",  "ND-G_2"), ordered = TRUE)

# Plot point estimates with 95 percent confidence intervals for the direct-effect comparison.
p1 <- ggplot(data, aes(x = direct, y = Interval, color = direct)) +
  geom_point(size = 2) +  # Solid points
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.5, size = 1) +  # Confidence intervals
  theme(axis.text.x = element_text(angle = 90, hjust = 1), aspect.ratio = 2/2) +  # Display x-axis labels vertically
  scale_color_manual(values = c("#A4C8FB", "#A3B4E8", "#A2A0D5", "#A19BC2", "#A18DAF", "#A1809C", "#A07489", "#9F6589"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  theme(legend.position = 'none')

tauhat <- c(0.0381, 0.0611, 0.0604, 0.0581, 0.0641, 0.0653, 0.0590, 0.0610)
se <- c(0.0447, 0.0292, 0.0270, 0.0241, 0.0258, 0.0250, 0.0236, 0.0233)






tauhat <- c(0.0381, 0.0611, 0.0604, 0.0581, 0.0660, 0.0686, 0.0592, 0.0603)
se <- c(0.0447, 0.0292, 0.0270, 0.0241, 0.0258, 0.0250, 0.0236, 0.0232)
lower_bounds <- tauhat - 1.96 * se
upper_bounds <- tauhat + 1.96 * se
midpoints <- (lower_bounds + upper_bounds) / 2

library(ggplot2)
data <- data.frame(spillover = c("HT", "Haj", "F", ("L"), "ND-F", "ND-G_1", "ND-L",  "ND-G_2"),
                   Interval = midpoints,       
                   lower_ci = lower_bounds,
                   upper_ci = upper_bounds)

data$spillover <- factor(data$spillover, levels = c("HT", "Haj", "F", ("L"), "ND-F", "ND-G_1", "ND-L",  "ND-G_2"), ordered = TRUE)

# Plot point estimates with 95 percent confidence intervals for the spillover-effect comparison.
p2 <- ggplot(data, aes(x = spillover, y = Interval, color = spillover)) +
  geom_point(size = 2) +  # Solid points
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.5, size = 1) +  # Confidence intervals
  theme(axis.text.x = element_text(angle = 90, hjust = 1), aspect.ratio = 2/2) +  # Display x-axis labels vertically
  scale_color_manual(values = c("#A4C8FB", "#A3B4E8", "#A2A0D5", "#A19BC2", "#A18DAF", "#A1809C", "#A07489", "#9F6589"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  theme(legend.position = 'none')



#grid.arrange(p1, p2, nrow = 1, widths = widths)
grid.arrange(p1, p2, nrow = 1)









































