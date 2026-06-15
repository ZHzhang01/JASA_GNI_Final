import matplotlib.pyplot as plt
import numpy as np
import os

methods = ["HT", "Haj", "F", "L", "ND-F", r"ND-$\phi_0(G_1)$", "ND-L", r"ND-$\phi_0(G_2)$"]
colors = ["#A4C8FB", "#A3B4E8", "#A2A0D5", "#A19BC2", "#A18DAF", "#A1809C", "#A07489", "#9F6589"]

tauhat1 = np.array([0.0082, 0.0146, 0.0170, 0.0168, 0.0167, 0.0161, 0.0164, 0.0175])
se1 = np.array([0.0272, 0.0225, 0.0218, 0.0197, 0.0211, 0.0205, 0.0195, 0.0193])
lower1 = tauhat1 - 1.96 * se1
upper1 = tauhat1 + 1.96 * se1

tauhat2 = np.array([0.0381, 0.0611, 0.0604, 0.0581, 0.0660, 0.0686, 0.0592, 0.0603])
se2 = np.array([0.0447, 0.0292, 0.0270, 0.0241, 0.0258, 0.0250, 0.0236, 0.0232])
lower2 = tauhat2 - 1.96 * se2
upper2 = tauhat2 + 1.96 * se2

x = np.arange(len(methods))

fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.8))

for ax, tauhat, lower, upper, xlabel in [
    (axes[0], tauhat1, lower1, upper1, "direct"),
    (axes[1], tauhat2, lower2, upper2, "spillover")
]:
    ax.set_facecolor("#EBEBEB")
    ax.set_axisbelow(True)
    ax.grid(True, which="major", color="white", linewidth=1.1)
    ax.grid(True, which="minor", color="white", linewidth=0.6)
    ax.minorticks_on()
    
    for i in range(len(methods)):
        ax.errorbar(
            x[i], tauhat[i],
            yerr=[[tauhat[i] - lower[i]], [upper[i] - tauhat[i]]],
            fmt="o",
            markersize=4,
            elinewidth=1.6,
            capsize=7,
            color=colors[i],
        )
    ax.axhline(0, linestyle="--", linewidth=1, color="black")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=90, ha="center", va="top")
    ax.set_xlabel(xlabel)
    ax.set_xlim(-0.5, len(methods) - 0.5)
    for spine in ax.spines.values():
        spine.set_visible(False)

fig.patch.set_facecolor("white")
fig.tight_layout(w_pad=2.5)

pdf_path = "/mnt/data/ci_intervals_direct_spillover_phi0_labels_ggplot_gray.pdf"
png_path = "/mnt/data/ci_intervals_direct_spillover_phi0_labels_ggplot_gray.png"

fig.savefig(pdf_path, bbox_inches="tight")
fig.savefig(png_path, dpi=300, bbox_inches="tight")
plt.close(fig)

print(pdf_path)
print(png_path)






# library(latex2exp)

# # install.packages("gridExtra")

# # Load the gridExtra package
# library(gridExtra)

# methods <- c("HT", "Haj", "F", ("L"), "ND-F", "ND-G_1", "ND-L",  "ND-G_2")
# tauhat <- c(0.0082, 0.0146, 0.0170, 0.0168, 0.0167 ,0.0161, 0.0164 ,0.0175)
# se <- c(0.0272, 0.0225, 0.0218, 0.0197, 0.0211, 0.0205, 0.0195, 0.0193)

# lower_bounds <- tauhat - 1.96 * se
# upper_bounds <- tauhat + 1.96 * se
# midpoints <- (lower_bounds + upper_bounds) / 2


# library(ggplot2)

# data <- data.frame(direct = c("HT", "Haj", "F", ("L"), "ND-F", "ND-G_1", "ND-L",  "ND-G_2"),
#                    Interval = midpoints,
#                    lower_ci = lower_bounds,
#                    upper_ci = upper_bounds)

# data$direct <- factor(data$direct, levels = c("HT", "Haj", "F", ("L"), "ND-F", "ND-G_1", "ND-L",  "ND-G_2"), ordered = TRUE)

# # Plot point estimates with 95 percent confidence intervals for the direct-effect comparison.
# p1 <- ggplot(data, aes(x = direct, y = Interval, color = direct)) +
#   geom_point(size = 2) +  # Solid points
#   geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.5, size = 1) +  # Confidence intervals
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), aspect.ratio = 2/2) +  # Display x-axis labels vertically
#   scale_color_manual(values = c("#A4C8FB", "#A3B4E8", "#A2A0D5", "#A19BC2", "#A18DAF", "#A1809C", "#A07489", "#9F6589"))+
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
#   theme(legend.position = 'none')

# tauhat <- c(0.0381, 0.0611, 0.0604, 0.0581, 0.0641, 0.0653, 0.0590, 0.0610)
# se <- c(0.0447, 0.0292, 0.0270, 0.0241, 0.0258, 0.0250, 0.0236, 0.0233)






# tauhat <- c(0.0381, 0.0611, 0.0604, 0.0581, 0.0660, 0.0686, 0.0592, 0.0603)
# se <- c(0.0447, 0.0292, 0.0270, 0.0241, 0.0258, 0.0250, 0.0236, 0.0232)
# lower_bounds <- tauhat - 1.96 * se
# upper_bounds <- tauhat + 1.96 * se
# midpoints <- (lower_bounds + upper_bounds) / 2

# library(ggplot2)
# data <- data.frame(spillover = c("HT", "Haj", "F", ("L"), "ND-F", "ND-G_1", "ND-L",  "ND-G_2"),
#                    Interval = midpoints,       
#                    lower_ci = lower_bounds,
#                    upper_ci = upper_bounds)

# data$spillover <- factor(data$spillover, levels = c("HT", "Haj", "F", ("L"), "ND-F", "ND-G_1", "ND-L",  "ND-G_2"), ordered = TRUE)

# # Plot point estimates with 95 percent confidence intervals for the spillover-effect comparison.
# p2 <- ggplot(data, aes(x = spillover, y = Interval, color = spillover)) +
#   geom_point(size = 2) +  # Solid points
#   geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.5, size = 1) +  # Confidence intervals
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), aspect.ratio = 2/2) +  # Display x-axis labels vertically
#   scale_color_manual(values = c("#A4C8FB", "#A3B4E8", "#A2A0D5", "#A19BC2", "#A18DAF", "#A1809C", "#A07489", "#9F6589"))+
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
#   theme(legend.position = 'none')



# #grid.arrange(p1, p2, nrow = 1, widths = widths)
# grid.arrange(p1, p2, nrow = 1)









































