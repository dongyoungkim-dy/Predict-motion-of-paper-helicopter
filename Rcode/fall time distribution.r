# library
library(ggplot2)
library(dplyr)
library(hrbrthemes)

# Build dataset with different distributions
data <- data.frame(
  type = c(rep("p1", length(p1))),
  value = c(p1 )
)

# Represent it
p <- data %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2")) +
  geom_density() +
  theme_ipsum() +
  labs(fill="")

ggarrange(p)

