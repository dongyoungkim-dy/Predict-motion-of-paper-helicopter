

#actual test data
falltime = read.csv("final_test_data.csv", header = TRUE)
f1 <- as.numeric(unlist(falltime[1]))
f2 <- as.numeric(unlist(falltime[2]))
f3 <- as.numeric(unlist(falltime[3]))
f4 <- as.numeric(unlist(falltime[4]))

#predictied fall time
falltime_quad_model1 = read.csv("06. quad_model1 fall time 10.67m_1clips.csv", header = TRUE)
falltime_quad_model2 = read.csv("03. quad_model2 fall time 10.67m_2clips.csv", header = TRUE)
falltime_quad_model3 = read.csv("04. quad_model3 fall time 6.82m_1clips.csv", header = TRUE)
falltime_quad_model4 = read.csv("05. quad_model4 fall time 6.82m_2clips.csv", header = TRUE)
p1 <- as.numeric(unlist(falltime_quad_model1))
p2 <- as.numeric(unlist(falltime_quad_model2))
p3 <- as.numeric(unlist(falltime_quad_model3))
p4 <- as.numeric(unlist(falltime_quad_model4))

title1 = "model1 cdf comparison of actual falltime and predicted fall time"
df <- data.frame(x = c(f1, p1), ggg=factor(rep(1:2, c(21, 3999))))
ggplot(df, aes(x, colour = ggg)) + 
  stat_ecdf()+
  scale_colour_hue(name="my legend", labels=c('fall time','predicted'))+
  labs(title=title1)

title1 = "model2 cdf comparison of actual falltime and predicted fall time"
df <- data.frame(x = c(f1, p1), ggg=factor(rep(1:2, c(21, 3999))))
ggplot(df, aes(x, colour = ggg)) + 
  stat_ecdf()+
  scale_colour_hue(name="my legend", labels=c('fall time','predicted'))+
  labs(title=title1)

title1 = "model3 cdf comparison of actual falltime and predicted fall time"
df <- data.frame(x = c(f1, p1), ggg=factor(rep(1:2, c(21, 3999))))
ggplot(df, aes(x, colour = ggg)) + 
  stat_ecdf()+
  scale_colour_hue(name="my legend", labels=c('fall time','predicted'))+
  labs(title=title1)

title1 = "model4 cdf comparison of actual falltime and predicted fall time"
df <- data.frame(x = c(f1, p1), ggg=factor(rep(1:2, c(21, 3999))))
ggplot(df, aes(x, colour = ggg)) + 
  stat_ecdf()+
  scale_colour_hue(name="my legend", labels=c('fall time','predicted'))+
  labs(title=title1)
################################


title1 = "model1 cdf comparison of actual falltime and predicted fall time"
a1 <- data.frame(x = c(f1, p1), ggg=factor(rep(1:2, c(21, 3999)))) %>%
  ggplot(aes(x, colour = ggg)) + 
  stat_ecdf()+
  scale_colour_hue(name="my legend", labels=c('fall time','predicted'))+
  labs(title=title1)

title1 = "model2 cdf comparison of actual falltime and predicted fall time"
a2 <- data.frame(x = c(f2, p2), ggg=factor(rep(1:2, c(21, 3999)))) %>%
  ggplot(aes(x, colour = ggg)) + 
  stat_ecdf()+
  scale_colour_hue(name="my legend", labels=c('fall time','predicted'))+
  labs(title=title1)

title1 = "model3 cdf comparison of actual falltime and predicted fall time"
a3 <- data.frame(x = c(f3, p3), ggg=factor(rep(1:2, c(21, 3999)))) %>%
  ggplot(aes(x, colour = ggg)) + 
  stat_ecdf()+
  scale_colour_hue(name="my legend", labels=c('fall time','predicted'))+
  labs(title=title1)

title1 = "model4 cdf comparison of actual falltime and predicted fall time"
a4 <- data.frame(x = c(f4, p4), ggg=factor(rep(1:2, c(21, 3999)))) %>%
  ggplot(aes(x, colour = ggg)) + 
  stat_ecdf()+
  scale_colour_hue(name="my legend", labels=c('fall time','predicted'))+
  labs(title=title1)

ggarrange(a1, a2, a3, a4)
