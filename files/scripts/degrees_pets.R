library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)

dat <- tibble::tribble(
  ~Year, ~Degrees, ~Dollars,
   2000,                                                                                      5853,                                                                        39.7,
   2001,                                                                                      5694,                                                                        41.9,
   2002,                                                                                      5695,                                                                        44.6,
   2003,                                                                                      5696,                                                                        46.8,
   2004,                                                                                      5942,                                                                        49.8,
   2005,                                                                                      6366,                                                                        53.1,
   2006,                                                                                      6649,                                                                        56.9,
   2007,                                                                                      7187,                                                                        61.8,
   2008,                                                                                      7798,                                                                        65.7,
   2009,                                                                                      8026,                                                                        67.1
  ) %>%
  mutate(Dollars= Dollars*100+2000) 

dat_long <- dat %>%
  gather(Variable, Value, -Year)
  

cor(dat[,-1])

ggplot(dat_long, 
       aes(x = Year, y = Value, color = Variable, shape = Variable)) +
  geom_point(size = 3) +
  geom_line(size = 2) +
  scale_y_continuous(name = "Biology/biomedical doctorates awarded (US)\nDegrees awarded (National Science Foundation)",
                     sec.axis = sec_axis(~., name = "Money spent on pets (US)\nBillions of dollars (Bureau of Economic Analysis)",
                                         label = function(x) round((x-2000)/100))) +
  theme_bw(base_size = 17) +
  scale_color_manual(values = c("#8B3E2F", "green4")) +
  ggtitle("", subtitle = "Correlation = 0.95")+
  theme(legend.position = "top")
