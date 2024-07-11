install.packages('devtools')
##devtools::install_github('ggcor')
devtools::install_local("D:/R/R-3.6.1/library/ggcor", force = TRUE)
library(ggcor)

spec <-read.csv("bacteria+fungia community.csv",row.names = 1,header = T)
env <-read.csv("phychemical property.csv",row.names = 1,header = T)
mantel <- mantel_test(spec, env, mantel.fun = 'mantel.randtest',
                      spec.dist.method = 'bray', env.dist.method = 'euclidean',
                      spec.select = list("Bacterial community"= 1:6880, "Fungal community" = 6881:8026))
mantel<-mantel_test(spec, env,
                    spec.select = list("Bacterial community" =1:6880,
                                       "Fungal community"  =6881:8026)) %>%
  mutate(cor = cut(r, breaks = c(-Inf, 0.5, Inf),
                 labels = c("<0.5", ">0.5"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">0.05"),
                       right = FALSE),
         lty = cut(r, breaks = c(-Inf, 0, Inf),
                   labels = c("r > 0", "r <= 0")))
mantel
#绘制组合图

p3 <- quickcor(env, cor.test = TRUE, type = 'lower') +
  geom_square(colour=NA)
p3
p3<-p3+anno_link(mantel, mapping = aes(color = p.value, size= cor),label.size = 5, pos = "right")+
  scale_fill_gradient2(low = "#FF7F00", mid = "white",
                       high = "#4DAF4A")+
    theme(legend.position = "right")+
    expand_axis(x = 25) +
  geom_mark(data = get_data(type = 'lower', show.diag = FALSE), sep = '\n', size = 3.0, sig.thres = 0.05)+
  guides(size=guide_legend(title="Mantel's r",override.aes=list(colour="grey35"),order=2),
         colour=guide_legend(title="Mantel's P",override.aes = list(size=3),order=3),
         fill=guide_colorbar(title="Correlations between \n environmental factors",order=1))

p3

ggsave('microbial community correlation with environmental factors.pdf', p3, height = 16, width = 10)
ggsave('microbial community correlation with environmental factors.png', p3, height = 16, width = 10)



