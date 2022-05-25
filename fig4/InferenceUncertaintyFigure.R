library(copula)
library(ggplot2)
library(gridExtra)
library(FAdist)
library(fitdistrplus)
library(scales)
library(metR)
source("utils.R")

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

set.seed(6)
#u = rCopula(100, normalCopula(0.6, dim=2))
u = rCopula(100, joeCopula(3, dim=2))

v1 = exp(qgamma3(u[,1], 6279, 1/261, -17.5))
v2 = qgev(u[,2], 3.584701e-02, 6.201633e+02, 1.165458e+03)

plot(v1,v2)

# fit same distribution to v1 and v2
LP3.params = gamma3MOM(log(v1))
v1.LP3.fit = fitdist(log(v1), "gamma3", order = c(1,2,3), memp=memp.centered,
                     start=list(shape=LP3.params$alpha, 
                                scale=1/LP3.params$beta,
                                thres=LP3.params$xi),
                     method="mme")
plot(v1.LP3.fit)

gev.params = gevMOM(v2)
v2.gev.fit = fitdist(v2, "gev", start=list(shape=gev.params$kappa, 
                                           scale=gev.params$alpha, 
                                           location=gev.params$xi), 
                     method="mle")
plot(v2.gev.fit)

# fit alternative distribution to v1 and v2
weibull3.params = Weibull3MOM(v1)
v1.weibull3.fit = fitdist(v1, "weibull3", start=list(shape=weibull3.params$kappa, 
                                                     scale=weibull3.params$alpha,
                                                     thres=weibull3.params$xi),
                          method="mle")
plot(v1.weibull3.fit)

LN3.params = lnorm3MOM(v2)
v2.LN3.fit = fitdist(v2, "lnorm3", 
                     start=list(shape=LN3.params$sigma, 
                                scale=LN3.params$mu, 
                                thres=LN3.params$tau),
                      method="mle")
plot(v2.LN3.fit)

# fit same copula to v1 and v2
fitCopula(joeCopula(), u, method="ml")
#Call: fitCopula(copula, data = data, method = "ml")
#Fit based on "maximum likelihood" and 100 2-dimensional observations.
#Copula: joeCopula 
#alpha 
#3.332 
#The maximized loglikelihood is 54.41 
#Optimization converged

# fit alternative copula to v1 and v2
fitCopula(normalCopula(), u, method="ml")
#Call: fitCopula(copula, data = data, method = "ml")
#Fit based on "maximum likelihood" and 100 2-dimensional observations.
#Copula: normalCopula 
#rho.1 
#0.7258 
#The maximized loglikelihood is 38.44 
#Optimization converged


# Make figure
empty = ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())

# make x and y vectors for plotting PDFs of v1 and v2, respectively
df = data.frame(v1,v2)
x = c(seq(0,2000,20))
y = c(seq(0,5000,50))

# Make grid for copula CDF contours
xgrid = meshgrid(x=c(seq(0,2000,20)), y=c(seq(0,5000,50)))$X
ygrid = meshgrid(x=c(seq(0,2000,20)), y=c(seq(0,5000,50)))$Y
coords = data.frame(v1=matrix(xgrid),v2=matrix(ygrid))
xgrid.u = pgamma3(log(xgrid), 6279, 1/261, -17.5)
ygrid.u = pgev(ygrid, 3.584701e-02, 6.201633e+02, 1.165458e+03)
coords.u = data.frame(u1=matrix(xgrid.u),u2=matrix(ygrid.u))
coords$probsTrue = pCopula(as.matrix(coords.u), joeCopula(3, dim=2))
coords$probsJoeFit= pCopula(as.matrix(coords.u), joeCopula(3.332, dim=2))
coords$probsNormalFit= pCopula(as.matrix(coords.u), normalCopula(0.7258, dim=2))

coords.u$probsTrue = coords$probsTrue
coords.u$probsJoeFit= coords$probsJoeFit
coords.u$probsNormalFit= coords$probsNormalFit


override.linetype <- c(1, 2, 3)
hist_v1_all = ggplot() + aes(x=v1) + geom_histogram(aes(y=..density..), colour="gray", fill="black", alpha=0.5) + 
  geom_line(aes(x,dgamma3(log(x), 6279, 1/261, -17.5)/x, color="True LP3 Distribution"),lwd=1.5, lty=1) +
  geom_line(aes(x,dgamma3(log(x), v1.LP3.fit$estimate[1], 
                          v1.LP3.fit$estimate[2], 
                          v1.LP3.fit$estimate[3])/x, color="Fitted LP3 Distribution"),lwd=1.5, lty=2) +
  geom_line(aes(x,dweibull3(x, v1.weibull3.fit$estimate[1], 
                            v1.weibull3.fit$estimate[2], 
                            v1.weibull3.fit$estimate[3]), color="Fitted Weibull3 Distribution"),lwd=1.5, lty=3) +
  scale_color_discrete(breaks=c("True LP3 Distribution","Fitted LP3 Distribution","Fitted Weibull3 Distribution")) +
  theme_bw() + theme(legend.title = element_blank(), legend.position=c(0.7,0.83), text=element_text(size=14)) +
  guides(colour=guide_legend(override.aes=list(linetype=override.linetype))) + scale_linetype(guide=FALSE)

hist_v2_all = ggplot() + aes(x=v2) + geom_histogram(aes(y=..density..), colour="gray", fill="black", alpha=0.5) +
  geom_line(aes(y, dgev(y, 3.584701e-02, 6.201633e+02, 1.165458e+03), 
                color="True GEV Distribution"),lwd=1.5, lty=1) +
  geom_line(aes(y, dgev(y, v2.gev.fit$estimate[1], 
                        v2.gev.fit$estimate[2], 
                        v2.gev.fit$estimate[3]), color="Fitted GEV Distribution"),lwd=1.5, lty=2) + 
  geom_line(aes(y, dlnorm3(y, v2.LN3.fit$estimate[1], 
                           v2.LN3.fit$estimate[2], 
                           v2.LN3.fit$estimate[3]), color="Fitted LN3 Distribution"),lwd=1.5, lty=3) +
  scale_color_discrete(breaks=c("True GEV Distribution","Fitted GEV Distribution","Fitted LN3 Distribution")) +
  theme_bw() + theme(legend.title = element_blank(), legend.position=c(0.75,0.83), text=element_text(size=14)) + coord_flip() +
  guides(colour=guide_legend(override.aes=list(linetype=override.linetype))) + scale_linetype(guide=FALSE)

copula_all = ggplot() + 
  geom_contour(aes(v1, v2, z=probsTrue, colour="True Joe Copula"), coords, lwd=1.5, 
               breaks=c(0.01,0.1,0.3,0.5,0.7,0.9,0.95,0.99), lty=1) + 
  geom_contour(aes(v1, v2, z=probsJoeFit, colour="Fitted Joe Copula"), 
               coords, lwd=1.5, breaks=c(0.01,0.1,0.3,0.5,0.7,0.9,0.95,0.99), lty=2) +
  geom_contour(aes(v1, v2, z=probsNormalFit, colour="Fitted Normal Copula"), 
               coords, lwd=1.5, breaks=c(0.01,0.1,0.3,0.5,0.7,0.9,0.95,0.99), lty=3) +
  geom_point(aes(v1, v2)) + 
  scale_color_discrete(breaks=c("True Joe Copula","Fitted Joe Copula","Fitted Normal Copula")) +
  geom_text_contour(aes(v1, v2, z=probsNormalFit), coords, breaks=c(0.01,0.1,0.3,0.5,0.7,0.9,0.95,0.99),
                    label.placer=label_placer_flattest(), stroke=0.2, skip=0) + theme_bw() +
  theme(legend.title = element_blank(), legend.position=c(0.25,0.83), text=element_text(size=14)) +
  guides(colour=guide_legend(override.aes=list(linetype=override.linetype))) + scale_linetype(guide=FALSE)

grid.arrange(hist_v1_all, empty, copula_all, hist_v2_all, ncol=2, nrow=2)#, 
             #widths=c(1.5,1), heights=c(1,1.5))






# copula fit only
ggplot() + 
  geom_contour(aes(u1, u2, z=probsTrue, colour="True Joe Copula"), coords.u, lwd=1.5, 
               breaks=c(0.01,0.1,0.3,0.5,0.7,0.9,0.95,0.99), lty=1) + 
  geom_contour(aes(u1, u2, z=probsJoeFit, colour="Fitted Joe Copula"), 
               coords.u, lwd=1.5, breaks=c(0.01,0.1,0.3,0.5,0.7,0.9,0.95,0.99), lty=2) +
  geom_contour(aes(u1, u2, z=probsNormalFit, colour="Fitted Normal Copula"), 
               coords.u, lwd=1.5, breaks=c(0.01,0.1,0.3,0.5,0.7,0.9,0.95,0.99), lty=3) +
  geom_point(aes(u[,1], u[,2])) + scale_color_discrete(breaks=c("True Joe Copula","Fitted Joe Copula","Fitted Normal Copula")) +
  geom_text_contour(aes(u1, u2, z=probsNormalFit), coords.u, breaks=c(0.01,0.1,0.3,0.5,0.7,0.9,0.95,0.99),
                    label.placer=label_placer_flattest(), stroke=0.2, skip=0) +
  theme(legend.title = element_blank(), text=element_text(size=14)) +
  guides(colour=guide_legend(override.aes=list(linetype=override.linetype))) + scale_linetype(guide=FALSE)


#alternative plot
hist_v1_true = ggplot() + aes(x=v1) + geom_histogram(aes(y=..density..), colour="gray", fill="black", alpha=0.5) + 
  geom_line(aes(x,dgamma3(log(x), 6279, 1/261, -17.5)/x, color="True LP3\nDistribution"),lwd=1.5) +
  theme(legend.title = element_blank())

hist_v2_true = ggplot() + aes(x=v2) + geom_histogram(aes(y=..density..), colour="gray", fill="black", alpha=0.5) +
  geom_line(aes(y, dgev(y, 3.584701e-02, 6.201633e+02, 1.165458e+03), 
                color="True GEV\nDistribution"),lwd=1.5) +
  theme(legend.title = element_blank()) + coord_flip() + scale_x_reverse()

copula_true = ggplot() + geom_contour_filled(aes(v1, v2, z=probsTrue), coords) + 
  geom_point(aes(v1, v2), color="white") + labs(title="True Joe Copula")

hist_v2_param_estim = ggplot() + aes(x=v2) + geom_histogram(aes(y=..density..), colour="gray", fill="black", alpha=0.5) +
  geom_line(aes(y, dgev(y, 3.584701e-02, 6.201633e+02, 1.165458e+03), 
                color="True GEV\nDistribution"),lwd=1.5) +
  geom_line(aes(y, dgev(y, v2.gev.fit$estimate[1], 
                        v2.gev.fit$estimate[2], 
                        v2.gev.fit$estimate[3]), color="Fitted GEV\nDistribution"),lwd=1.5) + 
  theme(legend.title = element_blank()) + coord_flip() + scale_x_reverse()

copula_param_estim = ggplot() + geom_contour_filled(aes(v1, v2, z=probsJoeFit), coords) + 
  geom_point(aes(v1, v2), color="white") + labs(title="Fitted Joe Copula")

copula_struct_estim = ggplot() + geom_contour_filled(aes(v1, v2, z=probsNormalFit), coords) + 
  geom_point(aes(v1, v2), color="white") + labs(title="Fitted Gaussian Copula")

hist_v2_struct_estim = ggplot() + aes(x=v2) + geom_histogram(aes(y=..density..), colour="gray", fill="black", alpha=0.5) +
  geom_line(aes(y, dgev(y, 3.584701e-02, 6.201633e+02, 1.165458e+03), 
                color="True GEV Distribution"),lwd=1.5) +
  geom_line(aes(y, dlnorm3(y, v2.LN3.fit$estimate[1], 
                           v2.LN3.fit$estimate[2], 
                           v2.LN3.fit$estimate[3]), color="Fitted LN3\nDistribution"),lwd=1.5) + 
  theme(legend.title = element_blank()) + coord_flip()

hist_v1_param_estim = ggplot() + aes(x=v1) + geom_histogram(aes(y=..density..), colour="gray", fill="black", alpha=0.5) + 
  geom_line(aes(x,dgamma3(log(x), 6279, 1/261, -17.5)/x, color="True LP3\nDistribution"),lwd=1.5) +
  geom_line(aes(x,dgamma3(log(x), v1.LP3.fit$estimate[1], 
                          v1.LP3.fit$estimate[2], 
                          v1.LP3.fit$estimate[3])/x, color="Fitted LP3\nDistribution"),lwd=1.5) +
  theme(legend.title = element_blank()) + scale_y_reverse()

hist_v1_struct_estim = ggplot() + aes(x=v1) + geom_histogram(aes(y=..density..), colour="gray", fill="black", alpha=0.5) + 
  geom_line(aes(x,dgamma3(log(x), 6279, 1/261, -17.5)/x, color="True LP3\nDistribution"),lwd=1.5) +
  geom_line(aes(x,dweibull3(x, v1.weibull3.fit$estimate[1], 
                            v1.weibull3.fit$estimate[2], 
                            v1.weibull3.fit$estimate[3]), color="Fitted Weibull3\nDistribution"),lwd=1.5) +
  theme(legend.title = element_blank()) + scale_y_reverse()

grid.arrange(empty, hist_v1_true, empty, empty, 
             hist_v2_true, copula_true, empty, empty, 
             hist_v2_param_estim, copula_param_estim, copula_struct_estim, hist_v2_struct_estim,
             empty, hist_v1_param_estim, hist_v1_struct_estim, empty, ncol=4, nrow=4, 
             widths=c(1,1.5,1.5,1), heights=c(1,1.5,1.5,1))