################# This script is for the comparison of relative runtime of K-ABC and AEA
################# setting: total amount of sampled networks for K-ABC is 4 times as many as AEA
################# karate AEA time: 1.573699 mins
################# faux AEA time: 17.97099 mins
faux_time <- readRDS("E:/UCI/Research/ABC for ERGMs/faux_time.rds")
karate_time <- readRDS("E:/UCI/Research/ABC for ERGMs/karate_time.rds")

##################
# correct for the units (in secs)
karate_time_correct <- data.frame(karate_time[,c(3,5)])
karate_time_correct[9:10,1] <- 60 * karate_time_correct[9:10,1]
karate_time_correct$relative <- karate_time_correct$Total.run.time/(1.573699*60)

faux_time_correct <- data.frame(faux_time[,c(3,5)])
faux_time_correct[,1] <- 60* faux_time_correct[,1]
faux_time_correct$relative <- faux_time_correct$Total.run.time/(17.97099*60)

#
run_time_summary <- data.frame( rbind(karate_time_correct, faux_time_correct),
                                dataset = rep(c("Karate Club","Faux Mesa High"), c(nrow(karate_time_correct),
                                                                    nrow(faux_time_correct))) )
names(run_time_summary)

#
ggplot(run_time_summary, aes(x=Number.of.cores,y=relative, color = dataset)) + geom_line() + geom_point() + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + xlab("Number of Cores") + ylab("Relative Time") + guides(fill=guide_legend(title="", label.theme = element_text(size = 20))) + theme(text = element_text(size=20), legend.title = element_blank(),
                                                                                                                                                                                                                     axis.text.x = element_text(angle = 0, hjust=1)) +
        scale_x_continuous(name="Number of cores", limits=c(1, 35), breaks=c(1,5,10,15,20,25,30,35)) + 
        scale_y_continuous(name="Relative Time", breaks = c(0,0.2,0.5,1,2,3)) + 
        geom_hline(yintercept = 1.0, linetype = "dashed") + geom_hline(yintercept = 0.2, linetype="dashed") + theme(
          legend.position = c(0.9, 0.9),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.text = element_text(size = 15)
        )

# general plot
# plot two curves
# par(mfrow=c(1,1))
# plot(faux_time_correct$Number.of.cores, faux_time_correct$relative,
#      col = "red", type = "o", xlab="Number of Cores", ylab="Relative Time", xlim = c(1,35),
#      cex.lab=1.5)
# lines(karate_time_correct$Number.of.cores, karate_time_correct$relative, col = "blue",
#       type = "o")