
plot_gam_diff <- function(plotdf){
  p<-ggplot(plotdf) + geom_ribbon(aes(x=x, ymin=lb.diff, ymax=ub.diff), alpha=0.5) +
    geom_path(aes(x=x, y=lb.diff), color="blue")+
    geom_path(aes(x=x, y=ub.diff), color="red")+
    geom_path(aes(x=x, y=point.diff), color="black") +
    geom_vline(aes(xintercept=q1)) +
    geom_vline(aes(xintercept=q3)) +
    xlab("Exposure") +
    ylab("GAM-estimated differences from 25th percentile of exposure")
  return(p)
}
