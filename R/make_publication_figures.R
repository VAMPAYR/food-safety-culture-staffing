# Publication figures from saved complete-case output; no models are fitted.
# Run from the repository directory: Rscript R/make_publication_figures.R
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/")
input <- file.path(root, "results", "complete_case")
output <- file.path(root, "figures")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

staff <- read.csv(file.path(input, "main_staffing_estimates.csv"), stringsAsFactors = FALSE)
staff <- staff[staff$term %in% c("B_c", "A_c"), ]
staff <- staff[match(c("6B_c", "6A_c", "4B_c", "4A_c"), paste0(staff$outcome, staff$term)), ]
stopifnot(nrow(staff) == 4, all(staff$valid), identical(staff$n, c(539L,539L,541L,541L)),
          all(staff$exp_lower < staff$exp_estimate), all(staff$exp_upper > staff$exp_estimate),
          all(abs(exp(staff$estimate) - staff$exp_estimate) < 1e-12))
staff$label <- rep(c("Busy-period coverage", "Absence coverage"), 2)
write.csv(staff, file.path(output, "figure1_plot_data.csv"), row.names = FALSE)

draw_forest <- function() {
  layout(matrix(1:3, nrow=1), widths=c(2.05,2.35,2.1))
  yy <- c(3.8,3.1,1.6,.9)
  par(mar=c(3.7,.2,1.1,0), family="serif", cex=.95)
  plot.new(); plot.window(xlim=c(0,1), ylim=c(.2,4.8), xaxs="i", yaxs="i")
  text(.02,4.55,"Perceived sick-worker exclusion",adj=0,font=2,cex=.91)
  text(.02,4.22,"Q6, n = 539",adj=0,cex=.9)
  text(.02,2.35,"Perceived handwashing",adj=0,font=2,cex=.91)
  text(.02,2.02,"Q4, n = 541",adj=0,cex=.9)
  text(.08,yy,staff$label,adj=0,cex=.94)
  par(mar=c(3.7,0,1.1,.1))
  plot.new(); plot.window(xlim=c(.6,2.5), ylim=c(.2,4.8), log="x", xaxs="i", yaxs="i")
  abline(v=1,lty=2,col="#888888")
  segments(staff$exp_lower,yy,staff$exp_upper,yy,lwd=1.5,col="#245A81")
  segments(staff$exp_lower,yy-.07,staff$exp_lower,yy+.07,lwd=1.2,col="#245A81")
  segments(staff$exp_upper,yy-.07,staff$exp_upper,yy+.07,lwd=1.2,col="#245A81")
  points(staff$exp_estimate,yy,pch=15,cex=1.05,col="#245A81")
  axis(1,at=c(.6,.8,1,1.5,2,2.5),labels=c("0.6","0.8","1.0","1.5","2.0","2.5"),cex.axis=.9)
  mtext("Adjusted odds ratio (log scale)",1,line=2.35,cex=.9)
  par(mar=c(3.7,.1,1.1,.15))
  plot.new(); plot.window(xlim=c(0,1), ylim=c(.2,4.8), xaxs="i", yaxs="i")
  text(.5,4.55,"Odds ratio (95% CI)",font=2,cex=.95)
  text(.5,yy,sprintf("%.2f (%.2f to %.2f)",staff$exp_estimate,staff$exp_lower,staff$exp_upper),cex=.95)
}

cells <- read.csv(file.path(input,"rq4_six_cells.csv"),stringsAsFactors=FALSE)
order_keys <- paste(rep(c(FALSE,TRUE),each=3),rep(c("inadequate","neutral","adequate"),2))
groups <- lapply(order_keys,function(key) {
  g <- cells[paste(cells$leadership_high,cells$absence_category)==key,]
  g <- g[order(g$outcome_category),]
  stopifnot(nrow(g)==4, identical(g$outcome_category,1:4),length(unique(g$n))==1,
            sum(g$outcome_n)==g$n[1], abs(sum(g$outcome_percent)-100)<1e-10,
            all(abs(g$outcome_percent -100*g$outcome_n/g$n)<1e-10))
  g
})
plot_cells <- do.call(rbind,groups)
stopifnot(sum(vapply(groups,function(g)g$n[1],numeric(1)))==558)
write.csv(plot_cells,file.path(output,"figure2_plot_data.csv"),row.names=FALSE)

draw_bars <- function() {
  layout(matrix(c(1,2),nrow=1),widths=c(1.8,4.7))
  ypos <- c(6.8,5.8,4.8,3.2,2.2,1.2)
  palette <- c("#CA8C34","#D4D4D4","#6CA9C9","#245A81")
  par(mar=c(7,.15,1.2,0),family="serif",cex=.95)
  plot.new();plot.window(xlim=c(0,1),ylim=c(.5,8),xaxs="i",yaxs="i")
  text(.02,7.7,"Leadership below 4",adj=0,font=2,cex=.98)
  text(.02,4.1,"Leadership at least 4",adj=0,font=2,cex=.98)
  for(i in seq_along(groups)) {
    g <- groups[[i]]
    coverage <- paste0(toupper(substr(g$absence_category[1],1,1)),substr(g$absence_category[1],2,99))
    text(.04,ypos[i],sprintf("%s (n = %d)",coverage,g$n[1]),adj=0,cex=.93)
  }
  par(mar=c(7,0,1.2,1.2))
  plot.new();plot.window(xlim=c(0,100),ylim=c(.5,8),xaxs="i",yaxs="i")
  for(i in seq_along(groups)) {
    p <- groups[[i]]$outcome_percent
    bounds <- c(0,cumsum(p))
    for(j in 1:4) {
      rect(bounds[j],ypos[i]-.33,bounds[j+1],ypos[i]+.33,col=palette[j],border="white",lwd=.6)
      if(p[j]>=7) text(mean(bounds[j:(j+1)]),ypos[i],sprintf("%.1f",p[j]),
                       cex=.84,col=if(j==4)"white" else "#151515")
    }
  }
  axis(1,at=seq(0,100,20),cex.axis=.9)
  mtext("Workers within each group (%)",1,line=2.1,cex=.94)
  par(xpd=NA)
  legend(0,-.75,legend=c("Disagreement","Neutral","Agreement","Strong agreement"),
         fill=palette,border=NA,bty="n",ncol=2,cex=.91,x.intersp=.6,y.intersp=1.1)
}

for(format in c("png","pdf")) {
  for(fig in 1:2) {
    stem <- if(fig==1)"figure1_staffing" else "figure2_leadership"
    dest <- file.path(output,paste0(stem,".",format))
    if(format=="png") png(dest,width=6.5,height=if(fig==1)3.7 else 5.3,units="in",res=300,pointsize=10,type="cairo")
    else pdf(dest,width=6.5,height=if(fig==1)3.7 else 5.3,pointsize=10,family="Times",useDingbats=FALSE)
    if(fig==1)draw_forest() else draw_bars()
    dev.off()
  }
}
writeLines(c(
  "# Publication figure validation",
  "",
  "Both figures use saved complete-case CSV output only. No model was refitted.",
  "Figure 1 retains exactly four valid B_c/A_c rows; Q6 uses n=539 and Q4 n=541.",
  "Exponentiated coefficients and interval order were checked to numerical tolerance.",
  "Figure 2 retains all 24 cells, with ordered categories 1 through 4 within six groups.",
  "Every group count equals the sum of its four counts; group denominators total 558.",
  "Every plotted percentage matches 100 times its count divided by its group denominator.",
  "Percentages below 7 are not printed inside segments; all values remain in plot data and Table 3.",
  "Source file MD5 values:",
  paste(names(tools::md5sum(file.path(input,c("main_staffing_estimates.csv","rq4_six_cells.csv")))),
        tools::md5sum(file.path(input,c("main_staffing_estimates.csv","rq4_six_cells.csv"))))) ,
  file.path(output,"figure_validation.md"))
cat("Publication figures and plot data saved to",output,"\n")
