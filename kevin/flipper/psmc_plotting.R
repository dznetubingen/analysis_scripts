# Plotting the result of PSMC analysis using R.
# By Shenglin Liu, Apr 3, 2016.


##-------Rescale the ith iteration result of PSMC, and make ready for plotting
# file: result file from PSMC
# i.iteration: the ith iteration
# mu: mutation rate
# s: bin size
# g: years per generation

psmc.result<-function(file,i.iteration=25,mu=1e-8,s=100,g=1)
{
  X<-scan(file=file,what="",sep="\n",quiet=TRUE)
  
  START<-grep("^RD",X)
  END<-grep("^//",X)
  
  X<-X[START[i.iteration+1]:END[i.iteration+1]]
  
  TR<-grep("^TR",X,value=TRUE)
  RS<-grep("^RS",X,value=TRUE)
  
  write(TR,"temp.psmc.result")
  theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
  N0<-theta0/4/mu/s
  
  write(RS,"temp.psmc.result")
  a<-read.table("temp.psmc.result")
  Generation<-as.numeric(2*N0*a[,3])
  Ne<-as.numeric(N0*a[,4])
  
  file.remove("temp.psmc.result")
  
  n.points<-length(Ne)
  YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
              Generation[n.points])*g
  Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
        Ne[n.points])
  
  data.frame(YearsAgo,Ne)
}


##-------Turn "ms" commandline into history, and make ready for plotting
# N0: diploid population size

history<-function(command,N0,g=1)
{
  X<-unlist(strsplit(command," "))
  X<-matrix(X,3,length(X)/3)
  Generation<-as.numeric(X[2,])*4*N0
  Generation<-c(10,Generation)
  Ne<-as.numeric(X[3,])*N0
  Ne<-c(N0,Ne)
  
  n.points<-length(Ne)
  Generation<-c(Generation,2*Generation[n.points])
  Ne<-c(Ne,Ne[n.points])
  
  n.points<-length(Ne)
  YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
              Generation[n.points])*g
  Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
        Ne[n.points])
  
  data.frame(YearsAgo,Ne)
}


##-------Plot population by population; real history can be added in ms commandline
# keywords: regular expression; beginning format of the psmc files' names
# labels: names of the output plots
# command: e.g., "-eN 0.01 0.1 -eN 0.06 1 -eN 0.2 0.5 -eN 1 1 -eN 2 2"

plotPsmc.popwise<-function(keywords, labels,
                           command=NA, N0=25000,
                           save.as="png", height=7, width=12,
                           i.iteration=25, mu=1e-8, s=100, g=1,
                           ylim=c(0,100000), xlim=c(200,500000),
                           col=rep("red",length(keywords)), col.hist="grey")
{
  n<-length(keywords)
  files<-grep("\\.psmc$",dir(),value=TRUE)
  dev.new(height=height,width=width)
  for(i in 1:n)
  {
    subfiles<-grep(paste("^",keywords[i],sep=""),files,value=TRUE)
    n.sub<-length(subfiles)
    plot(1,1,
         ylim=ylim,xlim=xlim,
         log="x",type="n",
         ylab="Ne",xlab="Generations ago")
    for(i.sub in 1:n.sub)
    {
      lines(psmc.result(subfiles[i.sub],i.iteration,mu,s,g),
            type="l",col=col[i],lwd=1)
    }
    if(!is.na(command))
    {
      lines(history(command,N0,g),
            type="l",col=col.hist,lwd=2)
    }
    savePlot(filename=paste(labels[i],save.as,sep="."),type=save.as)
  }
  dev.off()
}


##-------Plot all populations together

plotPsmc.allPops<-function(keywords, label, legend.names,
                           save.as="png", height=7, width=12,
                           i.iteration=25, mu=3.7e-8, s=100, g=2,
                           ylim=c(0,60000), xlim=c(200,500000),
                           col=rainbow(length(keywords)))
{
  n<-length(keywords)
  files<-grep("\\.psmc$",dir(),value=TRUE)
  dev.new(height=height,width=width)
  plot(1,1,
       ylim=ylim,xlim=xlim,
       log="x",type="n",
       main=label[1],ylab="Ne",xlab="Years ago")
  for(i in 1:n)
  {
    subfiles<-grep(paste("^",keywords[i],sep=""),files,value=TRUE)
    n.sub<-length(subfiles)
    for(i.sub in 1:n.sub)
    {
      lines(psmc.result(subfiles[i.sub],i.iteration,mu,s,g),
            type="l",col=col[i],lwd=2)
    }
  }
  legend("topright",legend=legend.names,col=col,lty=1,lwd=2)
  savePlot(filename=paste(label[1],save.as,sep="."),type=save.as)
  dev.off()
}