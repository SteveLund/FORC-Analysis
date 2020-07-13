require(akima)
require(fields)
require(dplyr)
require(stats)
require(lmenssp)
require(locfit)
require(rgl)
require(ggplot2)


### I
pres_plots<-function(dat1,dat2=NULL,method="",exclude_rng=NULL){
  ### This function is used to plot results of various fitting methods
  if(is.null(dat2))dat2<-dat1
  CA<-CL<-1.7
  par(mfcol=c(2,3))
  par(mai=c(.7,.7,1,1.2))
  nx<-ny<-300
  
  zz<-akima::interp(x=dat1$H_r,y=dat1$H,z=dat1$Moment,nx=nx,ny=ny)
  image.plot(zz$z,x=zz$x,y=zz$y,xlab=expression(mu[0]*H[R]),ylab=expression(mu[0]*H),main="Moment",zlim=range(ave.dat$Moment),cex.lab=CL,cex.axis=CA)
  
  zz<-akima::interp(x=dat1$H_r,y=dat1$H,z=dat1$pred,nx=nx,ny=ny)
  image.plot(zz$z,x=zz$x,y=zz$y,
             xlab=expression(mu[0]*H[R]*" (mT)"),ylab=expression(mu[0]*H*" (mT)"),
             main=paste(method,"
Smoothed Moment"),
             xlim=range(ave.dat$H_r),ylim=range(ave.dat$H),zlim=range(ave.dat$Moment),
             cex.lab=CL,cex.axis=CA)
  dat1$rnd.H_r<-2.5*round(dat1$H_r/2.5)
  
  h<-sort(unique(dat1$rnd.H))
  H_r<-sort(unique(dat1$rnd.H_r))
  zz<-rep(NA,length(h)*length(H_r))
  names(zz)<-paste(rep(H_r,each=length(h)),rep(h,length(H_r)))
  zz[paste(dat1$rnd.H_r,dat1$rnd.H)]<-(dat1$Moment-dat1$pred)/dat1$Moment_SE
  zmat<-matrix(zz,length(h),length(H_r))
  if(!is.null(exclude_rng))zmat[h>exclude_rng[1]&h<exclude_rng[2],]<-NA
  image.plot(t(zmat),x=H_r,y=h,
             xlab=expression(mu[0]*H[R]*" (mT)"),
             ylab=expression(mu[0]*H*" (mT)"),
             main=paste(method,"Standardized Residual"),
             xlim=range(H_r),
             ylim=range(h),
             zlim=quantile(zmat,c(.003,.997),na.rm=TRUE),
             cex.lab=CL,cex.axis=CA)
  
  par(mai=c(.7,.7,1,.2))
  norm.res<-(dat1$Moment-dat1$pred)/dat1$Moment_SE
  d<-density(norm.res)
  Ylim=range(d$y)
  Xlim=range(c(-3,3,d$x[d$y>(Ylim[2]/100)]))
  plot(density(norm.res,from=Xlim[1]-.05*diff(Xlim),to=Xlim[2]+.05*diff(Xlim)),
       xlim=Xlim,
       ylim=Ylim,
       lwd=2,
       lty=2,
       xlab="Standardized Residual",
       main="",
       cex.lab=CL,cex.axis=CA)
  grid()
  xx<-seq(-10,10,.01)
  #lines(xx,dt(xx,24),col=4,lwd=2)
  #    text(sum(c(.15,.85)*Xlim),(sum(c(.05,.95)*(Ylim))),expression(sigma),cex=CA)
  #text(sum(c(.15,.85)*Xlim),(sum(c(.15,.85)*(Ylim))),1.04,col=4,cex=CA)
  #    text(sum(c(.15,.85)*Xlim),(sum(c(.15,.85)*(Ylim))),signif(sd(norm.res),3),cex=CA)
  if("cv.pred"%in%names(dat1)){
    norm.cv.res<-(dat1$Moment-dat1$cv.pred)/dat1$Moment_SE
    lines(density(norm.cv.res,from=Xlim[1]-.05*diff(Xlim),to=Xlim[2]+.05*diff(Xlim)),lwd=2)
    legend("topleft",title= "CV (SD)",
           legend=c(paste0("Yes (",signif(sd(norm.cv.res),3),")"),
                    paste0("No (",signif(sd(norm.res),3),")")),
           lty=c(1,2),lwd=2,cex=1.3,bg="white")
  } 
  par(mai=c(.7,.7,1,1.2))
  

  dat3<-dat2[!is.na(dat2$dmomdhdhr),]
  zz<-akima::interp(x=dat3$H_r,y=dat3$H,z=-dat3$dmomdhdhr,nx=nx,ny=ny)
  image.plot(zz$z,x=zz$x,y=zz$y,xlab=expression(mu[0]*H[R]*" (mT)"),
             ylab=expression(mu[0]*H*" (mT)"),
             main=bquote(.(method)~":" ~ rho),
             xlim=range(dat3$H_r),ylim=range(dat3$H),
             zlim=c(-1.5e-5,2.75e-5),
             cex.lab=CL,cex.axis=CA)
  
  zz<-akima::interp(x=dat3$H_b,y=dat3$H_c,z=-dat3$dmomdhdhr,nx=nx,ny=ny)
  image.plot(zz$z,x=zz$x,y=zz$y,xlab=expression(mu[0]*H[b]*" (mT)"),
             ylab=expression(mu[0]*H[c]*" (mT)"),
             main=bquote(.(method)~":" ~ rho),
             xlim=range(dat3$H_b),ylim=range(dat3$H_c),
             zlim=c(-1.5e-5,2.75e-5),
             cex.lab=CL,cex.axis=CA)
  
  par(mfrow=c(1,3))
  par(mai=c(.7,.7,.6,.2))
  image.plot(zz$z,x=zz$x,y=zz$y,
             xlab=expression(mu[0]*H[b]*" (mT)"),ylab=expression(mu[0]*H[c]*" (mT)"),main=bquote(.(method)~":" ~ rho),
             xlim=range(dat3$H_b),ylim=range(dat3$H_c),zlim=c(-1.5e-5,2.75e-5),
             cex.lab=CL,cex.axis=CA)
  plot(zz$x,rowSums(zz$z,na.rm=TRUE)*mean(diff(zz$x)),
       xlab=expression(mu[0]*H[b]*" (mT)"),ylab="Density",main=paste(method),type="l",
       xlim=range(dat3$H_b),ylim=c(-.00015,.0017),cex.lab=CL,cex.axis=CA)
  plot(zz$y,colSums(zz$z,na.rm=TRUE)*mean(diff(zz$y)),xlab=expression(mu[0]*H[c]*" (mT)"),
       ylab="Density",main=paste(method),type="l",
       xlim=range(dat3$H_c),ylim=c(-2e-5,2e-4),cex.lab=CL,cex.axis=CA)
  
  mat1<-cbind(zz$x,rowSums(zz$z,na.rm=TRUE)*mean(diff(zz$y)))
  colnames(mat1)<-c("H_b","Density")
  
  mat2<-cbind(zz$y,colSums(zz$z,na.rm=TRUE)*mean(diff(zz$x)))
  colnames(mat2)<-c("H_c","Density")
  
  write.csv(mat1,
            file=paste0(File.name,"/Marginal Densities/",gsub(":","through",method),"_H_b.csv"),
            row.names=FALSE)
  write.csv(mat2,
            file=paste0(File.name,"/Marginal Densities/",gsub(":","through",method),"_H_c.csv"),
            row.names=FALSE)
}

res_plots<-function(dat1,dat2=NULL,method="",exclude_rng=NULL){
  ### This function is used to plot results of various fitting methods
  
  if(is.null(dat2))dat2<-dat1
  CA<-CL<-1.7
  par(mfcol=c(2,3))
  par(mai=c(.7,.7,1,1.2))
  nx<-ny<-300
  dat1$rnd.H_r<-2.5*round(dat1$H_r/2.5)
  h<-sort(unique(dat1$rnd.H))
  H_r<-sort(unique(dat1$rnd.H_r))
  zz<-rep(NA,length(h)*length(H_r))
  names(zz)<-paste(rep(H_r,each=length(h)),rep(h,length(H_r)))
  zz[paste(dat1$rnd.H_r,dat1$rnd.H)]<-(dat1$Moment-dat1$pred)/dat1$Moment_SE
  zmat<-matrix(zz,length(h),length(H_r))
  if(!is.null(exclude_rng))zmat[h>exclude_rng[1]&h<exclude_rng[2],]<-NA
  image.plot(t(zmat),x=H_r,y=h,
             xlab=expression(mu[0]*H[R]*" (mT)"),
             ylab=expression(mu[0]*H*" (mT)"),
             main=paste(method,"Standardized Residual"),
             xlim=range(H_r),
             ylim=range(h),
             zlim=quantile(zmat,c(.003,.997),na.rm=TRUE),
             cex.lab=CL,cex.axis=CA)
  
  par(mai=c(.7,.7,1,.2))
  norm.res<-(dat1$Moment-dat1$pred)/dat1$Moment_SE
  d<-density(norm.res)
  Ylim=range(d$y)
  Xlim=range(c(-3,3,d$x[d$y>(Ylim[2]/100)]))
  plot(density(norm.res,from=Xlim[1]-.05*diff(Xlim),to=Xlim[2]+.05*diff(Xlim)),
       xlim=Xlim,
       ylim=Ylim,
       lwd=2,
       lty=2,
       xlab="Standardized Residual",
       main="",
       cex.lab=CL,cex.axis=CA)
  grid()
  xx<-seq(-10,10,.01)
  #lines(xx,dt(xx,24),col=4,lwd=2)
  #    text(sum(c(.15,.85)*Xlim),(sum(c(.05,.95)*(Ylim))),expression(sigma),cex=CA)
  #text(sum(c(.15,.85)*Xlim),(sum(c(.15,.85)*(Ylim))),1.04,col=4,cex=CA)
  #    text(sum(c(.15,.85)*Xlim),(sum(c(.15,.85)*(Ylim))),signif(sd(norm.res),3),cex=CA)
  if("cv.pred"%in%names(dat1)){
    norm.cv.res<-(dat1$Moment-dat1$cv.pred)/dat1$Moment_SE
    lines(density(norm.cv.res,from=Xlim[1]-.05*diff(Xlim),to=Xlim[2]+.05*diff(Xlim)),lwd=2)
    legend("topleft",title= "CV (SD)",
           legend=c(paste0("Yes (",signif(sd(norm.cv.res),3),")"),
                    paste0("No (",signif(sd(norm.res),3),")")),
           lty=c(1,2),lwd=2,cex=1.3,bg="white")
  } else{
    legend("topleft",
           legend=c(paste0("No CV (",signif(sd(norm.res),3),")")),
           lty=c(2),lwd=2,cex=1.3,bg="white")
  }
  par(mai=c(.7,.7,1,1.2))
  

  dat3<-dat2[!is.na(dat2$dmomdhdhr),]
}

fit.1d<-function(dat,  CV=FALSE){
  ### This function is used to fit sequential cubic splines
  ### The input arguments are:
  ###     `dat' - the data frame for the FORC measurements being analyzed
  ###     `CV' - a binary indicator for whether or not cross-validation should be utilized`
  H<-"H"
  Moment<-"Moment"
  path<-dat$path
  ### Try fitting a single path using cubic splines
  dat$weighted.d2mom.dH.dHr<-dat$weighted.dmom.dH.pred<-dat$weighted.dmom.dH.cv.pred<-
    dat$d2mom.dH.dHr<-dat$dmom.dH.pred<-dat$dmom.dH.cv.pred<-dat$dmom.dH.resid<-
    dat$weighted.dmom.dH<-dat$weighted.cv.pred<-dat$weighted.pred<-
    dat$dmom.dH<-dat$cv.pred<-dat$pred<-rep(NA,nrow(dat))
  
  ### Estimate dmom.dH
  uni_hr<-unique(dat$H_r)
  weighted.spline.fit<-spline.fit<-list()
  spline.der<-function(spline.fit,x,delta.x){
    hi<-predict(spline.fit,x=x+delta.x)
    low<-predict(spline.fit,x=x-delta.x)
    .5*(hi-low)/delta.x
  }
  
  
  for(hr_ind in 1:length(uni_hr)){
    ind<-dat$H_r==uni_hr[hr_ind]
    invisible(spline.fit[[hr_ind]]<-sreg(dat$H[ind],dat$Moment[ind]))
    dat$dmom.dH[ind]<-spline.der(spline.fit[[hr_ind]],dat$rnd.H[ind],.0001)
    dat$pred[ind]<-spline.fit[[hr_ind]]$fitted.values
    k<-sample(rep(1:10,1000)[1:sum(ind)])
    if(CV){
      for(jj in 1:10){
        invisible(t.fit<-sreg(dat$H[ind][k!=jj],dat$Moment[ind][k!=jj]))
        dat$cv.pred[ind][k==jj]<-predict(t.fit,x=dat$H[ind][k==jj])
      }
    }
  }
  
  ### Evaluate d2mom.dHdHr
  ####Fit spline dmom.dH across H_r initial time
  spline.fit2<-spline.fit3<-list()
  h_table<-table(dat$rnd.H)
  for(ii in 1:sum(h_table>30)){
    H<-as.numeric(names(h_table)[ h_table>30])[ii]
    ind<-dat$rnd.H==H
    invisible(spline.fit2[[ii]]<-sreg(dat$H_r[ind],dat$dmom.dH[ind]))
    dat$dmom.dH.resid[ind]<-spline.fit2[[ii]]$residuals
    dat$dmom.dH.pred[ind]<-spline.fit2[[ii]]$fitted.values
    dat$d2mom.dH.dHr[ind]<-spline.der(spline.fit2[[ii]],dat$rnd.H_r[ind],.001)
    k<-sample(rep(1:10,1000)[1:sum(ind)])
    for(jj in 1:10){
      invisible(t.fit<-sreg(dat$H_r[ind][k!=jj],dat$dmom.dH[ind][k!=jj]))
      dat$dmom.dH.cv.pred[ind][k==jj]<-predict(t.fit,x=dat$H_r[ind][k==jj])
    }
  }
  
  if(nrow(dat)>5000){
    require(locfit)
    
    ## Construct weights to use when fitting Moment vs H
    dat$res<-dat$Moment-dat$pred
    dat$abs.res<-abs(dat$res)
    smooth.abs.res<-locfit(abs.res~lp(H,H_r,deg=3,nn=.001),data=dat)
    dat$weights<-1/fitted(smooth.abs.res)^2
    
    for(hr_ind in 1:length(uni_hr)){
      ind<-dat$H_r==uni_hr[hr_ind]
      invisible(weighted.spline.fit[[hr_ind]]<-sreg(dat$H[ind],dat$Moment[ind],weights = dat$weights[ind]))
      dat$weighted.dmom.dH[ind]<-spline.der(weighted.spline.fit[[hr_ind]],dat$rnd.H[ind],.0001)
      dat$weighted.pred[ind]<-weighted.spline.fit[[hr_ind]]$fitted.values
      k<-sample(rep(1:10,1000)[1:sum(ind)])
      if(CV){
        for(jj in 1:10){
          invisible(t.fit<-sreg(dat$H[ind][k!=jj],dat$Moment[ind][k!=jj],weights = dat$weights[ind][k!=jj]))
          dat$weighted.cv.pred[ind][k==jj]<-predict(t.fit,x=dat$H[ind][k==jj])
        }
      }
    }
    
    ### Refit splines to dmom.dH using weighted estimates
    for(ii in 1:sum(h_table>30)){
      H<-as.numeric(names(h_table)[ h_table>30])[ii]
      ind<-dat$rnd.H==H
      invisible(spline.fit3[[ii]]<-sreg(dat$H_r[ind],dat$weighted.dmom.dH[ind]))
      dat$weighted.dmom.dH.pred[ind]<-spline.fit3[[ii]]$fitted.values
    }
    
    ## Construct weights to use when fitting dMoment/dH vs H_r
    dat$res2<-dat$weighted.dmom.dH-dat$weighted.dmom.dH.pred
    dat$abs.res2<-abs(dat$res2)
    smooth.abs.res2<-locfit(abs.res2~lp(H,H_r,deg=3,nn=.001),data=dat)
    dat$weights2[!is.na(dat$weighted.dmom.dH.pred)]<-1/fitted(smooth.abs.res2)^2
    dat$weights2<- dat$weights2/mean( dat$weights2,na.rm=TRUE)
    for(ii in 1:sum(h_table>30)){
      H<-as.numeric(names(h_table)[ h_table>30])[ii]
      ind<-dat$rnd.H==H
      invisible(spline.fit3[[ii]]<-sreg(dat$H_r[ind],dat$weighted.dmom.dH[ind],weights=dat$weights2[ind]))
      dat$weighted.dmom.dH.pred[ind]<-spline.fit3[[ii]]$fitted.values
      dat$weighted.d2mom.dH.dHr[ind]<-spline.der(spline.fit3[[ii]],dat$rnd.H_r[ind],.001)
      k<-sample(rep(1:10,1000)[1:sum(ind)])
      if(CV){
        for(jj in 1:10){
          invisible(t.fit<-sreg(dat$H_r[ind][k!=jj],dat$weighted.dmom.dH[ind][k!=jj],weights=dat$weights2[ind][k!=jj]))
          dat$weighted.dmom.dH.cv.pred[ind][k==jj]<-predict(t.fit,x=dat$H_r[ind][k==jj])
        }
      }
    }
  }
  dat
}


for(File.name in c("FORC012514", "FORC011714","FORC02--14.fixed","FORC02--14.free")){
  setwd("C:/Users/spl/Documents/Magnetic Nanoparticle with Cindi Dennis/")
  orig.file<-File.name
  if (!file.exists(File.name)){
    dir.create(File.name)
    dir.create(file.path(File.name,"Marginal Densities"))
  }
  File.name.mod<-""
  ## Does magfield decrease within minor loops?
  BACKWARD<-0
  MagField<-"Magnetic Field (Oe)"
  Moment<-"Moment (emu)"
  MomSE<-"M. Std. Err. (emu)"
  if(File.name=="FORC02--14.fixed"){
    Moment<-"DC Moment Fixed Ctr (emu)"
    MomSE<-"DC Moment Err Fixed Ctr(emu)"
    File.name<-"FORC02--14"
  }
  if(File.name=="FORC02--14.free"){
    Moment<-"DC Moment Free Ctr (emu)"
    MomSE<-"DC Moment Err Free Ctr(emu)"
    File.name<-"FORC02--14"
  }
  
  #File.name="FORC022814"; BACKWARD<-1
  
  ### Load data
  source("LoadFORCdata2.R")
  
  ### Look at differences between consecutive reversal curves as a way to assess whether
  ### the data show evidence of an random additive offset to each reversal curve.
    dat2<-NULL
  for(p in 2:max(dat[,"path"]-1)){
    t.dat<-dat[dat[,"path"]==(p-1),]
    t.dat2<-dat[dat[,"path"]==(p+1),]
    invisible(f1<-sreg(t.dat[,"H"],t.dat[,"Moment"]))
    invisible(f2<-sreg(t.dat2[,"H"],t.dat2[,"Moment"]))
    p1<-predict(f1,x=dat[dat[,"path"]==p,"H"][-1])
    p2<-predict(f2,x=dat[dat[,"path"]==p,"H"][-1])
    w=(dat[dat[,"path"]==p,"H_r"][1]-dat[dat[,"path"]==(p-1),"H_r"][1])/
      (dat[dat[,"path"]==(p+1),"H_r"][1]-dat[dat[,"path"]==(p-1),"H_r"][1])
    zz<-  dat[dat[,"path"]==p,"Moment"][-1]-c(p2*w+p1*(1-w))
    dat2<-rbind(dat2,cbind(dat[dat[,"path"]==p,][-1,],"NeighborDiff"=zz))
  }
 
  
  par(mai=c(1,1,.2,.2))
  plot(dat2[,"H_r"],#+2.5/500*dat2[,"H"],
       dat2[,"NeighborDiff"],
       col=dat2[,"path"],
       xlab=expression(mu[0]*H[R]*" (mT)"),
       ylab=expression("Moment Residual  (x10"^"-3"*"A-m"^2*"/kg)"))


  
  layout(matrix(1:6,2,3))
  pdf(file.path(File.name,paste0(orig.file,"_Selected_Fitting_results_06132020.pdf")),height=4.75,width=12)
  
  remove_outlier_list<-0
  if(File.name=="FORC02--14")  remove_outlier_list<-0:1  ### This file contained some extreme outliers
  for(rm_out in remove_outlier_list){
    if(rm_out) ave.dat<-ave.dat[ave.dat$Moment<.13,]
    
    ### sequential 1-d cubic splines
    try({
      fit1d<-fit.1d(dat=ave.dat)
      fit1d$dmomdhdhr<-fit1d$d2mom.dH.dHr
      res_plots(dat1=fit1d,method="Unweighted Sequential 1-d Splines")
      
      source("Spline plots.R")
      
      fit1d$dmomdhdhr<-fit1d$weighted.d2mom.dH.dHr
      fit1d$pred<-fit1d$weighted.pred
      fit1d$cv.pred<-fit1d$weighted.cv.pred
      res_plots(dat1=fit1d,method="Weighted Sequential 1-d Splines")
      
    })
    ### Moving average
    source("fit_moving_ave.R")
    for(window_size in c(3,9)){
      MA=NULL
      try(MA<-fit_moving_ave(dat=ave.dat,window_size=window_size))
      try(res_plots(dat1=MA$surface,dat2=MA$derivative,method=paste("Moving Window of size",window_size)))
    }
    if(File.name=="FORC012514"){  ### Conduct 2d local polynomial fittings for FORC012514
  
      sd_compare<-data.frame(sd.res=NULL,
                             sd.cv.res=NULL,
                             deg=NULL,
                             nn=NULL,
                             exclude=NULL)
      degs<-2:3
      FIT_LIST<-list()
      for(deg in degs){
        if(deg==2) NN<-c(.0005,.00055,.0006,.00065,.0007,.0008,.0009,.001,.0011,.0012,.0015,.002,.0025,.005,.01)
        if(deg==3) NN<-c(.0008,.0009,.001,.0011,.0012,.0015,.002,.0025,.005,.01,.02)
        for(nn in NN){
          exclude_list<-0
          if((nn==.001&deg==3)|(nn%in%c(.0006,.00065)&deg==2)){exclude_list<-0:1}
          for(exclude_0 in exclude_list){
            print(c(nn,deg))
            t.df<-ave.dat;exclude="";
            exclude_rng=NULL
            if(exclude_0){
              t.df<-ave.dat[ave.dat$H<(-10.5)|ave.dat$H>13,]
              exclude<-"Ex -10.5:13"
              exclude_rng<-c(-10.5,13)
            } 
            try({
              mod<-locfit(Moment~lp(H,H_r,deg=deg,nn=nn),data=t.df)
              t.df$pred<-fitted(mod)
              deriv.mod<-locfit(Moment~lp(H,H_r,deg=deg,nn=nn),data=t.df,deriv=c(1,2))
              t.df$dmomdhdhr<-fitted(deriv.mod)
              {
                t.df$cv.pred<-NA
                ### Try using cross-validation by paths to detect over-fitting
                for(k in 1:max(t.df$path)){
                  mod<-locfit(Moment~lp(H,H_r,deg=deg,nn=nn),data=t.df[t.df$path!=k,])
                  t.df$cv.pred[t.df$path==k]<-predict(mod,newdata=t.df[t.df$path==k,])
                }
              }
              FIT_LIST[[length(FIT_LIST)+1]]<-t.df
              names(FIT_LIST)[length(FIT_LIST)]<-paste("deg",deg,",nn",nn,exclude_0)
              #save(FIT_LIST,file=file.path(File.name,paste0(orig.file,"_FIT_LIST.Rdata")))
              sd_compare<-rbind(sd_compare,
                                data.frame(sd.res=sd((t.df$Moment-t.df$pred)/t.df$Moment_SE),
                                           sd.cv.res=sd((t.df$Moment-t.df$cv.pred)/t.df$Moment_SE),
                                           deg=deg,
                                           nn=nn,
                                           exclude=exclude_0))
              pres_plots(t.df,method=paste("Locfit, deg =",deg,", nn =",nn,exclude),exclude_rng = exclude_rng)
            })
          }
        }
      }
      save(sd_compare,file=file.path(File.name,paste0(orig.file,"_sd_compare.Rdata")))
    }
  }
  dev.off()
}

File.name<-"FORC012514"
  setwd("C:/Users/spl/Documents/Magnetic Nanoparticle with Cindi Dennis/")
load("FORC012514/FORC012514_sd_compare.Rdata")
load("FORC012514/FORC012514_FIT_LIST.Rdata")
length(FIT_LIST)
nn<-substr(names(FIT_LIST),11,nchar(names(FIT_LIST))-2)%>%as.numeric
deg<-substr(names(FIT_LIST),5,5)%>%as.numeric

sd_compare<-  data.frame(sd.res=sapply(FIT_LIST,function(t.df) sd((t.df$Moment-t.df$pred)/t.df$Moment_SE)),
                         sd.cv.res=sapply(FIT_LIST,function(t.df) sd((t.df$Moment-t.df$cv.pred)/t.df$Moment_SE)),
                         deg=deg,
                         nn=nn,
                         exclude=substr(names(FIT_LIST),nchar(names(FIT_LIST)),nchar(names(FIT_LIST))))

sd_compare<-sd_compare[sd_compare$exclude==0,]

par(mfrow=c(1,1))
par(mai=c(1,1,.1,.1))
zzz<-rt(2e3,24)
plot(density(zzz),xlab="Standardized Residual",main="",
     ylab="Density",type="l",cex.lab=1.7,cex.axis=1.7,lwd=3)
grid()

pdf(file.path(File.name,paste0(orig.file,"_Best_Locfit_Fitting_results_05232020.pdf")),height=8,width=12)
res_plots(FIT_LIST[[14]],method="Quadratic nn = 0.002")
res_plots(FIT_LIST[[31]],method="Cubic nn = 0.003")
dev.off()

dat1<-FIT_LIST[[14]]

fit_deg2<-read_xlsx("NN analysis.xlsx",sheet=1)
fit_deg2$deg=2

fit_deg3<-read_xlsx("NN analysis.xlsx",sheet=2)
fit_deg3$deg=3

fits<-rbind(fit_deg2,fit_deg3)
names(fits)[1]<-"nn"

fits$sd.res<-fits$sd.cv.res<-NA
for(i in 1:nrow(fits)){
  ind<-which(sd_compare$nn==fits$nn[i]&sd_compare$deg==fits$deg[i])
  if(length(ind)>0){
    fits$sd.res[i]<-sd_compare$sd.res[ind[1]]
    fits$sd.cv.res[i]<-sd_compare$sd.cv.res[ind[1]]
  } else print(c(fits$deg[i],fits$nn[i]))
}

fits<-fits[,-which(names(fits)%in%c("Hc fitting comments","Hb fitting comments"))]
#[!is.na(fits$sd.cv.res)&fits$sd.cv.res<30,]
t.fits<-fits[fits$nn<.03&
               (fits$nn>.0007|fits$deg==2),]%>%
  dplyr::select(nn,deg,sd.cv.res,sd.res,Hc1,Hb1)%>%
  gather(key="key",value="value",sd.cv.res,sd.res,Hc1,Hb1)
t.fits$deg[t.fits$deg==2]<-"Quadratic"
t.fits$deg[t.fits$deg==3]<-"Cubic"
#t.fits$deg[t.fits$key=="Hb1"]<-"Cubic"
t.fits$deg<-factor(t.fits$deg,levels=c("Quadratic","Cubic"))
t.fits$key<-factor(t.fits$key,levels=c("sd.cv.res","sd.res","Hc1","Hb1"))
levels(t.fits$key)<-c("'SD with CV'","SD","H[C]","H[B]")
p<-ggplot(data=t.fits,aes(x=nn,y=value))+
  facet_grid(key~deg,scales="free_y",labeller = label_parsed)+
  scale_x_log10()+
  #scale_y_log10()+
  xlab("Smoothing Factor")+
  ylab("")+
  geom_line()+
  geom_point()+
  theme_bw()
print(p)

sd_compare<-sd_compare%>%arrange(deg,nn)
par(mfrow=c(1,1),mai=c(1,1.2,.2,.2))
plot(NULL,NULL,log="xy",
     xlim=range(sd_compare$nn),
     ylim=range(c(sd_compare$sd.res,sd_compare$sd.cv.res)),
     xlab="Smoothing Factor",
     ylab="Standard Deviation of
Standardized Residuals",cex.lab=1.3,cex.axis=1.3)
grid()
for(deg in 2:3){
  ind<-sd_compare$deg==deg&sd_compare$exclude==0
  lines(sd_compare$nn[ind],sd_compare$sd.cv.res[ind],col=deg-1,lwd=2)
  lines(sd_compare$nn[ind],sd_compare$sd.res[ind],lty=2,col=deg-1,lwd=2)
}
legend("top",legend=c("Quadratic (CV)","Quadratic (No CV)","Cubic (CV)","Cubic (No CV)"),lty=c(1,2,1,2),lwd=2,col=c(1,1,2,2),cex=1.28)






