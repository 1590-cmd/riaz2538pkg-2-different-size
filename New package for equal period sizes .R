  
#################################################################################
# Algorithm Coded with R to Generate Efficient Circular Balanced Repeated Measurements Designs and Their Useful Classes
# for period of equal size(P)

# Algorithm from paper:

#  Muhammad Riaz, Mahmood ul Hassan,abid khan
#   H. M. Kashif Rasheed and Rashid Ahmed*

#
# Coded by Riaz et al., 2021-2022 
# Version 2.1.0  (2023-04-20)
#################################################################################
Designcheck<-function(v,p,i){
  #v=9
  #p=4
  #if(i<=0)stop("i must be a positive integer")
  if(v%%2!=0){
    if((v-1)%%p==0)
    {
      i=(v-1)/p
      if(i<=0)stop("i must be a positive integer")
      cat("MCBRMD IS POSSIBLE IN EQUAL PERIOD SIZES FOR GIVEN", "V=",v, "i=",i, "AND", "P=",p, "AND ALSO RUN D=1")
    }else
      if(v%%p==0)
      {
        i=v/p
        if(i<=0)stop("i must be a positive integer")
        cat("MCSBRMD IS POSSIBLE IN EQUAL PERIOD SIZES FOR GIVEN", "V=",v,"i=",i,  "AND", "P=",p, "AND ALSO RUN D=2")
      }else
        cat("Neither MCBRMD nor MCSBRMD is possible for", "V=",v,  "AND", "P=",p)
  }
  if(v%%2==0){
    
    if((v-2)%%p==0){
      i=(v-2)/p
      if(i<=0)stop("i must be a positive integer")
      cat("MCPBRMD IS POSSIBLE IN EQUAL PERIOD SIZES FOR GIVEN", "V=",v,"i=",i,  "AND", "P=",p, "AND ALSO RUN D=3")
    }else
      if((v-1)%%p==0){
        i=(v-1)/p
        if(i<=0)stop("i must be a positive integer")
        cat("MCSPRMD IS POSSIBLE IN EQUAL PERIOD SIZES FOR GIVEN", "V=",v,"i=",i,  "AND", "P=",p, "AND ALSO RUN D=4")
      }else    
        if(v%%p==0){
          i=v/p
          if(i<=0)stop("i must be a positive integer")
          cat("MCWBRMD IS POSSIBLE IN EQUAL PERIOD SIZES FOR GIVEN", "V=",v,"i=",i,  "AND", "P=",p, "AND ALSO RUN D=5")
        }else
          if((v+1)%%p==0){
            i=(v+1)/p
            if(i<=0)stop("i must be a positive integer")
            cat("MCSBGRMD IS POSSIBLE IN EQUAL PERIOD SIZES FOR GIVEN", "V=",v,"i=",i,  "AND", "P=",p, "AND ALSO RUN D=6")
          }else
            
            cat("No Design IS POSSIBLE IN EQUAL PERIOD SIZES FOR GIVEN", "V=",v,  "AND", "P=",p)
    
  }
}

################################################################
# Division of adjusted A in i groups to get the set(s) of shifts
################################################################
grouping1<-function(A,p,v,i){
  bs<-c()
  z=0;f=1
  A1=A
  while(f<=i){
    
    for(y in 1:5000){
      comp<-sample(1:length(A1),p)
      com<-A1[comp]
      cs<-sum(com)
      if(cs%%v==0){
        bs<-rbind(bs,com)
        A1<-A1[-comp]
        z<-z+1
        f=f+1
      }
      if(z==i) break
    }
    if(z<i) {bs<-c();z=0;f=1;A1=A}  
    
  }
  
  
  bs1<-t(apply(bs,1,sort))
  bs1<-cbind(bs1,rowSums(bs),rowSums(bs)/v)
  rownames(bs1)<-paste("G",1:i, sep="")
  colnames(bs1)<-c(paste(1:p, sep=""),"sum" ,"sum/v")
  
  bs2<-t(apply(bs,1,sort))
  bs2<-delmin(bs2)
  list(B1=list(bs2),B2=list(bs1),B3=A1)
}


#######################################################################
# Obtaing set(s) of shifts by deleting smallest value of each group
#######################################################################

delmin<-function(z){
  fs<-c()
  n<-nrow(z)
  c<-ncol(z)-1
  for(i in 1:n){
    z1<-z[i,]
    z2<-z1[z1!=min(z1)]
    fs<-rbind(fs,z2)
  }
  rownames(fs)<-paste("S",1:n, sep="")
  colnames(fs)<-rep("",c)
  return(fs)
}


#################################################################################
# Selection of adjusted A and the set(s) of shifs to obtain Circular Balanced, 
#Strongly Balanced, partialy balanced, strongly partialy balanced, Weakly Balanced,
#strongly generalized balanced Repeated Measurements Designs for period of equal size. 
#################################################################################

# D=1: Minimal CBRMDs 
# D=2: Minimal CSBRMDs 
# D=3: Minimal CPBRMDs
# D=4: Minimal CSPBRMDs
# D=5: Minimal CWBRMDs
# D=6: Minimal CSBGRMDs 


#   P: Period size
#   i: Number of set of shifts for P


CGSBRMD_equalsize<-function(v,p,i,D=1){
  
  if(p<=2) stop("P= Period size: Period size must be greater than 2")
  if(i<=0) stop("i= Must be a positive integer")
  
  
  setClass( "stat_test", representation("list"))
  
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following are required sets of shifts to obtain the 
minimal CBRMDs,CSBRMDs and  CWBRMDs for", "v=" ,object[[3]][1], "and","P=",object[[3]][2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
  })
  
  if(v%%2!=0 & D==1){
    v=p*i+1
    A<-c(1,2:(v-1))
    A1<-grouping1(A,p,v,i)
    A2<-c(v,p);names(A2)<-c("V","P")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
  
  if(v%%2!=0 & D==2){
    v=p*i
    A<-c(0,1:(v-1))
    A1<-grouping1(A,p,v,i)
    A2<-c(v,p);names(A2)<-c("V","p")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
  
  if(v%%2==0 & D==3){
    v=p*i+2
    A<-c(1:((v-2)/2),((v+2)/2),((v+4)/2):(v-1))
    A1<-grouping1(A,p,v,i)
    A2<-c(v,p);names(A2)<-c("V","P")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
  
  if(v%%2==0 & D==4){
    v=p*i+1
    A<-c(0,1:((v-2)/2),((v+2)/2),((v+4)/2):(v-1))
    A1<-grouping1(A,p,v,i)
    A2<-c(v,p);names(A2)<-c("V","P")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
  
  
  if(v%%2==0 & D==5){  
    v=p*i
    A<-c(1:(v-1),(v/2))
    A1<-grouping1(A,p,v,i)
    A2<-c(v,p);names(A2)<-c("V","P")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
  
  if(v%%2==0 & D==6){  
    v=p*i-1
    A<-c(0,1:(v-1),(v/2))
    A1<-grouping1(A,p,v,i)
    A2<-c(v,p);names(A2)<-c("V","P")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
  
  
  
  new("stat_test", x)
  
}

##################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from CGSBRMD_equalsize
# The output is called using the design_CGSBRMD to generate design

design_CGSBRMD<-function(H){
  
  setClass( "CGSBRMD_design", representation("list"))
  setMethod("show", "CGSBRMD_design", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following is minimal CGSBRMD for", "v=" ,object$R[1], "and","p=",object$R[2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    for(i in 1:length(ss)){
      W<-ss[[i]]
      nr<-dim(W)[1]
      for(j in 1:nr){
        print(object$Design[[i]][[j]])
        cat("\n\n")
      }}
  })  
  
  v<-H$R[1]
  p<-H$R[2]
  ss<-H$S  
  treat<-(1:v)-1
  fn<-(1:v)
  G<-list()
  
  
  for(j in 1:length(ss)){ 
    W<-ss[[j]]
    nr<-dim(W)[1]
    nc<-dim(W)[2]
    D<-list()
    
    for(i in 1:nr){
      dd<-c()
      d1<-matrix(treat,(nc+1),v,byrow = T)
      ss1<-cumsum(c(0,W[i,]))
      dd2<-d1+ss1
      dd<-rbind(dd,dd2)
      rr<-dd[which(dd>=v)]%%v
      dd[which(dd>=v)]<-rr
      colnames(dd)<-paste("B",fn, sep="")
      rownames(dd)<-rep("",(nc+1))
      fn<-fn+v
      D[[i]]<-dd
    }
    G[[j]]<-D
    
  }
  
  x<-list(Design=G,R=H$R)
  new("CGSBRMD_design", x)
}



##################################################################################
# Examples: Using CBRMD,CSBRMD AND CWBRMD_equalsize function to obtain the set(s) of shifts
# for construction of Circular Balanced, Strongly Balance and weakly balanced Repeated Measure 
# Design for equal period sizes (p)
##################################################################################


# example#
(H<-Designcheck(v=24,p=5))

# example#1
(H<-CGSBRMD_equalsize(v=24,p=5,i=5,D=6))
H$G
(D<-design_CGSBRMD(H))

# example #2
(H<-CGSBRMD_equalsize(v=21,p=5,i=4,D=1))
H$G
design_CGSBRMD(H)

# example #3
(H<-CGSBRMD_equalsize(v=22,p=7,i=2,D=3))
(H<-CGSBRMD_equalsize(v=23,p=7,i=2,D=4))
(H<-CGSBRMD_equalsize(v=18,p=6,i=2,D=5))
(H<-CGSBRMD_equalsize(v=14,p=9,i=3,D=6))



# example #4
(H<-CGSBRMD_equalsize(v=22,p=7,i=3,D=7))
design_CGSBRMD(H)


# example#1
(H<-CGSBRMD_equalsize(v=23,p=7,i=3,D=8))
(D<-design_CGSBRMD(H))

# example #2
(H<-CGSBRMD_equalsize(v=10,p=5,i=2,D=4))
design_CGSBRMD(H)

# example #3
(H<-CGSBRMD_equalsize(v=9,p=3,i=3,D=4))
(H<-CGSBRMD_equalsize(v=6,p=3,i=2,D=4))
H$G
(H<-CGSBRMD_equalsize(p=4,i=2,D=3))
(H<-CGSBRMD_equalsize(p=3,i=2,D=3))
H$G


# example #4
(H<-CGSBRMD_equalsize(p=8,i=7,D=2))
design_CGSBRMD(H)













