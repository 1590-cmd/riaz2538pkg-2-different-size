gr1<-function(v,p1,p2,i){
  #v=9
  #p=4
  #if(i<=0)stop("i must be a positive integer")
  if(v%%2!=0){
    if((v-1-p2)%%p1==0)
    {
      i=(v-1-p2)/p1
      if(i<=0)stop("i must be a positive integer")
      cat("MCBRMD IS POSSIBLE with one set of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=1")
    }else
      if((v-1-2*p2)%%p1==0)
      {
        i=(v-1-2*p2)/p1
        if(i<=0)stop("i must be a positive integer")
        cat("MCBRMD IS POSSIBLE with two sets of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=1")
      }else
        if((v-p2)%%p1==0)
        {
          i=(v-p2)/p1
          if(i<=0)stop("i must be a positive integer")
          cat("MCSBRMD IS POSSIBLE with one set of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=2")
        }else
          if((v-2*p2)%%p1==0)
          {
            i=(v-2*p2)/p1
            if(i<=0)stop("i must be a positive integer")
            cat("MCSBRMD IS POSSIBLE with two sets of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=2")
          }else
            cat("No Design IS POSSIBLE IN two different PERIOD SIZES FOR GIVEN", "V=",v, "P1=",p1,  "AND", "P2=",p2)  
  }   
    if(v%%2==0){
      if((v-2-p2)%%p1==0)
      {
        i=(v-2-p2)/p1
        if(i<=0)stop("i must be a positive integer")
        cat("MCPBRMD-I IS POSSIBLE with one set of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=3")
      }else
        if((v-2-2*p2)%%p1==0)
        {
          i=(v-2-2*p2)/p1
          if(i<=0)stop("i must be a positive integer")
          cat("MCPBRMD-I IS POSSIBLE with two sets of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=3")
        }else
          
            if((v-1-p2)%%p1==0)
            {
              i=(v-1-p2)/p1
              if(i<=0)stop("i must be a positive integer")
              cat("MCSPBRMD-I IS POSSIBLE with one set of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=4")
            }else
              if((v-1-2*p2)%%p1==0)
              {
                i=(v-1-2*p2)/p1
                if(i<=0)stop("i must be a positive integer")
                cat("MCSPBRMD-I IS POSSIBLE with two sets of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=4")              
              }else
                if((v-p2)%%p1==0)
                {
                  i=(v-p2)/p1
                  if(i<=0)stop("i must be a positive integer")
                  cat("MCWBRMD-I IS POSSIBLE with one set of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=5")
                }else
            if((v-2*p2)%%p1==0)
            {
              i=(v-2*p2)/p1
              if(i<=0)stop("i must be a positive integer")
              cat("MCWBRMD-I IS POSSIBLE with two sets of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=5")
            }else
              if((v+1-p2)%%p1==0)
              {
                i=(v+1-p2)/p1
                if(i<=0)stop("i must be a positive integer")
                cat("MCSBGRMD-I IS POSSIBLE with one set of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=6")
              }else
                if((v+1-2*p2)%%p1==0)
                {
                  i=(v+1-2*p2)/p1
                  if(i<=0)stop("i must be a positive integer")
                  cat("MCSBGRMD-I IS POSSIBLE with two sets of p2", "V=",v, "P1=",p1, "P2=",p2, "AND", "i=",i, "AND ALSO RUN D=6")
                }else
                  cat("No Design IS POSSIBLE IN two different PERIOD SIZES FOR GIVEN", "V=",v, "P1=",p1,  "AND", "P2=",p2)
}
}
#################################################################################
# Minimal CBRMD,CSBRMD,CWBRMD-1 and CWBRMD-2_equalsize:Minimal Circular Balanced, 
# Strongly Balanced and Weakly Balanced Repeated Measurements Designs Design 
# for period of two different sizes (p1 and p2)

# Algorithm from paper:

#  Muhammad Riaz, M. H. Tahir, Zawar Hussain,
# Abid khan and Rashid Ahmed*#
# Coded by Riaz et al., 2021-2022 
# Version 2.1.0  (2022-04-20)
#################################################################################



#################################################################################
# Selection of i groups of period size p from adjusted A. The set of remaining 
# (Unselected) elements are saved in the object named as B2. 
#################################################################################
grouping1<-function(A,p,v,i){
  bs<-c()
  z=0;f=1
  A1=A
  while(f<=i){
    
    for(y in 1:5000){
      comp<-sample(1:length(A1),p[1])
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
  list(B1=bs,B2=A1)
}

#################################################################################
# Selection of i group of size p1 from adjusted A and division of required 
# number of groups of size p2 from B2. 
#################################################################################
grouping2<-function(A,p,v,i,sp2){
  bs1<-c()
  j=i+sp2
  z=0;f=1
  A1=A
  while(f<=j){
    s<-grouping1(A1,p[1],v,i)
    A2<-s$B2
    z=i;f=f+i
    if(j-z==1){A2<-c(0,A2)}
    for(y in 1:1000){
      comp<-sample(1:length(A2),p[2])
      com<-A2[comp]
      cs<-sum(com)
      if(cs%%v==0){
        bs1<-rbind(bs1,com)
        A2<-A2[-comp]
        z<-z+1
        f=f+1
        if(j-z==1){A2=c(0,A2)}
      }
      if(z==j) break
    }
    
    
    if(z<j) {bs1<-c();z=0;f=1;A1=A}  
    
  }
  
  
  gs1<-t(apply(s$B1,1,sort))
  gs1<-cbind(gs1,rowSums(gs1),rowSums(gs1)/v,(i+sp2)*v)
  rownames(gs1)<-paste("G",1:i, sep="")
  colnames(gs1)<-c(paste(1:p[1], sep=""),"sum" ,"sum/v","n")
  
  gs2<-t(apply(bs1,1,sort))
  gs2<-cbind(gs2,rowSums(gs2),rowSums(gs2)/v,(i+sp2)*v)
  rownames(gs2)<-paste("G",(nrow(gs1)+1):(nrow(gs1)+sp2), sep="")
  colnames(gs2)<-c(paste(1:p[2], sep=""),"sum" ,"sum/v","n")
  
  
  fs1<-t(apply(s$B1,1,sort))
  fs1<-delmin(fs1)
  rownames(fs1)<-paste("S",1:i, sep="")
  colnames(fs1)<-rep("",(p[1])-1)
  
  
  fs2<-t(apply(bs1,1,sort))
  fs2<-delmin(fs2)
  rownames(fs2)<-paste("S",(i+1):(i+sp2), sep="")
  colnames(fs2)<-rep("",(p[2]-1))
  
  list(B1=list(fs1,fs2),B3=list(gs1,gs2),B4=A2)
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
  return(fs)
}

#################################################################################
# Selection of adjusted A and the set(s) of shifs to obtain Circular Balanced, Strongly Balanced and Weakly Balanced 
# Repeated Measurements Designs for period of two different sizes. 
#################################################################################

# D=1: Minimal CBRMDs 
# D=2: Minimal CSBRMDs 
# D=3: Minimal CPBRMDs-I
# D=4: Minimal CSPBRMDs-I
# D=5: Minimal CWBRMDs-I
# D=6: Minimal CSBGRMDs-I 
#   p: Vector of three different period sizes 
#   i: Number of sets of shifts for p1
# Sp2: Number of sets of shifts for p2

CGSBRMD_2diffsize<-function(v,p,i,D=1,sp2=1){
  
  if(length(p)>2 | length(p)<2){stop("Must be length(p)=2 ")}
  if(any(p<=2)!=0) stop("p=Period sizes: Each period size must be greater than 2")
  if(i<=0) stop("i= Must be a positive integer")
  if(p[1]==p[2]) stop("Must be fullfill this condition: p1 not equal to p2")
  setClass( "stat_test", representation("list"))
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following are required sets of shifts to obtain the minimal CBRMDs,CSBRMDs 
and  CWBRMDs for", "v=" ,object$R[1], ",","p1=",object$R[2],",", "p2=",object$R[3],
"and","n=",(i+sp2)*v,"\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
    cat("\n")
    print(object$S[[2]])
  })
  
  if(sp2==1){
    if(v%%2!=0 & D==1){ 
      v=p[1]*i+p[2]+1
      A<-c(1,2:(v-1))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    if(v%%2!=0 & D==2){ 
      v=p[1]*i+p[2]
      A<-c(0,1:(v-1))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    
    if(v%%2==0 & D==3){
      v=p[1]*i+p[2]+2
      A<-c(1:((v-2)/2),((v+2)/2),((v+4)/2):(v-1))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    
    if(v%%2==0 & D==4){
      v=p[1]*i+p[2]+1
      A<-c(0,1:((v-2)/2),((v+2)/2),((v+4)/2):(v-1))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    
    if(v%%2==0 & D==5){  
      v=p[1]*i+p[2]
      A<-c(1:(v-1),(v/2))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    
    if(v%%2==0 & D==6){  
      v=p[1]*i+p[2]-1
      A<-c(0,1:(v-1),(v/2))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    
  }
  if(sp2==2){
    if(v%%2!=0 & D==1){ 
      v=p[1]*i+2*p[2]+1
      A<-c(1,2:(v-1))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    if(v%%2!=0 & D==2){ 
      v=p[1]*i+2*p[2]
      A<-c(0,1:(v-1))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    
    if(v%%2==0 & D==3){
      v=p[1]*i+2*p[2]+2
      A<-c(1:((v-2)/2),((v+2)/2),((v+4)/2):(v-1))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    
    if(v%%2==0 & D==4){
      v=p[1]*i+2*p[2]+1
      A<-c(0,1:((v-2)/2),((v+2)/2),((v+4)/2):(v-1))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    
    if(v%%2==0 & D==5){  
      v=p[1]*i+2*p[2]
      A<-c(1:(v-1),(v/2))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    
    if(v%%2==0 & D==6){  
      v=p[1]*i+2*p[2]-1
      A<-c(0,1:(v-1),(v/2))
      A1<-c(grouping2(A,p,v,i,sp2))
      A2<-c(v,p);names(A2)<-c("V","p1","p2")
      x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
    }
    
    
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

###############################################################################
# Examples: Using CBRMD,CSBRMD AND CWBRMD_equalsize function to obtain the set(s) of shifts
# for construction of Circular Balanced, Strongly Balance and weakly balanced Repeated Measure 
# Design for period of two different sizes (p1 and p2)
###############################################################################


(H<-gr1(v=13,p1=3,p2=4))

# Examples for all cases D=1-case-1
p=c(5,4);i=3;D=1;sp2=1
(H<-CGSBRMD_2diffsize(p,i,D,sp2))
(H<-CGSBRMD_2diffsize(p=c(3,4),v=13,i=3,D=2,sp2=1))
H$G
(H<-CGSBRMD_2diffsize(p=c(5,3),v=11,i=3,D=8,sp2=2))
(design_CGSBRMD(H))

# Examples for all cases D=2-case-1
p=c(5,4);i=3;D=2;sp2=2
(H<-CGSBRMD_2diffsize(p,i,D,sp2))
(H<-CGSBRMD_2diffsize(p=c(5,3),v=23,i=3,D=1,sp2=2))
H$G
(H<-CGSBRMD_2diffsize(p=c(7,3),v=24,i=3,D=3,sp2=1))
(H<-CGSBRMD_2diffsize(p=c(7,5),v=26,i=3,D=3,sp2=1))
(design_CGSBRMD(H))

# Examples for all cases D=1-case-2
(H<-CGSBRMD_2diffsize(p=c(13,15),v=33,i=2,D=4,sp2=1))
(H<-CGSBRMD_2diffsize(p=c(3,6),v=11,i=5,D=2,sp2=1))
(design_CGSBRMD(H))

# Examples for all cases D=2-case-2
(H<-CGSBRMD_2diffsize(p=c(7,5),v=11,i=3,D=5,sp2=1))
(H<-CGSBRMD_2diffsize(p=c(8,5),v=50,i=5,D=5,sp2=2))
(design_CGSBRMD(H))


# Example with message

H<-CGSBRMD_2diffsize(p=c(10,6),v=11,i=4,D=6,sp2=1)
H<-CGSBRMD_2diffsize(p=c(10,8),v=11,i=4,D=6,sp2=2)
H<-CGSBRMD_2diffsize(p=c(10,8),v=11,i=3,D=7,sp2=1)
H<-CGSBRMD_2diffsize(p=c(10,8),v=11,i=3,D=7,sp2=2)
H<-CGSBRMD_2diffsize(p=c(7,5),v=11,i=3,D=8,sp2=1)




# Examples for all cases D=3-case-1
p=c(5,4);i=3;D=1;sp2=1
(H<-CGSBRMD_2diffsize(p,i,D,sp2))
(H<-CGSBRMD_2diffsize(p=c(6,4),v=24,i=3,D=3,sp2=1))
(H<-CGSBRMD_2diffsize(p=c(4,3),v=13,i=2,D=3,sp2=1))
(H<-CGSBRMD_2diffsize(p=c(7,3),v=32,i=5,D=3,sp2=1))
(design_CGSBRMD(H))

# Examples for all cases D=4-case-1
p=c(5,4);i=3;D=2;sp2=1
(H<-CGSBRMD_2diffsize(p,i,D,sp2))
(H<-CGSBRMD_2diffsize(p=c(7,4),v=32,i=4,D=4,sp2=1))
(H<-CGSBRMD_2diffsize(p=c(7,3),v=24,i=3,D=4,sp2=1))
(H<-CGSBRMD_2diffsize(p=c(7,5),v=26,i=3,D=4,sp2=1))
(design_CGSBRMD(H))

# Examples for all cases D=3-case-2
(H<-CGSBRMD_2diffsize(p=c(8,4),v=34,i=3,D=3,sp2=2))
H$G
(H<-CGSBRMD_2diffsize(p=c(9,6),i=5,D=3,sp2=2))
(design_CGSBRMD(H))

# Examples for all cases D=4-case-2
(H<-CGSBRMD_2diffsize(p=c(7,5),i=3,D=4,sp2=2))
(H<-CGSBRMD_2diffsize(p=c(8,5),v=50,i=5,D=4,sp2=2))
(design_CGSBRMD(H))


# Example with message

H<-CGSBRMD_2diffsize(p=c(10,6),i=4,D=3,sp2=2)
H<-CGSBRMD_2diffsize(p=c(10,8),i=4,D=3,sp2=2)
H<-CGSBRMD_2diffsize(p=c(10,8),i=3,D=3,sp2=2)
H<-CGSBRMD_2diffsize(p=c(10,8),i=3,D=3,sp2=1)
H<-CGSBRMD_2diffsize(p=c(7,5),i=3,D=3,sp2=1)


