#!/usr/bin/Rscript --slave
# Date: 8th Feb 2014 #
# Reading command line arguments #

SpliceNet <- function(IP,OP,A)
{
input<-IP
output<-OP
outpval<-paste(output,"Pval.txt",sep="-")
alpha<-as.numeric(A)
library(pracma)
library(MSBVAR)
library(tseries)
library(base)
IT<-read.table(input, header = FALSE, sep = "\t")
IM<-as.matrix(IT)
# Remore Column Exon corr
IM<-IM[,-2] 
# Copy ID's
ID<-IM[,1] 
IDU<-unique(ID) 
# Remove ID's
IM<-IM[,-1]
D<-length(ID)
RM<-matrix(data=NA,nrow=(length(IDU)),ncol=(length(IDU)))
PM<-matrix(data=NA,nrow=length(IDU),ncol=length(IDU))
for (i in 1:length(IDU))
{
	RM[i,i]=1;
	PM[i,i]=0;
}

sp1=1 
ep1=1
ii=0
for (i in sp1:D)
{
	if((ID[i]!=ID[i+1])|| (i == D))
	{
		ii=ii+1
		ep1=i
		X1<-IM[sp1:ep1,]
		#cat(ID[sp1])
		#cat("\n")
		#cat(sp1)
		#cat("\t")
		#cat(ep1)
		#cat("\n")
		#prmatrix(X1, right = TRUE, quote = FALSE)
        X1<-t(X1)
		if(sp1==ep1)
		{
			X1<-t(X1)
		}
		#prmatrix(X1, right = TRUE, quote = FALSE)
		sp2=1
		ep2=1
		#if(ep1 < D)
		{
		jj=0;
		for(j in sp2:D)
		{
			if((ID[j]!=ID[j+1]) || (j == D))
			{
				jj=jj+1;
				ep2=j
				#if( (sp1 != sp2) && (ep1 != ep2))
				{
				X2<-IM[sp2:ep2,]
				#cat(ID[sp2])
				#cat("\n")
				#cat(sp2)
                #cat("\t")
                #cat(ep2)
				#cat("\n")
                #prmatrix(X2, right = TRUE, quote = FALSE)
				X2<-t(X2)
				if(sp2==ep2)
				{
					X2<-t(X2)
				}
				#prmatrix(X2, right = TRUE, quote = FALSE)
				#Code for LDT #
				n1=nrow(X1)
				n2=nrow(X2)
				p1=ncol(X1)
				p2=ncol(X2)
				N1=n1-1
				N2=n2-1
				#prmatrix(cbind(X1,X2),right = TRUE, quote = FALSE)
				X3<-cbind(X1,X2)
				X3<-matrix(X3,nrow=nrow(X3),ncol=ncol(X3))
				#prmatrix(X3, right = TRUE, quote = FALSE)
				mode(X3) <- 'numeric'
				Sigma=cov(X3)
				A11=N1*Sigma[1:p1,1:p1]
				A22=N1*Sigma[(p1+1):(p1+p2),(p1+1):(p1+p2)]
				A12=N1*Sigma[1:p1,(p1+1):(p1+p2)]
				A21=N1*Sigma[(p1+1):(p1+p2),1:p1]
				A11=as.matrix(A11)
				A22=as.matrix(A22)
				A12=matrix(A12,p1,p2)
				A21=matrix(A21,p2,p1)
				M=A21%*%inv(A11)%*%A12%*%inv(A22)
				#r1N=p2/p1
				#r2N=p2/(N1-1-p1)
				r1=p2/p1
				r2=p2/(n1-1-p1)
				hn=sqrt(r1+r2-r1*r2)
				eg=r2/(r1+r2)
				vg=2*(hn^2)*(r1^2)*(r2^2)/((r1+r2)^4)
				TS=0
				TS=(vg^(-0.5))*(sum(diag(M))-p2*eg)
				if ((is.infinite(TS))|(is.nan(TS))|(is.null(TS)))
				{
						TS=-10000
				}
				result=NULL
				result[1]=TS
				pvalue=1-pnorm(TS)
				#cat(TS)
				cvi=1-alpha
				CV=qnorm(cvi)
				result[2]=CV
				#cat(CV)
				if (TS >= CV)
				{
        			result[3]=1
					#cat(ID[sp1]);cat("\t");cat(ID[sp2]);cat("\t");cat("1\t");cat(pvalue);cat("\n");
				}
				else
				{
        			result[3]=0
					#cat(ID[sp1]);cat("\t");cat(ID[sp2]);cat("\t");cat("0\t");cat(pvalue);cat("\n");
				}
				RM[ii,jj]<-result[3]
				PM[ii,jj]<-pvalue
				sp2=j+1
				#jj=jj+1
				} # for sp1 = sp2
			} #J If
		} #J For exta if
		} # Extra if
	sp1=i+1
	#ii=ii+1
	} # I If
} # I for

RM<-rbind(IDU,RM)
PM<-rbind(IDU,PM)
IDU<-c("I-MATRIX",IDU)
RM<-cbind(IDU,RM)
PM<-cbind(IDU,PM)
#for (i in 1:(length(IDU)+1))
#{
#	for (j in 1:(length(IDU)+1))
#	{
#		#cat(RM[i,j])
#		#cat("\t")
#	}
#	cat("\n")
#}

write(RM, file = output, ncolumns = ncol(RM), append = FALSE, sep = "\t")
write(PM, file = outpval, ncolumns = ncol(PM), append = FALSE, sep = "\t")
} 

SpliceNetDN <- function(NN,CN)
{
input1<-NN
output2<-CN
IT1<-read.table(input1, header = FALSE, sep = "\t")
IM1<-as.matrix(IT1)
ID1<-IM1[,1] 
IM1<-IM1[,-1]
IM1<-IM1[-1,]
ID1<-ID1[-1]

IT2<-read.table(input2, header = FALSE, sep = "\t")
IM2<-as.matrix(IT2)
ID2<-IM2[,1] 
IM2<-IM2[,-1]
IM2<-IM2[-1,]
ID2<-ID2[-1]

D<-(length(ID1))
key<-c("","")
temp<-c("")
for (i in 1:D)
{
	for (j in 1:D)
	{
		
		if((is.na(IM1[i,j])) && (is.na(IM2[i,j]))) {}
		else
		{
			if(IM1[i,j] != IM2[i,j])
			{
				key<-sort(c(ID1[i],ID1[j]))
				val<-paste(key[1],"\t",key[2],"\t",IM1[i,j],"\t",IM2[i,j],"\n",sep="")
				temp<-c(temp,val)
				#cat(key[1]);cat("\t");cat(key[2]);cat("\t");cat(IM1[i,j]);cat("\t");cat(IM2[i,j]);cat("\n");
			}
		}
		
	}
}
#cat(length(temp));cat("\n");
temp<-unique(temp)
#cat(length(temp));cat("\n");
for( i in 1:length(temp))
{
	cat(temp[i])
}
}

