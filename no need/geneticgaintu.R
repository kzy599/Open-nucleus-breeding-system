rm(list = ls())
gc()
library(data.table)
library(ggplot2)
library(ggsci)
library(stringr)
library(showtext)
setwd("/home/GTDisk1/kangziyi/onbsresult/")

dt10<- fread("dtallfinal.csv",sep = ",")


dt10[program%flike%"S",ebv_calculation:="ssGBLUP"]
dt10[!program%flike%"S",ebv_calculation:="PBLUP"]


#======================
dt10[program%flike%"1",gradient:="ONBS10%"]
dt10[program%flike%"3",gradient:="ONBS30%"]
dt10[program%flike%"5",gradient:="ONBS50%"]
dt10[program%flike%"C",gradient:="CNBS"]
dt10[program%flike%"0",gradient:="non-overlap"]
#=====================

pname<- unique(dt10$program)
for (i in 1:length(pname)) {
  if(pname[i]%in%c("OP100","OS100")){
    ppname= pname[i]
  }else if(pname[i]%flike%"1"){
   ppname<- paste(substr(pname[i],1,2),"_10",sep = "")
 }else if(pname[i]%flike%"3"){
   ppname<- paste(substr(pname[i],1,2),"_30",sep = "")
 }else if(pname[i]%flike%"5"){
   ppname<- paste(substr(pname[i],1,2),"_50",sep = "")
 }  else{
   ppname<- pname[i]
 }
  dt10[program==pname[i],program:=ppname]
}
gc()
pname<- unique(dt10$program)
colnames(dt10)[2] = "gg"
#第一步需要自己执行i=1，之后2到length可以for循环 因为内存不够
for(nG in 0:20){
  for (i in 1:length(pname)) {
    
      dt_temporary<- as.data.table(cbind(nG,
                               unique(dt10[gg==nG&program==pname[i],program]),
                               unique(dt10[gg==nG&program==pname[i],ebv_calculation]),
                               unique(dt10[gg==nG&program==pname[i],gradient]),
                               mean(dt10[gg==nG&program==pname[i],gv]),
                               sd(dt10[gg==nG&program==pname[i],gv]),
                               mean(dt10[gg==nG&program==pname[i],inb]),
                               sd(dt10[gg==nG&program==pname[i],inb]),
                               mean(dt10[gg==nG&program==pname[i],inb_plink]),
                               sd(dt10[gg==nG&program==pname[i],inb_plink]),
                               mean(dt10[gg==nG&program==pname[i],inbprate]),
                               sd(dt10[gg==nG&program==pname[i],inbprate]),
                               mean(dt10[gg==nG&program==pname[i],inbprates]),
                               sd(dt10[gg==nG&program==pname[i],inbprates]),
                               mean(dt10[gg==nG&program==pname[i],inbrate]),
                               sd(dt10[gg==nG&program==pname[i],inbrate]),
                               mean(dt10[gg==nG&program==pname[i],Ne]),
                               sd(dt10[gg==nG&program==pname[i],Ne]),
                               mean(dt10[gg==nG&program==pname[i],ld]),
                               sd(dt10[gg==nG&program==pname[i],ld]),
                               mean(dt10[gg==nG&program==pname[i],Va]),
                               sd(dt10[gg==nG&program==pname[i],Va]),
                               mean(dt10[gg==nG&program==pname[i],Vd]),
                               sd(dt10[gg==nG&program==pname[i],Vd]),
                               mean(dt10[gg==nG&program==pname[i],Vg]),
                               sd(dt10[gg==nG&program==pname[i],Vg]),
                               mean(dt10[gg==nG&program==pname[i],pgv_MP]),
                               sd(dt10[gg==nG&program==pname[i],pgv_MP]),
                               mean(dt10[gg==nG&program==pname[i],pgv_NP]),
                               sd(dt10[gg==nG&program==pname[i],pgv_NP]),
                               mean(dt10[gg==nG&program==pname[i],prel]),
                               sd(dt10[gg==nG&program==pname[i],prel]),
                               mean(dt10[gg==nG&program==pname[i],nCoancestor]),
                               sd(dt10[gg==nG&program==pname[i],nCoancestor]),
                               mean(dt10[gg==nG&program==pname[i],identicalp]),
                               sd(dt10[gg==nG&program==pname[i],identicalp]),
                               mean(dt10[gg==nG&program==pname[i],rgv10]),
                               sd(dt10[gg==nG&program==pname[i],rgv10]),
                               mean(dt10[gg==nG&program==pname[i],rgv20]),
                               sd(dt10[gg==nG&program==pname[i],rgv20]),
                               mean(dt10[gg==nG&program==pname[i],pgvrat]),
                               sd(dt10[gg==nG&program==pname[i],pgvrat])))
      colnames(dt_temporary) <- c("nGeneration","program","ebv_calcualtion","gradients","gv","gvse","inb","inbse","inbp","inbpse","prate","pratese","prates","pratesse","irate","iratese","Ne","Nese","ld","ldse",
                                  "Va","Vase","Vd","Vdse","Vg","Vgse","MPgv","MPgvse","NPgv","NPgvse","prel","prelse","nCo",'nCose',"identicalp","identicalpse","rgv10","rgv10se","rgv20","rgv20se","pgvr","pgvrse")
      if(nG==0&i==1){
        dt<- dt_temporary
      }else{
        dt<- rbind(dt,dt_temporary)
      }
    
  }
}

dt$gv <- as.numeric(dt$gv)
dt$gvse <- as.numeric(dt$gvse)
dt$Ne <- as.numeric(dt$Ne)
dt$Nese <- as.numeric(dt$Nese)
dt$ld<- as.numeric(dt$ld)
dt$ldse <- as.numeric(dt$ldse)
dt$nGeneration = as.numeric(dt$nGeneration)
dt$inb= as.numeric(dt$inb)
dt$inbse= as.numeric(dt$inbse)
dt$inbpse= as.numeric(dt$inbpse)
dt$inbp= as.numeric(dt$inbp)
dt$prate = as.numeric(dt$prate)
dt$pratese = as.numeric(dt$pratese)
dt$prates = as.numeric(dt$prates)
dt$pratesse = as.numeric(dt$pratesse)
dt$irate = as.numeric(dt$irate)
dt$iratese = as.numeric(dt$iratese)
dt$Va = as.numeric(dt$Va)
dt$Vase = as.numeric(dt$Vase)
dt$Vd = as.numeric(dt$Vd)
dt$Vdse = as.numeric(dt$Vdse)
dt$Vg = as.numeric(dt$Vg)
dt$Vgse = as.numeric(dt$Vgse)
dt$prel = as.numeric(dt$prel)
dt$prelse = as.numeric(dt$prelse)
dt$nCo = as.numeric(dt$nCo)
dt$nCose = as.numeric(dt$nCose)
dt$identicalp = as.numeric(dt$identicalp)
dt$identicalpse = as.numeric(dt$identicalpse)
dt$MPgv = as.numeric(dt$MPgv)
dt$MPgvse = as.numeric(dt$MPgvse)
dt$NPgv = as.numeric(dt$NPgv)
dt$Npgvse = as.numeric(dt$NPgvse)
dt$pgvr = as.numeric(dt$pgvr)
dt$pgvrse = as.numeric(dt$pgvrse)
dt$rgv10 = as.numeric(dt$rgv10)
dt$rgv10se = as.numeric(dt$rgv10se)
dt$rgv20 = as.numeric(dt$rgv20)
dt$rgv20se = as.numeric(dt$rgv20se)
dt[gradients=="ONBS10%",gradients:="ONBS (10%)"]
dt[gradients=="ONBS30%",gradients:="ONBS (30%)"]
dt[gradients=="ONBS50%",gradients:="ONBS (50%)"]
###
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.text = element_text(face = "bold",size = 12),#坐标轴刻度标签加粗
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          #legend.position=c(0.5, 0.95),#图例在绘图区域的位置
          #legend.position="none",
          legend.position="right",
          #legend.direction = "horizontal",
          legend.direction = "vertical",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}
###

#===============genetic gain==============
for (p in c(10,20)) {
  dt_plot= dt[nGeneration==p,]
  dt_plot= dt_plot[!gradients=="non-overlap",]
  if(p == 10){
    dt_plot[,nsbg:="Short-term selection"]
    dt_plot1= dt_plot
  }else{
    dt_plot[,nsbg:="Long-term selection"]
    dt_plot2= dt_plot
  }
}
dt_plot = rbind(dt_plot1,dt_plot2)
dt_plot$nsbg= factor(dt_plot$nsbg,levels = c("Short-term selection","Long-term selection"))
dt_plot[gradients=="CNBS",stype:="CNBS"]
dt_plot[gradients%flike%"ONBS",stype:="ONBS"]
#dt_plot[ebv_calcualtion=="ssGBLUP",ebv_calcualtion:="genomic"]
#dt_plot[ebv_calcualtion=="PBLUP",ebv_calcualtion:="traditional"]
#dt_plot$ebv_calcualtion= factor(dt_plot$ebv_calcualtion,levels = c("traditional","genomic"))
P<- ggplot(data = dt_plot,aes(x = gradients, y = gv, fill = ebv_calcualtion))+
  geom_bar(stat = "identity",position = "dodge",)+
  ylab("Genetic gain")+
  geom_errorbar(aes(ymax = gv+gvse, ymin = gv-gvse),
                position = position_dodge(0.9), width = 0.15)+##不要下标的话把ymin=mean-se去掉
  #ylim(0,35)+
  #ylim(0,70)+
  #scale_fill_brewer(palette = "Set1")+
theme_zg()+facet_wrap(~nsbg)
P
ggsave("Aall.png", P , width = 15, height = 5, dpi = 300)

dt_plot[gradients=="CNBS",gradients:="0%"]
dt_plot[gradients=="ONBS (10%)",gradients:="10%"]
dt_plot[gradients=="ONBS (30%)",gradients:="30%"]
dt_plot[gradients=="ONBS (50%)",gradients:="50%"]
P<- ggplot(data = dt_plot[nsbg == "Short-term selection"],aes(x = stype, y = gv, fill = gradients))+
  geom_bar(stat = "identity",position = "dodge",)+
  ylab("Genetic gain")+
  geom_errorbar(aes(ymax = gv+gvse, ymin = gv-gvse),
                position = position_dodge(0.9), width = 0.15)+##不要下标的话把ymin=mean-se去掉
  #ylim(0,35)+
  #ylim(0,70)+
  #scale_fill_brewer(palette = "Set1")+
  theme_zg()+facet_wrap(~ebv_calcualtion)
P
ggsave("Aall_ppt.png", P , width = 15, height = 5, dpi = 300)

#P <- ggarrange(P1 +ylab("Genetic gain")+xlab("10 generations")+ylim(0,70) , P2 + rremove("ylab")+xlab("20 generations")+ylim(0,70) , # remove axis labels from plots
 #                   #labels = c("a","b"),
  #                  ncol = 1, nrow =2,
   #                 widths = c(1,1),
    #                common.legend = TRUE, legend = "right",
     #               align = "v", 
      #              font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
#ggsave("A20normal.tiff", P , width = 10, height = 5, dpi = 300)

#dt[gradients=="ONBS10%",gradients:="ONBS (10%)"]
#dt[gradients=="ONBS30%",gradients:="ONBS (30%)"]
#dt[gradients=="ONBS50%",gradients:="ONBS (50%)"]

##===============continues plot=========================

dt_plot= dt[!gradients=="non-overlap",]
#dt_plot= dt[ebv_calcualtion=="ssGBLUP",]
#dt_plot= dt[ebv_calcualtion=="PBLUP",]
dt_plot[program%flike%"C",breeding_sys:="CNBS"]
dt_plot[!program%flike%"C",breeding_sys:="ONBS"]
#dt_plot= dt_plot[!gradients=="non-overlap",]
dt_plot[gradients=="CNBS",gradients:=NA]
dt_plot[nGeneration==0,inb:=0]
dt_plot[nGeneration==0,inbp:=0]
#dt_plot[nGeneration==0,ld:=0]
dt_plot[gradients=="ONBS (10%)",gradients:="10"]
dt_plot[gradients=="ONBS (30%)",gradients:="30"]
dt_plot[gradients=="ONBS (50%)",gradients:="50"]


dt_plot[program%flike%"10",program:="ONBS (10%)"]
dt_plot[program%flike%"30",program:="ONBS (30%)"]
dt_plot[program%flike%"50",program:="ONBS (50%)"]
dt_plot[program%flike%"C",program:="CNBS"]


#dt_plot[program%flike%"OP_",program:=paste("OP",substr(dt_plot[program%flike%"OP_",program],4,5),sep = "")]
#dt_plot[program%flike%"OS_",program:=paste("OS",substr(dt_plot[program%flike%"OS_",program],4,5),sep = "")]

#======if the plot is about the relationship between nucleus and multipliers, run the additional code below==========
dt_plot= dt_plot[program%flike%"ONBS",]
dt_plot = dt_plot[nGeneration%in%c(2:20),]
#=======================================================================================


dt_plot[program=="ONBS (10%)"&ebv_calcualtion=="ssGBLUP",.(program,nGeneration,pgvr)]


#dt_plot[ebv_calcualtion=="ssGBLUP",ebv_calcualtion:="genomic"]
#dt_plot[ebv_calcualtion=="PBLUP",ebv_calcualtion:="traditional"]
#dt_plot$ebv_calcualtion= factor(dt_plot$ebv_calcualtion,levels = c("traditional","genomic"))

P1<- ggplot(data = dt_plot,aes(x=nGeneration,y=inbp,group=program,color = program))+
  geom_point()+
  geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  xlab("Generation")+
  ylab("Genomic inbreedingr")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+
  #geom_ribbon(aes(ymin = pgvr -pgvrse, ymax = pgvr+pgvrse), alpha = 0.4)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
  geom_errorbar(aes(ymin=inbp-inbpse,
                    ymax=inbp+inbpse),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(color="Breeding scheme")+scale_color_aaas()+facet_wrap(~ ebv_calcualtion)
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(linetype="Breeding system",shape="Proportion of\nthe introduced\nMP individuals (%)")
P1
P1+ylim(0,150)
#P2
ggsave("inbreeding.pdf", P1, width = 10, height = 5, dpi = 300)
ggsave("nCo.tiff", P1, width = 10, height = 5, dpi = 300)
ggsave("nCo.png", P1+ylim(0,150) , width = 10, height = 5, dpi = 300)
ggsave("Ne.png", P1, width = 10, height = 5, dpi = 300)

ggsave("图片.png", P1, width = 10, height = 5, dpi = 300)



#P2<- ggplot(data = dt_plot,aes(x=nGeneration,y=inbp,group=program,shape = gradients,linetype = breeding_sys))+
 # geom_point()+
#  geom_line()+
#  scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
#  xlab("Generation")+
#  ylab("Genomic inbreeding")+
#  theme_bw()+
#  theme(panel.grid.major = element_line(colour = NA),
 #       panel.background = element_rect(fill = "transparent",colour = NA),
  #      plot.background = element_rect(fill = "transparent",colour = NA),
   #     panel.grid.minor = element_blank(),
    #    text = element_text(family = "STXihei"),
     #   legend.position = "right",
      #  legend.background = element_rect(colour = "black"))+
  #geom_ribbon(aes(ymin = inbp - inbpse, ymax = inbp+inbpse), alpha = 0.3)+
  #scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(fill="Breeding scheme")
#geom_errorbar(aes(ymin=inbp-inbpse,
 #                ymax=inbp+inbpse),
  #          width=0.05)+
#scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(linetype="Breeding system",shape="Proportion of\nthe introduced\nMP individuals (%)")
#P1
#P2

figure <- ggarrange(P1 +ylab("inbreeding")+ylim(0,0.3) , P2 + rremove("ylab")+ylim(0,0.3), P3 + rremove("ylab")+ylim(0,0.3) , # remove axis labels from plots
                labels = c("a","b","c"),
                   ncol = 3, nrow =1,
                    widths = c(1,1),
                    common.legend = TRUE, legend = "right",
                    align = "v", 
                    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
figure
ggsave("anewinbx.png", figure , width = 10, height = 5, dpi = 300)



#============================inb for discussion===================================
#p = ssGBLUP or PBLUP
tidydt = function(p){
  dt_plot= dt[ebv_calcualtion==p,]
  dt_plot[program%flike%"C",breeding_sys:="CNBS"]
  dt_plot[!program%flike%"C",breeding_sys:="ONBS"]
  dt_plot[nGeneration==0,inb:=0]
  dt_plot[nGeneration==0,inbp:=0]
  dt_plot[gradients=="ONBS (50%)",gradients:="overlap"]
  
  if(p == "ssGBLUP"){
    dt_plot = dt_plot[program%in%c("OS_50","OS0","OS100"),]
    dt_plot[program=="OS_50",num:="8"]
    dt_plot[program=="OS0",num:="8"]
    dt_plot[program=="OS100",num:="16"]
  }else{
    dt_plot = dt_plot[program%in%c("OP_50","OP0","OP100"),]
    dt_plot[program=="OP_50",num:="8"]
    dt_plot[program=="OP0",num:="8"]
    dt_plot[program=="OP100",num:="16"]
  }
  return(dt_plot)
}

dt_plot = tidydt(p = "ssGBLUP")
dt_plot = tidydt(p = "PBLUP")
p2<- ggplot(data = dt_plot,aes(x=nGeneration,y=inbp,group=program,shape=num,linetype=gradients))+
  geom_point()+
  geom_line()+
  scale_shape_manual(values = c(0,3,15),na.translate=FALSE)+
  xlab("Generation")+
  ylab("Genomic inbreeding")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+
  scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))+labs(shape="The number of NP families \n used to construct MP",linetype="The overlap between \n NP and MP parents")
p2
ggsave("discussinb.tiff", p2 , width = 10, height = 5, dpi = 300)

(dt_plot[nGeneration==20&program%in%c("OP_50","OP0","OP100"),inbp]-dt_plot[nGeneration==20&program=="OP_50",inbp])/dt_plot[nGeneration==20&program=="OP_50",inbp]

(dt_plot[nGeneration==20&program%in%c("OS_50","OS0","OS100"),inbp]-dt_plot[nGeneration==20&program=="OS_50",inbp])/dt_plot[nGeneration==20&program=="OS_50",inbp]


figure1 <- ggarrange(p1 +xlab("ssGBLUP")+ylim(0,0.3) , p2 + rremove("ylab")+xlab("PBLUP")+ylim(0,0.3) , # remove axis labels from plots
                    #labels = c("a","b"),
                    ncol = 2, nrow =1,
                    widths = c(1,1),
                    common.legend = TRUE, legend = "right",
                    align = "v", 
                    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

#========================================================================================================================================================================

dt10<- fread("dtallfinal.csv",sep = ",")

dtop = dt10[program%flike%"OP",]
dtos = dt10[program%flike%"OS",]
dtonbs =rbind(dtop,dtos) 
dtcp = dt10[program%flike%"CP",]
dtcs = dt10[program%flike%"CS",]
dtcnbs = rbind(dtcp,dtcs)
dtpblup = rbind(dtop,dtcp)
dtssgblup = rbind(dtos,dtcs)
dtmodel= dtop[gg==20,]
dtmodel= dtos[gg==20,]
dtmodel= dtonbs[gg==20,]
dtmodel= dtcs[gg==20,]
dtmodel= dtcp[gg==20,]
dtmodel= dtcnbs[gg==20,]
dtmodel = dtpblup[gg==20,]
dtmodel = dtssgblup[gg==20,]

#忽略了世代的影响，只有一个变量program
mymodel<- aov(gv~program,data= dtos[gg==20&program%in%c("OS_50","OS_10","OS_30"),])

mymodel<- aov(gv~program,data= dtop[gg==20&program%in%c("OP_50","OP_10","OP_30"),])


summary(mymodel)
t.test(x = dtos[gg==20&program=="OS_50",inbprate],y=dtos[gg==20&program=="OS100",inbprate]) 
t.test(x = dtop[gg==20&program=="OP_50",inbprate],y=dtop[gg==20&program=="OP100",inbprate])  

t.test(x = dtop[gg==20&program=="OP_50",conrate],y=dtcp[gg==20&program=="CP",conrate])  

t.test(x = dtos[gg==20&program=="OS_10",conrate],y=dtcs[gg==20&program=="CS",conrate])  

pva = c()
for(g in 1:20){
  chk= t.test(x = dtop[gg==g&program=="OP_10",gv],y=dtop[gg==g&program=="OP_50",gv])  
  pva = c(pva,chk$p.value)
}

t.test(x = dtop[gg==20&program=="OP_10",gv],y=dtop[gg==20&program=="OP_50",gv])  
t.test(x = dtop[program=="OP_10",gv],y=dtop[program=="OP_30",gv])


t.test(x = dtos[gg==20&program=="OS_10",gv],y=dtos[gg==20&program=="OS_30",gv])
t.test(x = dtos[program=="OS_10",gv],y=dtos[program=="OS_50",gv])


t.test(x = dtos[gg==20&program=="OS_50",Vg],y=dtos[gg==20&program=="OS_10",Vg])
t.test(x = dtop[gg==20&program=="OP_50",ld],y=dtop[gg==20&program=="OP_10",ld])
t.test(x = dtop[gg==20&program=="OP100",gv],y=dtop[gg==20&program=="OP_50",gv])
t.test(x = dtos[gg==20&program=="OS_10",inb_plink],y=dtos[gg==20&program=="OS_50",inb_plink])


dtmodel= dt20[gg==20&program%flike%"OP",]

#genetic gain
tempvalue = (dt[nGeneration==10&program%in%c("OS_10","OS_30","OS_50"),gv] - dt[nGeneration==10&program=="CS",gv])/dt[nGeneration==10&program=="CS",gv]
round(tempvalue*100,2)
tempvalue = (dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),gv] - dt[nGeneration==20&program=="CS",gv])/dt[nGeneration==20&program=="CS",gv]

(dt[nGeneration==10&program%in%c("OS_10","OS_30","OS_50"),gv] - dt[nGeneration==10&program=="CP",gv])/dt[nGeneration==10&program=="CP",gv]
(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),gv] - dt[nGeneration==20&program=="CP",gv])/dt[nGeneration==20&program=="CP",gv]


tempvalue = (dt[nGeneration==10&program%in%c("OP_10","OP_30","OP_50"),gv] - dt[nGeneration==10&program=="CP",gv])/dt[nGeneration==10&program=="CP",gv]
round(tempvalue*100,2)
tempvalue = (dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),gv] - dt[nGeneration==20&program=="CP",gv])/dt[nGeneration==20&program=="CP",gv]


(dt[nGeneration==10&program%in%c("OS_10","OS_30","OS_50"),gv] - dt[nGeneration==10&program%in%c("OP_10","OP_30","OP_50"),gv])/dt[nGeneration==10&program%in%c("OP_10","OP_30","OP_50"),gv]
(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),gv] - dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),gv])/dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),gv]

(dt[nGeneration==20&program=="CS",gv] - dt[nGeneration==20&program=="CP",gv])/dt[nGeneration==20&program=="CP",gv]
(dt[nGeneration==10&program=="CS",gv] - dt[nGeneration==10&program=="CP",gv])/dt[nGeneration==10&program=="CP",gv]
#rae of inbreeding have already calculated during the process of abstract data
#pname= unique(dt$program)
#NG=10
#NG=20
#for(p in 1:length(pname)){
#  dtneed= dt[program==pname[p],]
#  inbrate = c(numeric(length = NG))
##  inbprate = c(numeric(length = NG))
#  for(g in 2:(NG+1)){
#    inbrate[g] = (dtneed$inb[g] - dtneed$inb[(g-1)])/(1-dtneed$inb[(g-1)])
#    inbprate[g] = (dtneed$inbp[g] - dtneed$inbp[(g-1)])/(1-dtneed$inbp[(g-1)])
#  }
#  if(p==1){
  #  prate= mean(inbprate)
  #  rate = mean(inbrate)
  #}else{
  #  prate= c(prate,mean(inbprate))
  #  rate = c(rate,mean(inbrate))
  #}
#  if(p==length(pname)){
#    finalrate20= data.table(pname,prate,rate)
#  }
#}
(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),inb] - dt[nGeneration==20&program=="CS",inb])/dt[nGeneration==20&program=="CS",inb]
tempvalue =(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),inbp] - dt[nGeneration==20&program=="CS",inbp])/dt[nGeneration==20&program=="CS",inbp]

round(tempvalue*100,2)

(dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),inb] - dt[nGeneration==20&program=="CP",inb])/dt[nGeneration==20&program=="CP",inb]
tempvalue =(dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),inbp] - dt[nGeneration==20&program=="CP",inbp])/dt[nGeneration==20&program=="CP",inbp]

(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),inb] - dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),inb])/dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),inb]
(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),inbp] - dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),inbp])/dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),inbp]
(dt[nGeneration==20&program=="CS",inbp] -  dt[nGeneration==20&program=="CP",inbp])/ dt[nGeneration==20&program=="CP",inbp]

dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50","OP_10","OP_30","OP_50","CS","CP"),.(program,prate,pratese)]

#ld
tempvalue =(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),ld] - dt[nGeneration==20&program=="CS",ld])/dt[nGeneration==20&program=="CS",ld]
tempvalue =(dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),ld] - dt[nGeneration==20&program=="CP",ld])/dt[nGeneration==20&program=="CP",ld]
tempvalue =(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),ld] - dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),ld])/dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),ld]
tempvalue =(dt[nGeneration==20&program=="CS",ld] -  dt[nGeneration==20&program=="CP",ld])/ dt[nGeneration==20&program=="CP",ld]
round(tempvalue*100,2)
#Ne
(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),Ne] - dt[nGeneration==20&program=="CS",Ne])/dt[nGeneration==20&program=="CS",Ne]
(c(67,60,57)-86)/86

(dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),Ne] - dt[nGeneration==20&program=="CP",Ne])/dt[nGeneration==20&program=="CP",Ne]
((c(46,43,42)-79)/79 * 100 ) %>% round(2)

dt[nGeneration==20&program=="CP",.(Ne,Nese)]
dt[nGeneration==20&program=="CS",.(Ne,Nese)]


(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),Ne] - dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),Ne])/dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),Ne]
(c(67,60,57) - c(46,43,42)) /c(46,43,42)

(dt[nGeneration==20&program=="CS",Ne] -  dt[nGeneration==20&program=="CP",Ne])/ dt[nGeneration==20&program=="CP",Ne]
(86-79)/79

#genetic variance
tempvalue =(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),Va] - dt[nGeneration==20&program=="CS",Va])/dt[nGeneration==20&program=="CS",Va]
round(tempvalue*100,2)
tempvalue =(dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),Va] - dt[nGeneration==20&program=="CP",Va])/dt[nGeneration==20&program=="CP",Va]


(dt[nGeneration==20&program%in%c("OS_10","OS_30","OS_50"),Va] - dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),Va])/dt[nGeneration==20&program%in%c("OP_10","OP_30","OP_50"),Va]


(dt[nGeneration==20&program=="CS",Va] -  dt[nGeneration==20&program=="CP",Va])/ dt[nGeneration==20&program=="CP",Va]


#The gv of introduced MP individuals compared with the NP individuals of the NP parents 
dt_onbs= dt[program%in%c("OS_10","OS_30","OS_50","OP_10","OP_30","OP_50")&nGeneration%in%c(2:20),]
ckck= (dt_onbs$MPgv - dt_onbs$NPgv)/dt_onbs$NPgv
ckck = dt_onbs[,pgvr]
max(ckck)
min(ckck)
mean(ckck)

dtonbs[gg%in%c(2:20),pgvrat] %>% max()
dtonbs[gg%in%c(2:20),pgvrat] %>% min()

dt_onbs10 = dt[program%in%"OP_10"&nGeneration%in%c(2:20),]
dt_onbs30 = dt[program%in%"OP_50"&nGeneration%in%c(2:20),]

mean((dt_onbs10$MPgv - dt_onbs30$MPgv)/dt_onbs30$MPgv)

dt_onbs= dt[program%in%"OP_30"&nGeneration%in%c(2:20),]
ckck= (dt_onbs$MPgv - dt_onbs$NPgv)/dt_onbs$NPgv
ckck = dt_onbs[,pgvr]
ckck
max(ckck)
min(ckck)
mean(ckck)

#co-ancestor
dt_onbs= dt[program%in%c("OS_10","OS_30","OS_50","OP_10","OP_30","OP_50","CP","CS")&nGeneration==5,]
ckck= (dt_onbs$nCo - 150)/150
min(ckck)
max(ckck)
dt_onbs= dt[program%in%c("OS_10","OS_30","OS_50","OP_10","OP_30","OP_50")&nGeneration%in%c(6:20),]
dt_onbs[program=="OS_10",nCo]

#shrot-term inbreeding rate
sinbrate= dt[program%in%c("OS_10","OS_30","OS_50","OP_10","OP_30","OP_50","OS100","OS0","OP100","OP0"),.(program,prates)] %>% unique()

#long-term inbreeding rate
linbrate= dt[program%in%c("OS_10","OS_30","OS_50","OP_10","OP_30","OP_50","OS100","OS0","OP100","OP0"),.(program,prate)] %>% unique()
#

#relationship between MP individuals and NP individuals
#dt_plot= dt[program%in%c("OP_10","OP_30","OP_50")&nGeneration%in%c(2:20),]
#dt_plot= dt[program%in%c("OS_10","OS_30","OS_50")&nGeneration%in%c(2:20),]
#dt_plot[program%flike%"C",breeding_sys:="CNBS"]
#dt_plot[!program%flike%"C",breeding_sys:="ONBS"]
#dt_plot= dt_plot[!gradients=="non-overlap",]
#dt_plot[gradients=="CNBS",gradients:=NA]
#dt_plot[nGeneration==0,inb:=0]
#dt_plot[nGeneration==0,inbp:=0]
#dt_plot[nGeneration==0,ld:=0]
#dt_plot[gradients=="ONBS (10%)",gradients:="10"]
#dt_plot[gradients=="ONBS (30%)",gradients:="30"]
#dt_plot[gradients=="ONBS (50%)",gradients:="50"]
#P2<- ggplot(data = dt_plot,aes(x=nGeneration,y=prel,group=program,shape = gradients,linetype=breeding_sys))+
 # geom_point()+
  #geom_line()+
  #scale_shape_manual(values = c(0,3,15,1),na.translate=FALSE)+
  #xlab("Generation")+
#  ylab("inbreeding")+
#  theme_bw()+
#  theme(panel.grid.major = element_line(colour = NA),
 #       panel.background = element_rect(fill = "transparent",colour = NA),
#        plot.background = element_rect(fill = "transparent",colour = NA),
#        panel.grid.minor = element_blank(),
 #       text = element_text(family = "STXihei"),
#        legend.position = "right",
#        legend.background = element_rect(colour = "black"))+
  #geom_errorbar(aes(ymin=inbp-inbpse,
  #                 ymax=inbp+inbpse),
  #            width=0.05)+
#  scale_x_continuous(limits = c(2,20),breaks = seq(2,20,1))+labs(linetype="Breeding system",shape="Proportion of\nthe introduced\nMP individuals (%)")
#P1
#P2

#figure <- ggarrange(P1 +ylab("Genomic relationship")+ylim(0,0.4) , P2 + rremove("ylab")+ylim(0,0.4) , # remove axis labels from plots
#                    labels = c("a","b"),
 #                   ncol = 2, nrow =1,
#                    widths = c(1,1),
 #                   common.legend = TRUE, legend = "right",
#                    align = "v", 
 #                   font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
#figure
#ggsave("prel.tiff", figure , width = 10, height = 5, dpi = 300)
