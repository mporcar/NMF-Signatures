# remove.packages("data.table")                         # First remove the current version
# install.packages("data.table", type = "source",repos = "http://Rdatatable.github.io/data.table") #
# install.packages("data.table")
country="Spain"
if(country == "Spain" )name_file<-"NQ-154707_v3_ES" #Spain
if(country == "Argentina" )name_file<-"NQ-154707_v3_AR" #Argentina
if(country == "Mexico" )name_file<-"NQ-154707_v3_MX" #Mexico
if(country == "Portugal" )name_file<-"NQ-154707_v3_PT" #Portugal
if(country == "Colombia" )name_file <- "NQ-154707_v3_CO" #Colombia
if(country == "Peru" )name_file <- "NQ-154707_v3_PE" #Peru
if(country == "Chile" )name_file <-"NQ-154707_v3_CL" # Chile
if(country == "Brasil" )name_file <-"NQ-154707_v3_BR" #Brasil
dt<-fread(paste0("data/",name_file,".csv"), na.strings = "NULL") 
dz <- fread("data/z_funnel_captacion.csv")
dt<-merge(dt,dz[,.(panelist,CampaignType)],by="panelist",all.x = T)
rm(dz)
#filter by the first recruitment survey
dt<-dt[survey_1_sm2a_nombre=="global.recruitment@netquest.com"]

#Prepare data------------------
setkey(dt,panelist)
dt[, p_sexo:=as.factor(p_sexo)]
dt[, p_edad:=as.factor(p_edad)]
dt[, claseSocial:=as.factor(claseSocial)]
dt[, regio_geografica:=as.factor(regio_geografica)]
dt[, survey_1_sm2a_nombre := as.factor(survey_1_sm2a_nombre)]
dt[, survey_2_sm2a_nombre := as.factor(survey_2_sm2a_nombre)]
dt[, survey_3_sm2a_nombre := as.factor(survey_3_sm2a_nombre)]
dt[, survey_4_sm2a_nombre := as.factor(survey_4_sm2a_nombre)]
dt[, survey_5_sm2a_nombre := as.factor(survey_5_sm2a_nombre)]
dt[, survey_6_sm2a_nombre := as.factor(survey_6_sm2a_nombre)]

##### The missing data is imputed to avoid programing problems
dt[is.na(survey_1_inv_date),survey_1_inv_date:="1900-01-01 00:00:00"]
dt[is.na(survey_2_inv_date),survey_2_inv_date:="1900-01-01 00:00:00"]
dt[is.na(survey_3_inv_date),survey_3_inv_date:="1900-01-01 00:00:00"]
dt[is.na(survey_4_inv_date),survey_4_inv_date:="1900-01-01 00:00:00"]
dt[is.na(survey_5_inv_date),survey_5_inv_date:="1900-01-01 00:00:00"]
dt[is.na(survey_6_inv_date),survey_6_inv_date:="1900-01-01 00:00:00"]
dt[is.na(survey_1_participation_date),survey_1_participation_date:="1900-01-01 00:00:00"]
dt[is.na(survey_2_participation_date),survey_2_participation_date:="1900-01-01 00:00:00"]
dt[is.na(survey_3_participation_date),survey_3_participation_date:="1900-01-01 00:00:00"]
dt[is.na(survey_4_participation_date),survey_4_participation_date:="1900-01-01 00:00:00"]
dt[is.na(survey_5_participation_date),survey_5_participation_date:="1900-01-01 00:00:00"]
#### In survey 6 there are dates with just ""
dt[is.na(survey_6_participation_date)| survey_6_participation_date == "",survey_6_participation_date:="1900-01-01 00:00:00"]
dt[is.na(p_fecha_registro) | p_fecha_registro == "",p_fecha_registro:="1900-01-01 00:00:00"]
dt=dt[survey_1==1 & survey_2==1 & survey_3==1 & survey_4==1 & survey_5==1 & survey_6==1 ]
dt[,delay_part1:=difftime(dt$survey_1_participation_date,dt$survey_1_inv_date,units = "secs")]
dt[,delay_part2:=difftime(dt$survey_2_participation_date,dt$survey_2_inv_date,units = "secs")]
dt[,delay_part3:=difftime(dt$survey_3_participation_date,dt$survey_3_inv_date,units = "secs")]
dt[,delay_part4:=difftime(dt$survey_4_participation_date,dt$survey_4_inv_date,units = "secs")]
dt[,delay_part5:=difftime(dt$survey_5_participation_date,dt$survey_5_inv_date,units = "secs")]
dt[,delay_part6:=difftime(dt$survey_6_participation_date,dt$survey_6_inv_date,units = "secs")]

dt[, time_inv1 := tstrsplit(survey_1_inv_date, ' ')[2]]
dt[, time_inv2 := tstrsplit(survey_2_inv_date, ' ')[2]]
dt[, time_inv3 := tstrsplit(survey_3_inv_date, ' ')[2]]
dt[, time_inv4 := tstrsplit(survey_4_inv_date, ' ')[2]]
dt[, time_inv5 := tstrsplit(survey_5_inv_date, ' ')[2]]
dt[, time_inv6 := tstrsplit(survey_6_inv_date, ' ')[2]]


# delay2=as.numeric(delay2[delay2>0])
# delay3=as.numeric(delay3[delay3>0])
# delay4=as.numeric(delay4[delay4>0])
# delay5=as.numeric(delay6[delay5>0])
# delay6=as.numeric(delay6[delay6>0])
# delay1=round(delay1/3600,4)
# delay2=round(delay2/3600,4)
# delay3=round(delay3/3600,4)
# delay4=round(delay4/3600,4)
# delay5=round(delay5/3600,4)
# delay6=round(delay6/3600,4)
# delay1=as.numeric(delay3[delay1<10*24])
# delay2=as.numeric(delay4[delay2<10*24])
# delay3=as.numeric(delay3[delay3<10*24])
# delay4=as.numeric(delay4[delay4<10*24])
# delay5=as.numeric(delay5[delay5<10*24])
# delay6=as.numeric(delay6[delay6<10*24])
# plot(density(delay3,na.rm= TRUE),xlim=c(0,150))
# lines(density(delay4,na.rm= TRUE))
# lines(density(delay5,na.rm= TRUE))
# lines(density(delay6,na.rm= TRUE))
# data=as.data.table(rbind(cbind.data.frame(delay=delay1,survey="survey1"),
#                          cbind.data.frame(delay=delay2,survey="survey2"),
#                          cbind.data.frame(delay=delay3,survey="survey3"),
#                          cbind.data.frame(delay=delay4,survey="survey4"),
#                          cbind.data.frame(delay=delay5,survey="survey5"),
#                          cbind.data.frame(delay=delay6,survey="survey6")))
# 
# require(gridExtra)
# plot1<-ggplot(data[survey %in% c("survey1","survey3")], aes(delay, fill = survey, colour = survey)) +
#   geom_density(alpha = 0.1) +
#   xlim(0, 150)+
#   ylim(0,0.07)
# plot2<-ggplot(data[survey %in% c("survey2","survey4")], aes(delay, fill = survey, colour = survey)) +
#   geom_density(alpha = 0.1) +
#   xlim(0, 150)+
#   ylim(0,0.07)
# plot3<-ggplot(data[survey %in% c("survey5","survey6")], aes(delay, fill = survey, colour = survey)) +
#   geom_density(alpha = 0.1) +
#   xlim(0, 150)+
#   ylim(0,0.07)
# grid.arrange(plot1, plot2,plot3,ncol=3)

dt[, time_par1 := tstrsplit(survey_1_participation_date, ' ')[2]]
dt[, time_par2 := tstrsplit(survey_2_participation_date, ' ')[2]]
dt[, time_par3 := tstrsplit(survey_3_participation_date, ' ')[2]]
dt[, time_par4 := tstrsplit(survey_4_participation_date, ' ')[2]]
dt[, time_par5 := tstrsplit(survey_5_participation_date, ' ')[2]]
dt[, time_par6 := tstrsplit(survey_6_participation_date, ' ')[2]]

col <- c("survey_1_participation_date","survey_2_participation_date","survey_3_participation_date",
         "survey_4_participation_date","survey_5_participation_date","survey_6_participation_date",
         "survey_1_inv_date","survey_2_inv_date","survey_3_inv_date","survey_4_inv_date",
         "survey_5_inv_date","survey_6_inv_date","p_fecha_registro"
)
dt[, (col) := NULL]



#transform hours to Itime
dt[, time_inv1 := as.ITime(time_inv1,format ="%H:%M:%S")]
dt[, time_inv2 := as.ITime(time_inv2,format ="%H:%M:%S")]
dt[, time_inv3 := as.ITime(time_inv3,format ="%H:%M:%S")]
dt[, time_inv4 := as.ITime(time_inv4,format ="%H:%M:%S")]
dt[, time_inv5 := as.ITime(time_inv5,format ="%H:%M:%S")]
dt[, time_inv6 := as.ITime(time_inv6,format ="%H:%M:%S")]



id <- c("panelist", "p_sexo", "p_edad" , "claseSocial", "regio_geografica", "CampaignType")

survey <- paste0(rep("survey_",6),1:6)
name_survey <-paste0(rep("survey_",6),1:6,rep("_sm2a_nombre",6))
col <- c(id,survey)
dt1<-melt(dt[,(col), with = FALSE], id = id)
dt1[, Nsurvey := gsub("survey_","",variable)]

#transform the data for the delays of the participation
col <-c(id,paste0(rep("delay_part",6),1:6)) 
dtdiff<-melt(dt[,(col), with = FALSE], id = id)
dtdiff[,`:=` (Nsurvey=gsub("delay_part","",variable),
              variable = NULL,
              diffpar = value,
              value = NULL)]


#time in hours
col <- c(id, paste0(rep("time_inv",6),1:6))
dthour <-melt(dt[,(col), with = FALSE], id = id)
dthour[,`:=` (Nsurvey=gsub("time_inv","",variable),
              variable = NULL,
              hour_inv = value,
              value = NULL)]

# we merge the datasets
dataset<-merge(dthour, dtdiff, by = c(id,"Nsurvey"), all.x = TRUE)
dataset=dataset[diffpar>0]

dataset[,Hinv:=hour(hour_inv)]

cuts<-c(6,12,18)
lab <-c("0-6","6-12","12-18","18-24")
#order the dataset findInterval
dataset<-dataset[order(Hinv)]
#finding the intervals 
dataset[,hinv:=findInterval(Hinv,cuts)]
#put the labels in a new column
dataset[,interval:= lab[1+hinv]]
#change the column names
dataset[,`:=`(hinv = NULL)]
dataset[,delay_part:=round(diffpar/3600,4)]
