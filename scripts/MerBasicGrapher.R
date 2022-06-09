#maybe okay 012220


#rkbypass expression
fWT_Data<-csvSheller("EWT.csv","filteredfiles", targetColumn = "polfire")
fM4_Data<-csvSheller("EWT.csv","filteredfiles", targetColumn = "polfire")
fMD_Data<-csvSheller("EWT.csv","filteredfiles", targetColumn = "polfire")
fTM_Data<-csvSheller("EWT.csv","filteredfiles", targetColumn = "polfire")
f9I_Data<-csvSheller("EWT.csv","filteredfiles", targetColumn = "polfire")
f5I_Data<-csvSheller("EWT.csv","filteredfiles", targetColumn = "polfire")
fMI_Data<-csvSheller("EWT.csv","filteredfiles", targetColumn = "polfire")
fCM_Data<-csvSheller("EWT.csv","filteredfiles", targetColumn = "polfire")
f5M_Data<-csvSheller("EWT.csv","filteredfiles", targetColumn = "polfire")



mergedDynamicData<-do.call("cbind", list(fWT_Data,fM4_Data,fMD_Data,fTM_Data,f9I_Data,f5I_Data,fMI_Data,fCM_Data,f5M_Data))


colnames(mergedDynamicData)<-c("FWTL1","FWTL2","FWTL3","FWTL4","FWTL5","FWTL6","FWTL7",
                               "FM4L1","FM4L2","FM4L3","FM4L4","FM4L5","FM4L6","FM4L7",
                               "FMDL1","FMDL2","FMDL3","FMDL4","FMDL5","FMDL6","FMDL7",
                               "FTML1","FTML2","FTML3","FTML4","FTML5","FTML6","FTML7",
                               "F9iL1","F9iL2","F9iL3","F9iL4","F9iL5","F9iL6","F9iL7",
                               "F5iL1","F5iL2","F5iL3","F5iL4","F5iL5","F5iL6","F5iL7",
                               "FMiL1","FMiL2","FMiL3","FMiL4","FMiL5","FMiL6","FMiL7",
                               "FCML1","FCML2","FCML3","FCML4","FCML5","FCML6","FCML7",
                               "F5ML1","F5ML2","F5ML3","F5ML4","F5ML5","F5ML6","F5ML7"
                               )





mergedDynamicData[is.na(mergedDynamicData)] = 0





markerexpgenotype<-c(rep("1: wildtype",5),rep("2: 970m4",5),rep("3: 970md",5),rep("4: tm",5),rep("5: 970i",5),rep("6: 950i",5),rep("7: 970M4i",5),rep("8: clv3 null",5),rep("9: 950m",5))
simgenotype<-c(rep("1: wildtype",140),rep("2: 970m4",140),rep("3: 970md",140),rep("4: tm",140),rep("5: 970i",140),rep("6: 950i",140),rep("7: 970M4i",140),rep("8: clv3 null",140),rep("9: 950m",140))



simgenotype<-melt(simgenotype)
colnames(simgenotype)<-c("simgenotype")
mergedDynamicData<-melt(mergedDynamicData)
mergedDynamicData<-cbind(mergedDynamicData,simgenotype)



ggplot(mergedDynamicData, aes(x=variable, y=value, colour=simgenotype)) + ylab("polybind") + xlab("genotype")+ ggtitle("published directpolybind test 3 direct sim time for wt")  +
  theme(plot.title = element_text(hjust = 0.5))+ geom_jitter(shape=16, position=position_jitter(0.2))+coord_cartesian(ylim=c(0,5400))





