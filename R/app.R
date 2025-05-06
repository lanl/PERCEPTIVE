#' @export
PERCEPTIVE<-function(){
library(purrr)
library(shiny)
library(reactable)
library(shinyFiles)
library(DT)
library(htmltools)
library(fontawesome)
library(reactablefmtr)
library(fresh)
library(shinydashboard)
library(shinydashboardPlus)
library(ggplot2)





server <- function(input, output, session){



  character_match <<- function(string, edit_string, match, drop = NA){
    # convert to array
    string = strsplit(string, "")[[1]]
    edit_string = strsplit(edit_string, "")[[1]]

    if(!is.na(drop)){
      edit_string = edit_string[edit_string != drop]
    }

    if(length(string) != length(edit_string)){
      stop("string and edit_string are different lengths")
    }

    output = rep("_", length(edit_string))
    is_match = edit_string == match
    output[is_match] = string[is_match]

    output = paste0(output, collapse = "")
    return(output)
  }


  similarity<<- function(histone, species){
    if (histone=="H1")
    {organismtable<-H1
    pos<-2
    ac<-3
    ph<-4
    me<-5
    ub<-6
    UP<-7}
    if (histone=="H2A")
    {organismtable<-H2A
    pos<-9
    ac<-10
    ph<-11
    me<-12
    ub<-13
    UP<-14}
    if (histone=="H2B")
    {organismtable<-H2B
    pos<-16
    ac<-17
    ph<-18
    me<-19
    ub<-20
    UP<-21}
    if (histone=="H3")
    {organismtable<-H3
    pos<-23
    ac<-24
    ph<-25
    me<-26
    ub<-27
    UP<-28}
    if (histone=="H4")
    {organismtable<-H4
    pos<-30
    ac<-31
    ph<-32
    me<-33
    ub<-34
    UP<-35}
    if (species=="sc")
    {modeltable<-sc}
    if (species=="hh")
    {modeltable<-hh}

    outtable<-rbind(modeltable[pos,2], modeltable[ac,2], modeltable[ph,2], modeltable[me,2], modeltable[ub,2])
    names<-c(paste0(modeltable[pos-1,1]," ", modeltable[pos,1] ), modeltable[ac,1], modeltable[ph,1], modeltable[me,1], modeltable[ub,1])
    for(i in 1:dim(organismtable)[1])
    {

      if(nchar(modeltable[pos,2])>nchar(organismtable[i,2]))
      {
        s1 = modeltable[pos,2]
        s2 = organismtable[i,2]
      }else{
        s2 = modeltable[pos,2]
        s1 = organismtable[i,2]
      }


      out = adist(s1, s2, counts = TRUE)
      edit_string = drop(attr(out, "trafos"))
      match<-character_match(s1, edit_string, "M", "I")
      sub<-character_match(s1, edit_string, "S", "I")
      delinsert<-character_match(s1, edit_string, "D", "I")

      outtable<-rbind(outtable, organismtable[i,2], match,sub,delinsert)
      names<-c(names, paste0("Protein sequence for transcript ID: ", organismtable[i,1]), "Match", "Substitutions", "Insertion/Deletions" )

    }
    outtable<-cbind(names, outtable)
    bar<-vector()
    num<-vector()
    for (c in 1:floor(nchar(paste(outtable[1,2]))/10))
    {


      if (c>10)
      {
        bar<-c(bar, c("&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp", "|"))
        num<-c(num, c("&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp",c*10))
      } else
      {
        bar<-c(bar, c("&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","|"))
        num<-c(num, c("&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp","&nbsp",c*10))
      }
    }





    outtable<-rbind(c("Residue #", paste(num, collapse="")),c("", paste(bar, collapse="")), outtable)


    head1<-paste(outtable[4,1], "<br>", outtable[5,1], "<br>", outtable[6,1], "<br>", outtable[7,1], "<br>", outtable[1,1], "<br>", outtable[2,1], "<br>", outtable[3,1], "<br>")
    head2<-paste(outtable[4,2], "<br>", outtable[5,2], "<br>", outtable[6,2], "<br>", outtable[7,2], "<br>", outtable[1,2], "<br>", outtable[2,2], "<br>", outtable[3,2], "<br>")
    outtable<-outtable[-(1:7),]
    colnames(outtable)[1]<-head1
    colnames(outtable)[2]<-head2
    outtable<<-outtable
  }






  volumes = getVolumes() ()
  observe({
    shinyDirChoose(input, 'folder', roots=volumes)

  }

  )
  output$dir <-dir<- renderText(parseDirPath(roots = volumes, input$folder)
  )


  output$source <- renderText("You must select a data source for analysis!")
  output$description <- renderText("InterProScan HMM Hits")

  observeEvent(input$startanalysis,
               {

                    if(isTRUE(dir.exists(parseDirPath(roots = volumes, input$folder))))
                    {dir<-parseDirPath(roots = volumes, input$folder)
                    output$dir <-dir<- renderText(parseDirPath(roots = volumes, input$folder)
                    )
                    blastp<<-blastp<-read.table(paste0(dir(),"/blastbpstricteval.table"), sep="\t", header=FALSE)
                    Interpro<<-Interpro<-read.table(paste0(dir(),"/InterproIdentified.txt"), sep="\t", header=FALSE)
                    protname<<-protname<-read.delim(paste0(dir(),"/Speciesandprotein.txt"), header=FALSE)
                    IP_Table<-read.csv(paste0(dir(),"/ListofINTERPRONUMBERS.csv"), header=TRUE)
                    blastp[,2]<<-protname
                    Interprolist<-unique(Interpro[,12])
                    counts<-as.data.frame(table(Interpro[,12]))
                    totalhits<-as.numeric(counts[,2])
                    counts<-as.character(droplevels(counts[,1]))
                    IP_list<-as.character(IP_Table[,1])
                    matching<-0
                    IP_table_w_hits<-cbind(IP_Table[,1], rep(NA, dim(IP_Table)[1]), IP_Table[,2:7])


                    for (i in 1:length(counts))
                    {
                      x<-which(IP_list[]==counts[i])
                      if(!is_empty(x))
                      {
                        IP_table_w_hits[x,2]<-totalhits[i]
                      } else

                        rm(x)
                    }
                    IP_table_w_hits<-IP_table_w_hits[,-6]

                    IP_table_w_hits<<-IP_table_w_hits<-IP_table_w_hits[order (decreasing = TRUE, IP_table_w_hits[,2]),]

                    IP_table_w_hits_details<<-IP_table_w_hits_details<-cbind(IP_table_w_hits, details = NA)

                    for (y in 1:dim(IP_table_w_hits_details)[1])
                    {
                      if (!is.na(IP_table_w_hits_details[y,2]))
                      {
                        IP_table_w_hits_details[y,8]<- as.character(tags$div(tags$button("Show details")))}
                    }

                    IP_table_w_hits_details<-IP_table_w_hits_details

                    AT<-as.data.frame(read.csv(paste0(dir(),"/associationTable.csv"), header = FALSE))
                    icontable<-as.vector(IP_table_w_hits_details[,1])
                    colortable<-as.vector(IP_table_w_hits_details[,1])
                  for(a in 1:dim(IP_table_w_hits_details)[1])
                  {
                    if(length(which(IP_table_w_hits_details[a,1]==AT[,2]))>0)
                    {
                  icontable[a]<-AT[which(IP_table_w_hits_details[a,1]==AT[,2]),3]
                  colortable[a]<-AT[which(IP_table_w_hits_details[a,1]==AT[,2]),4]

                    }
                  }
                    icontable[grepl("IPR", icontable)]<-"question"
                    colortable[grepl("IPR", colortable)]<-"black"
                    IP_table_w_hits_details<-cbind(IP_table_w_hits_details, icontable, colortable)
                    colnames(IP_table_w_hits_details)[1]<- c("Accession")
                    colnames(IP_table_w_hits_details)[2]<- c("HMM Protein Homology Hits")
                    colnames(IP_table_w_hits_details)[3]<- c("NAME")
                    colnames(IP_table_w_hits_details)[4]<- c("HMM Database")
                    colnames(IP_table_w_hits_details)[5]<- c("Class")
                    colnames(IP_table_w_hits_details)[6]<- c("Integrated Signatures")
                    colnames(IP_table_w_hits_details)[7]<- c("Go Terms")
                    colnames(IP_table_w_hits_details)[8]<- c("details")
                    colnames(IP_table_w_hits_details)[9]<- c("icontable")
                    colnames(IP_table_w_hits_details)[10]<- c("colortable")
                    IP_table_w_hits_details<<-IP_table_w_hits_details



                    output$IPhits<-renderReactable({reactable(IP_table_w_hits_details, pagination = FALSE, searchable = TRUE,
                              columns = list(
                                Accession = colDef(width=90),
                                "HMM Protein Homology Hits" = colDef(width=120, na="0",
                                                                     cell = icon_sets(icon_ref = "icontable", icon_color_ref ="colortable", IP_table_w_hits_details),
                                                                     style = function(value) {
                                                                       if (is.na(value)) {
                                                                         color <- "lightgrey"
                                                                       } else
                                                                       {
                                                                         color <- "darkblue"
                                                                       }
                                                                       list(color = color, fontWeight = "bold")
                                                                     }
                                ),
                                NAME = colDef(minWidth=200),
                                'HMM Database' = colDef(width=100),
                                Class = colDef(minWidth=135),
                                'Integrated Signatures' = colDef(minWidth=125),
                                'Go Terms' = colDef(minWidth=175),
                                icontable = colDef(show=FALSE),
                                colortable = colDef(show=FALSE),
                                details = colDef(
                                  name = "",
                                  sortable = TRUE,
                                  html = TRUE)

                              ),
                              onClick = JS("function(rowInfo, column) {
    // Only handle click events on the 'details' column
    if (column.id !== 'details') {
      return
    }

    // Send the click event to Shiny, which will be available in input$show_details
    // Note that the row index starts at 0 in JavaScript, so we add 1
    if (window.Shiny) {
      Shiny.setInputValue('show_details', { index: rowInfo.index + 1 }, { priority: 'event' })
    }
  }")



                    )})

                    output$source <- renderText(parseDirPath(roots = volumes, input$folder))
                    H1<-as.data.frame(0)
                    H1<-tryCatch(read.table(paste0(dir(),"/H2A"), header=FALSE), error=function(e) {NULL})
                    if(is.numeric(dim(H1)))
                    {
                    header<-which(grepl(">", H1[,1])==TRUE)
                  if(length(header)>1)
                  {
                    x<-rbind(paste0(c(H1[((header[1]+1):(header[1+1]+-1)),1]), collapse=""))
                    for (t in 2:(length(header)-1))
                    {
                      x<-rbind(x,paste0(c(H1[((header[1]+1):(header[1+1]+-1)),1]), collapse=""))
                    }
                    x<-rbind(x,paste0(c(H1[(header[length(header)]+1):dim(H1)[1],1]), collapse=""))
                  }else
                  {
                    x<-rbind(paste0(c(H1[((header[1]+1):dim(H1[1])[1]),1]), collapse=""))
                  }
                    H2A<-x
                    H2A<<-cbind(H1[header,],H2A)
                    rm(H1)
                    }else
                    {
                    H2A<<-matrix("NULL", 1,2)
                    }
                    H1<-as.data.frame(0)
                    H1<-tryCatch(read.table(paste0(dir(),"/H2B"), header=FALSE), error=function(e) {NULL})
                    if(is.numeric(dim(H1)))
                    {

                    header<-which(grepl(">", H1[,1])==TRUE)
                    if(length(header)>1)
                    {
                      x<-rbind(paste0(c(H1[((header[1]+1):(header[1+1]+-1)),1]), collapse=""))
                      for (t in 2:(length(header)-1))
                      {
                        x<-rbind(x,paste0(c(H1[((header[1]+1):(header[1+1]+-1)),1]), collapse=""))
                      }
                      x<-rbind(x,paste0(c(H1[(header[length(header)]+1):dim(H1)[1],1]), collapse=""))
                    }else
                    {
                      x<-rbind(paste0(c(H1[((header[1]+1):dim(H1[1])[1]),1]), collapse=""))
                    }
                    H2B<-x
                    H2B<<-cbind(H1[header,],H2B)
                    rm(H1)
                    }else{H2B<<-matrix("NULL", 1,2)

                    }
                    H1<-as.data.frame(0)
                    H1<-tryCatch(read.table(paste0(dir(),"/H3"), header=FALSE), error=function(e) {NULL})
                    if(is.numeric(dim(H1)))
                    {
                    header<-which(grepl(">", H1[,1])==TRUE)
                    if(length(header)>1)
                    {
                      x<-rbind(paste0(c(H1[((header[1]+1):(header[1+1]+-1)),1]), collapse=""))
                      for (t in 2:(length(header)-1))
                      {
                        x<-rbind(x,paste0(c(H1[((header[1]+1):(header[1+1]+-1)),1]), collapse=""))
                      }
                      x<-rbind(x,paste0(c(H1[(header[length(header)]+1):dim(H1)[1],1]), collapse=""))
                    }else
                    {
                      x<-rbind(paste0(c(H1[((header[1]+1):dim(H1[1])[1]),1]), collapse=""))
                    }
                    H3<-x
                    H3<<-cbind(H1[header,],H3)
                    rm(H1)
                    }else{H3<<-matrix("NULL", 1,2)

                    }
                    H1<-as.data.frame(0)
                    H1<-tryCatch(read.table(paste0(dir(),"/H4"), header=FALSE), error=function(e) {NULL})
                    if(is.numeric(dim(H1)))
                    {
                    header<-which(grepl(">", H1[,1])==TRUE)
                    if(length(header)>1)
                    {
                      x<-rbind(paste0(c(H1[((header[1]+1):(header[1+1]+-1)),1]), collapse=""))
                      for (t in 2:(length(header)-1))
                      {
                        x<-rbind(x,paste0(c(H1[((header[1]+1):(header[1+1]+-1)),1]), collapse=""))
                      }
                      x<-rbind(x,paste0(c(H1[(header[length(header)]+1):dim(H1)[1],1]), collapse=""))
                    }else
                    {
                      x<-rbind(paste0(c(H1[((header[1]+1):dim(H1[1])[1]),1]), collapse=""))
                    }
                    H4<-x
                    H4<<-cbind(H1[header,],H4)
                    rm(H1)
                    }else{H4<<-matrix("NULL", 1,2)

                    }
                    H1<-as.data.frame(0)
                    H1<-tryCatch(read.table(paste0(dir(),"/H1"), header=FALSE), error=function(e) {NULL})
                    if(is.numeric(dim(H1)))
                    {

                    header<-which(grepl(">", H1[,1])==TRUE)
                        if(length(header)>1)
                        {
                          x<-rbind(paste0(c(H1[((header[1]+1):(header[1+1]+-1)),1]), collapse=""))
                          for (t in 2:(length(header)-1))
                          {
                            x<-rbind(x,paste0(c(H1[((header[1]+1):(header[1+1]+-1)),1]), collapse=""))
                          }
                          x<-rbind(x,paste0(c(H1[(header[length(header)]+1):dim(H1)[1],1]), collapse=""))
                        }else
                        {
                          x<-rbind(paste0(c(H1[((header[1]+1):dim(H1[1])[1]),1]), collapse=""))
                        }
                    H1.0<-x
                    H1.0<-cbind(H1[header,],H1.0)
                    H1<<-H1.0
                    rm (H1.0)
                    rm(x)
                    }else{
                      H1<<-matrix("NULL",1,2)

                    }

                    hh<-read.csv(paste0(dir(),"/humanhistones.csv"))
                    hh<-hh[,-1]
                    hh<-apply(hh, 1, paste, collapse="")
                    hh<<-cbind(c("Histone: H1.0", "Sequence", "acetylation", "phosphorylation", "methylation", "ubiquitylation", "UniProtID:",
                                "Histone: H2A.1", "Sequence", "acetylation", "phosphorylation", "methylation", "ubiquitylation", "UniProtID:",
                                "Histone: H2B.1", "Sequence", "acetylation", "phosphorylation", "methylation", "ubiquitylation", "UniProtID:",
                                "Histone: H3.1", "Sequence", "acetylation", "phosphorylation", "methylation", "ubiquitylation", "UniProtID:",
                                "Histone: H4", "Sequence", "acetylation", "phosphorylation", "methylation", "ubiquitylation", "UniProtID:"),hh)



                    sc<-read.csv(paste0(dir(),"/scerhistones.csv"))
                    sc<-sc[,-1]
                    sc<-apply(sc, 1, paste, collapse="")
                    sc<<-cbind(c("Histone: H1", "Sequence", "acetylation", "phosphorylation", "methylation", "ubiquitylation", "UniProtID:",
                                "Histone: H2A.1", "Sequence", "acetylation", "phosphorylation", "methylation", "ubiquitylation", "UniProtID:",
                                "Histone: H2B.1", "Sequence", "acetylation", "phosphorylation", "methylation", "ubiquitylation", "UniProtID:",
                                "Histone: H3", "Sequence", "acetylation", "phosphorylation", "methylation", "ubiquitylation", "UniProtID:",
                                "Histone: H4", "Sequence", "acetylation", "phosphorylation", "methylation", "ubiquitylation", "UniProtID:"),sc)
                    sp<<-"hh"
                    his<<-"H1"
                    ### generate predictive matrix ###

                    H2BHH<-similarity(histone = "H2B", species = "hh")
                    H3HH<-similarity(histone = "H3", species = "hh")
                    H4HH<-similarity(histone = "H4", species = "hh")
                    H2BSC<-similarity(histone = "H2B", species = "sc")
                    H3SC<-similarity(histone = "H3", species = "sc")
                    H4SC<-similarity(histone = "H4", species = "sc")
                    ###H2bub###

                    H2BHHhits<-0
                    l<-2
                    for(f in 1:(dim(H2BHH)[1]/4))
                    {

                    if(unlist(strsplit(H2BHH[(l),2],""))[121]=="K")
                    {
                    H2BHHhits<-1
                    }
                    l<-l+4
                    }

                    H2BSChits<-0
                    l<-2
                    for(f in 1:(dim(H2BSC)[1]/4))
                    {

                      if(unlist(strsplit(H2BSC[(l),2],""))[124]=="K")
                      {
                        H2BSChits<-1
                      }
                      l<-l+4
                    }
                    ###H3Marks###
                    H3HHk4<-0
                    H3HHk36<-0
                    H3HHk79<-0
                    H3HHk9<-0
                    H3HHk27<-0
                    H3HHs10<-0
                    H3HHk14<-0
                    H3HHk18<-0
                    H3HHk23<-0
                    H3HHk56<-0
                    l<-2
                    for(f in 1:(dim(H3HH)[1]/4))
                    {

                      if(unlist(strsplit(H3HH[(l),2],""))[5]=="K")
                      {
                        H3HHk4<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[37]=="K")
                      {
                        H3HHk36<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[80]=="K")
                      {
                        H3HHk79<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[10]=="K")
                      {
                        H3HHk9<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[28]=="K")
                      {
                        H3HHk27<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[11]=="S")
                      {
                        H3HHs10<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[15]=="K")
                      {
                        H3HHk14<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[19]=="K")
                      {
                        H3HHk18<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[24]=="K")
                      {
                        H3HHk23<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[57]=="K")
                      {
                        H3HHk56<-1
                      }
                      l<-l+4
                    }

                    H3SCk4<-0
                    H3SCk36<-0
                    H3SCk79<-0
                    H3SCk9<-0
                    H3SCk27<-0
                    H3SCs10<-0
                    H3SCk14<-0
                    H3SCk18<-0
                    H3SCk23<-0
                    H3SCk56<-0
                    l<-2
                    for(f in 1:(dim(H3SC)[1]/4))
                    {

                      if(unlist(strsplit(H3SC[(l),2],""))[5]=="K")
                      {
                        H3SCk4<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[37]=="K")
                      {
                        H3SCk36<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[80]=="K")
                      {
                        H3SCk79<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[10]=="K")
                      {
                        H3SCk9<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[28]=="K")
                      {
                        H3SCk27<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[11]=="S")
                      {
                        H3SCs10<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[15]=="K")
                      {
                        H3SCk14<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[19]=="K")
                      {
                        H3SCk18<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[24]=="K")
                      {
                        H3SCk23<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[57]=="K")
                      {
                        H3SCk56<-1
                      }
                      l<-l+4
                    }
                    ###H4 MARKS ###
                    H4HHk5<-0
                    H4HHk8<-0
                    H4HHk12<-0
                    H4HHk14<-0
                    H4HHk16<-0
                    H4HHk20<-0

                    l<-2
                    for(f in 1:(dim(H4HH)[1]/4))
                    {

                      if(unlist(strsplit(H4HH[(l),2],""))[6]=="K")
                      {
                        H4HHk5<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[9]=="K")
                      {
                        H4HHk8<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[13]=="K")
                      {
                        H4HHk12<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[15]=="K")
                      {
                        H4HHk14<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[17]=="K")
                      {
                        H4HHk16<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[21]=="K")
                      {
                        H4HHk20<-1
                      }
                      l<-l+4
                    }

                    H4SCk5<-0
                    H4SCk8<-0
                    H4SCk12<-0
                    H4SCk14<-0
                    H4SCk16<-0
                    H4SCk20<-0
                    l<-2
                    for(f in 1:(dim(H4SC)[1]/4))
                    {

                      if(unlist(strsplit(H4SC[(l),2],""))[6]=="K")
                      {
                        H4SCk5<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[9]=="K")
                      {
                        H4SCk8<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[13]=="K")
                      {
                        H4SCk12<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[15]=="K")
                      {
                        H4SCk14<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[17]=="K")
                      {
                        H4SCk16<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[21]=="K")
                      {
                        H4SCk20<-1
                      }
                      l<-l+4
                    }

                    ### left adjacent residue ###

                    l<-2
                    H2BHHlefthits<-0
                    for(f in 1:(dim(H2BHH)[1]/4))
                    {

                      if(unlist(strsplit(H2BHH[(l),2],""))[120]=="K")
                      {
                        H2BHHlefthits<-1
                      }
                      l<-l+4
                    }

                    H2BSClefthits<-0
                    l<-2
                    for(f in 1:(dim(H2BSC)[1]/4))
                    {

                      if(unlist(strsplit(H2BSC[(l),2],""))[123]=="K")
                      {
                        H2BSClefthits<-1
                      }
                      l<-l+4
                    }
                    ###H3Marks###
                    H3HHleftk4<-0
                    H3HHleftk36<-0
                    H3HHleftk79<-0
                    H3HHleftk9<-0
                    H3HHleftk27<-0
                    H3HHlefts10<-0
                    H3HHleftk14<-0
                    H3HHleftk18<-0
                    H3HHleftk23<-0
                    H3HHleftk56<-0
                    l<-2
                    for(f in 1:(dim(H3HH)[1]/4))
                    {

                      if(unlist(strsplit(H3HH[(l),2],""))[4]=="K")
                      {
                        H3HHleftk4<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[36]=="K")
                      {
                        H3HHleftk36<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[79]=="K")
                      {
                        H3HHleftk79<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[9]=="K")
                      {
                        H3HHleftk9<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[27]=="K")
                      {
                        H3HHleftk27<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[10]=="S")
                      {
                        H3HHlefts10<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[14]=="K")
                      {
                        H3HHleftk14<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[18]=="K")
                      {
                        H3HHleftk18<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[23]=="K")
                      {
                        H3HHleftk23<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[56]=="K")
                      {
                        H3HHleftk56<-1
                      }
                      l<-l+4
                    }

                    H3SCleftk4<-0
                    H3SCleftk36<-0
                    H3SCleftk79<-0
                    H3SCleftk9<-0
                    H3SCleftk27<-0
                    H3SClefts10<-0
                    H3SCleftk14<-0
                    H3SCleftk18<-0
                    H3SCleftk23<-0
                    H3SCleftk56<-0
                    l<-2
                    for(f in 1:(dim(H3SC)[1]/4))
                    {

                      if(unlist(strsplit(H3SC[(l),2],""))[4]=="K")
                      {
                        H3SCleftk4<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[36]=="K")
                      {
                        H3SCleftk36<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[79]=="K")
                      {
                        H3SCleftk79<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[9]=="K")
                      {
                        H3SCleftk9<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[27]=="K")
                      {
                        H3SCleftk27<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[10]=="S")
                      {
                        H3SClefts10<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[14]=="K")
                      {
                        H3SCleftk14<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[18]=="K")
                      {
                        H3SCleftk18<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[23]=="K")
                      {
                        H3SCleftk23<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[56]=="K")
                      {
                        H3SCleftk56<-1
                      }
                      l<-l+4
                    }
                    ###H4 MARKS ###
                    H4HHleftk5<-0
                    H4HHleftk8<-0
                    H4HHleftk12<-0
                    H4HHleftk14<-0
                    H4HHleftk16<-0
                    H4HHleftk20<-0
                    l<-2
                    for(f in 1:(dim(H4HH)[1]/4))
                    {

                      if(unlist(strsplit(H4HH[(l),2],""))[5]=="K")
                      {
                        H4HHleftk5<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[8]=="K")
                      {
                        H4HHleftk8<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[12]=="K")
                      {
                        H4HHleftk12<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[14]=="K")
                      {
                        H4HHleftk14<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[16]=="K")
                      {
                        H4HHleftk16<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[20]=="K")
                      {
                        H4HHleftk20<-1
                      }
                      l<-l+4
                    }

                    H4SCleftk5<-0
                    H4SCleftk8<-0
                    H4SCleftk12<-0
                    H4SCleftk14<-0
                    H4SCleftk16<-0
                    H4SCleftk20<-0
                    l<-2
                    for(f in 1:(dim(H4SC)[1]/4))
                    {

                      if(unlist(strsplit(H4SC[(l),2],""))[5]=="K")
                      {
                        H4SCleftk5<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[8]=="K")
                      {
                        H4SCleftk8<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[12]=="K")
                      {
                        H4SCleftk12<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[14]=="K")
                      {
                        H4SCleftk14<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[16]=="K")
                      {
                        H4SCleftk16<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[20]=="K")
                      {
                        H4SCleftk20<-1
                      }
                      l<-l+4
                    }
                    ####

                    ### Right adjacent residue ###

                    ###H2bub###
                    l<-2
                    H2BHHrighthits<-0
                    for(f in 1:(dim(H2BHH)[1]/4))
                    {

                      if(unlist(strsplit(H2BHH[(l),2],""))[122]=="K")
                      {
                        H2BHHrighthits<-1
                      }
                      l<-l+4
                    }
                    l<-2
                    H2BSCrighthits<-0
                    for(f in 1:(dim(H2BSC)[1]/4))
                    {

                      if(unlist(strsplit(H2BSC[(l),2],""))[125]=="K")
                      {
                        H2BSCrighthits<-1
                      }
                      l<-l+4
                    }
                    ###H3Marks###
                    H3HHrightk4<-0
                    H3HHrightk36<-0
                    H3HHrightk79<-0
                    H3HHrightk9<-0
                    H3HHrightk27<-0
                    H3HHrights10<-0
                    H3HHrightk14<-0
                    H3HHrightk18<-0
                    H3HHrightk23<-0
                    H3HHrightk56<-0
                    l<-2
                    for(f in 1:(dim(H3HH)[1]/4))
                    {

                      if(unlist(strsplit(H3HH[(l),2],""))[6]=="K")
                      {
                        H3HHrightk4<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[38]=="K")
                      {
                        H3HHrightk36<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[81]=="K")
                      {
                        H3HHrightk79<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[11]=="K")
                      {
                        H3HHrightk9<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[29]=="K")
                      {
                        H3HHrightk27<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[12]=="S")
                      {
                        H3HHrights10<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[16]=="K")
                      {
                        H3HHrightk14<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[20]=="K")
                      {
                        H3HHrightk18<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[25]=="K")
                      {
                        H3HHrightk23<-1
                      }
                      if(unlist(strsplit(H3HH[(l),2],""))[58]=="K")
                      {
                        H3HHrightk56<-1
                      }
                      l<-l+4
                    }

                    H3SCrightk4<-0
                    H3SCrightk36<-0
                    H3SCrightk79<-0
                    H3SCrightk9<-0
                    H3SCrightk27<-0
                    H3SCrights10<-0
                    H3SCrightk14<-0
                    H3SCrightk18<-0
                    H3SCrightk23<-0
                    H3SCrightk56<-0
                    l<-2
                    for(f in 1:(dim(H3SC)[1]/4))
                    {

                      if(unlist(strsplit(H3SC[(l),2],""))[6]=="K")
                      {
                        H3SCrightk4<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[38]=="K")
                      {
                        H3SCrightk36<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[81]=="K")
                      {
                        H3SCrightk79<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[11]=="K")
                      {
                        H3SCrightk9<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[29]=="K")
                      {
                        H3SCrightk27<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[12]=="S")
                      {
                        H3SCrights10<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[16]=="K")
                      {
                        H3SCrightk14<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[20]=="K")
                      {
                        H3SCrightk18<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[25]=="K")
                      {
                        H3SCrightk23<-1
                      }
                      if(unlist(strsplit(H3SC[(l),2],""))[58]=="K")
                      {
                        H3SCrightk56<-1
                      }
                      l<-l+4
                    }
                    ###H4 MARKS ###
                    H4HHrightk5<-0
                    H4HHrightk8<-0
                    H4HHrightk12<-0
                    H4HHrightk14<-0
                    H4HHrightk16<-0
                    H4HHrightk20<-0
                    l<-2
                    for(f in 1:(dim(H4HH)[1]/4))
                    {

                      if(unlist(strsplit(H4HH[(l),2],""))[7]=="K")
                      {
                        H4HHrightk5<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[10]=="K")
                      {
                        H4HHrightk8<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[14]=="K")
                      {
                        H4HHrightk12<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[16]=="K")
                      {
                        H4HHrightk14<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[18]=="K")
                      {
                        H4HHrightk16<-1
                      }
                      if(unlist(strsplit(H4HH[(l),2],""))[22]=="K")
                      {
                        H4HHrightk20<-1
                      }
                      l<-l+4
                    }

                    H4SCrightk5<-0
                    H4SCrightk8<-0
                    H4SCrightk12<-0
                    H4SCrightk14<-0
                    H4SCrightk16<-0
                    H4SCrightk20<-0
                    l<-2
                    for(f in 1:(dim(H4SC)[1]/4))
                    {

                      if(unlist(strsplit(H4SC[(l),2],""))[7]=="K")
                      {
                        H4SCrightk5<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[10]=="K")
                      {
                        H4SCrightk8<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[14]=="K")
                      {
                        H4SCrightk12<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[16]=="K")
                      {
                        H4SCrightk14<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[18]=="K")
                      {
                        H4SCrightk16<-1
                      }
                      if(unlist(strsplit(H4SC[(l),2],""))[22]=="K")
                      {
                        H4SCrightk20<-1
                      }
                      l<-l+4
                    }


                    ###

                    ###END###

                    Functions<-as.data.frame(read.csv(paste0(dir(),"/Interprofunctions.csv"), header = TRUE))
                    functionhits<-rep(NA,(dim(Functions)[1]))
                    for(y in 1:dim(Functions)[1])
                    {
                    functionhits[y]<-IP_table_w_hits[which(IP_table_w_hits[,1]==Functions[y,1]),2][1]
                    }
                    Functions<-cbind(Functions, functionhits)
                    Functions[which(is.na(Functions[,25])==TRUE),25]<-0

                    writer<-rep(0,20)
                    for(ww in 3:22)
                    {
                        if(length(which(Functions$functionhits[which(is.na(Functions[,ww])==FALSE & Functions$X.2==1)]>0))==0)
                        {
                        writer[ww-2]<-0
                        }else if(length(which(Functions$Specific[which(Functions$functionhits[which(is.na(Functions[,ww])==FALSE & Functions$X.2==1)]>0)]>0))>0){
                          writer[ww-2]<-0.5
                        }else{
                          writer[ww-2]<-0.25
                        }
                    }
                    eraser<-rep(0,20)
                    for(ww in 3:22)
                    {
                      if(length(which(Functions$functionhits[which(is.na(Functions[,ww])==FALSE & Functions$X.2==(-1))]>0))==0)
                      {
                        eraser[ww-2]<-0
                      }else if(length(which(Functions$Specific[which(Functions$functionhits[which(is.na(Functions[,ww])==FALSE & Functions$X.2==(-1))]>0)]>0))>0){
                        eraser[ww-2]<-0.5
                      }else{
                        eraser[ww-2]<-0.25
                      }
                    }
                    modifier<-writer+eraser
                    modifier[which(modifier[1:4]==0.5)]<-1
                    modifiernum<-modifier
                    descri<-matrix("NONE",20,8)
                    for(ww in 3:22)
                    {
                    if(modifier[ww-2]>0)
                    {
                    descri[ww-2,5]<-sum(Functions$functionhits[which(is.na(Functions[,ww])==FALSE)])
                    descri[ww-2,8]<-paste(Functions$X.1[is.na(Functions[,ww])==FALSE & Functions$functionhits>0], collapse = " | ")

                    }
                    }
                    writer[which(writer==0)]<-NA
                    writer[which(writer==0.25)]<-"None-specific"
                    writer[which(writer==0.5)]<-"Specific"

                    eraser[which(eraser==0)]<-NA
                    eraser[which(eraser==0.25)]<-"None-specific"
                    eraser[which(eraser==0.5)]<-"Specific"

                    descri[,6]<-writer
                    descri[,7]<-eraser
                  residue<-as.vector(as.numeric(c(NA,NA,NA,NA,sum(H3HHk4,H3SCk4), sum(H3HHk36,H3SCk36),sum(H3HHk9,H3SCk9),sum(H3HHk9,H3SCk9), sum(H3HHk27,H3SCk27),sum(H4HHk5,H4SCk5),
                    sum(H4HHk12,H4SCk12),sum(H3HHk9,H3SCk9),sum(H3HHk27,H3SCk27),sum(H4HHk16,H4SCk16),sum(H4HHk5,H4SCk5), sum(H4HHk20,H4SCk20), sum(H3HHk56,H3SCk56)
                    ,sum(H3HHk14,H3SCk14),sum(H4HHk8,H4SCk8),sum(H3HHk79,H3SCk79))))
                  residuenum<-residue
                  residue[which(residue>0)]<-"Yes"
                  residue[which(residue==0)]<-"No"


                  residueleft<-as.vector(as.numeric(c(NA,NA,NA,NA,sum(H3HHleftk4,H3SCleftk4), sum(H3HHleftk36,H3SCleftk36),sum(H3HHleftk9,H3SCleftk9),sum(H3HHleftk9,H3SCleftk9),sum(H3HHleftk27,H3SCleftk27),sum(H4HHleftk5,H4SCleftk5),
                                                  sum(H4HHleftk12,H4SCleftk12),sum(H3HHleftk9,H3SCleftk9),sum(H3HHleftk27,H3SCleftk27),sum(H4HHleftk16,H4SCleftk16),sum(H4HHleftk5,H4SCleftk5), sum(H4HHleftk20,H4SCleftk20),sum(H3HHleftk56,H3SCleftk56)
                                                  ,sum(H3HHleftk14,H3SCleftk14),sum(H4HHleftk8,H4SCleftk8),sum(H3HHleftk79,H3SCleftk79))))
                  residuenumleft<-residueleft
                  residueleft[which(residueleft>0)]<-"Yes"
                  residueleft[which(residueleft==0)]<-"No"

                  residueright<-as.vector(as.numeric(c(NA,NA,NA,NA, sum(H3HHrightk4,H3SCrightk4), sum(H3HHrightk36,H3SCrightk36),sum(H3HHrightk9,H3SCrightk9),sum(H3HHrightk9,H3SCrightk9),sum(H3HHrightk27,H3SCrightk27),sum(H4HHrightk5,H4SCrightk5),
                                                  sum(H4HHrightk12,H4SCrightk12),sum(H3HHrightk9,H3SCrightk9),sum(H3HHrightk27,H3SCrightk27),sum(H4HHrightk16,H4SCrightk16),sum(H4HHrightk5,H4SCrightk5),sum(H4HHrightk20,H4SCrightk20), sum(H3HHrightk56,H3SCrightk56)
                                                  ,sum(H3HHrightk14,H3SCrightk14),sum(H4HHrightk8,H4SCrightk8),sum(H3HHrightk79,H3SCrightk79))))
                  residuenumright<-residueright
                  residueright[which(residueright>0)]<-"Yes"
                  residueright[which(residueright==0)]<-"No"
                  descri[,2]<-residue
                  descri[,3]<-residueleft
                  descri[,4]<-residueright
                  probability<-rep(0,20)
                  probability[1:4]<-modifiernum[1:4]
                  probability[5:20]<-modifiernum[5:20]*0.5
                  for (ty in 5:20)
                    if(residuenum[ty]>0)
                    {
                    probability[ty]<-probability[ty]+0.5
                    }else if (residuenumleft[ty]>0)
                  {probability[ty]<-probability[ty]+0.25

                  }else if (residuenumright[ty]>0)
                  {
                    probability[ty]<-probability[ty]+0.25
                  }
                  H2bub<-c("Zero (0%)", "No", "No", "No", "Not Measured", "Not Measured", "Not Measured", "Not Measured")
                  if(sum(H2BHHhits, H2BSChits)>0)
                  {
                    H2bub[1]<-"Good (62.5%)"
                  } else if (sum(H2BHHlefthits, H2BSClefthits)>0)
                  {
                    H2bub[1]<-"Poor (37.5%)"
                  }else if (sum(H2BHHrighthits, H2BSCrighthits)>0)
                  {
                    H2bub[1]<-"Poor (37.5%)"
                  }

                  if(sum(H2BHHhits, H2BSChits)>0)
                  {
                    H2bub[2]<-"Yes"
                  }
                  if (sum(H2BHHlefthits, H2BSClefthits)>0)
                  {
                    H2bub[3]<-"Yes"
                  }
                  if (sum(H2BHHrighthits, H2BSCrighthits)>0)
                  {
                    H2bub[4]<-"Yes"
                  }
                  if(H2bub[1]!="Good (62.5%)")
                  {probability[20]<-probability[20]-0.125}

                  probablity<-as.vector(probability)
                  probablity[which(probability==0)]<-"Zero (0%)"
                  probablity[which(probability==0.250)]<-"Negligible (25%)"
                  probablity[which(probability==0.375)]<-"Poor (37.5%)"
                  probablity[which(probability==0.5)]<-"Fair (50%)"
                  probablity[which(probability==0.625)]<-"Good (62.5%)"
                  probablity[which(probability==0.75)]<-"Better (75%)"
                  probablity[which(probability==0.875)]<-"Excellent (87.5%)"
                  probablity[which(probability==1)]<-"Best (100%)"
                  descri[,1]<-probablity
                  descri<-as.data.frame(descri)

                  descri<-rbind(descri, H2bub)
                  descri[which(descri[,7]=="NONE"),7]<-NA
                  colnames(descri)<-c("Probability", "Correct Residue", "Correct Residue left of canonical", "Correct Residue right of canonical", "# of HMM hits", "Specific Writer", "Specific Eraser", "Hit Descriptions")
                  rownames(descri)<-c("4mc",	"5mc",	"6ma",	"5hmc",	"H3K4me",	"H3K36me",	"H3K9me1/2", "H3K9me3",	"H3K27me",	"H4K5ac",	"H4K12ac",	"H3K9ac",	"H3K27ac",	"H4K16ac",	"H4K5me",	"H4K20me",	"H3K56ac",	"H3K14ac",	"H4K8ac",	"H3K79me", "H2Bub")
                  descri<<-descri


                  output$Predictions<-renderReactable({reactable(descri, pagination = FALSE, searchable = FALSE,
                                                            columns = list(
                                                              ".rownames" = colDef(width=90),
                                                              "Probability" = colDef(width=100,
                                                                                     style = function(value) {
                                                                                       if (grepl("Zero (0%)", value))
                                                                                       {color <- "#CA300F"
                                                                                       }else if (grepl("Negligible (25%)", value))
                                                                                       {color <- "#F76D4F"
                                                                                       }else if (grepl("Poor (37.5%)", value))
                                                                                       {color <- "#FFAB3F"
                                                                                       }else if (grepl("Fair (50%)", value))
                                                                                       {color <- "#FFc23F"
                                                                                       }else if (grepl("Good (62.5%)", value))
                                                                                       {color <- "#BDE74D"
                                                                                       }else if (grepl("Better (75%)", value))
                                                                                       {color <- "#81D63D"
                                                                                       }else if (grepl("Excellent (87.5%)", value))
                                                                                       {color <- "#0Fa300"
                                                                                       }else if (grepl("Best (100%)", value))
                                                                                       {color <- "#0c8200"
                                                                                       }else
                                                                                       {color <- "black"}
                                                                                       list(color = color, fontWeight = "bold", borderRight = "1px solid #eee" )
                                                                                     }
                                                                                    ),
                                                              "Correct Residue" = colDef(width=90,
                                                                                         style = function(value1) {
                                                                                           if (grepl("No", value1))
                                                                                           {color1 <- "#CA300F"
                                                                                           }else if (grepl("Yes", value1))
                                                                                           {color1 <- "#387C03"
                                                                                           }else
                                                                                           { color1 <- "black"}
                                                                                           list(color = color1, fontWeight = "bold", borderRight = "1px solid #eee" )
                                                                                         }
                                                                                         ),
                                                              "Correct Residue left of canonical" = colDef(width=90,
                                                                                                           style = function(value1) {
                                                                                                             if (grepl("No", value1))
                                                                                                             {color1 <- "#CA300F"
                                                                                                             }else if (grepl("Yes", value1))
                                                                                                             {color1 <- "#387C03"
                                                                                                             }else
                                                                                                             { color1 <- "black"}
                                                                                                             list(color = color1, fontWeight = "bold", borderRight = "1px solid #eee" )
                                                                                                           }
                                                              ),
                                                              "Correct Residue right of canonical" = colDef(width=90,
                                                                                                            style = function(value1) {
                                                                                                              if (grepl("No", value1))
                                                                                                              {color1 <- "#CA300F"
                                                                                                              }else if (grepl("Yes", value1))
                                                                                                              {color1 <- "#387C03"
                                                                                                              }else
                                                                                                              { color1 <- "black"}
                                                                                                              list(color = color1, fontWeight = "bold", borderRight = "1px solid #eee" )
                                                                                                            }
                                                              ),
                                                              "# of HMM hits" = colDef(width=90,
                                                                                       style = function(value1) {
                                                                                         if (grepl("NONE", value1))
                                                                                         {color1 <- "#E3E3E3"
                                                                                         }else
                                                                                         { color1 <- "black"}
                                                                                         list(color = color1, fontWeight = "bold", borderRight = "1px solid #eee" )
                                                                                       }
                                                              ),
                                                              "Specific Writer" = colDef(width=90,
                                                                                                               style = function(value1) {
                                                                                                                 if (grepl("Specific", value1))
                                                                                                                 {color1 <- "black"
                                                                                                                 }else
                                                                                                                 { color1 <- "#E3E3E3"}
                                                                                                                 list(color = color1, fontWeight = "bold", borderRight = "1px solid #eee" )
                                                                                                               }
                                                              ),
                                                              "Specific Eraser" = colDef(width=90,
                                                                                         style = function(value1) {
                                                                                           if (grepl("Specific", value1))
                                                                                           {color1 <- "black"
                                                                                           }else
                                                                                           { color1 <- "#E3E3E3"}
                                                                                           list(color = color1, fontWeight = "bold", borderRight = "1px solid #eee" )
                                                                                         }
                                                              ),
                                                              "Hit Descriptions" = colDef()



                                                            )



                  )})








                    outtable<-similarity(histone = his, species = sp)
                    outtable<<-as.data.frame(outtable)
                    tnames<-colnames(outtable)
                    columns <- setNames(
                      list(
                        colDef(width=300,
                               sticky = "left",
                               html=TRUE,
                               align="right",
                               # Add a right border style to visually distinguish the sticky column

                               headerStyle = list(borderRight = "1px solid #eee",fontWeight = "bold", fontSize=20 ,fontFamily= "monospace"),

                               style = function(value) {
                                if (grepl("Deletion", value)) {
                                   color <- "firebrick"
                                 } else if (grepl("Substitution", value)) {
                                   color <- "maroon"
                                 } else if (grepl("Match", value)) {
                                   color <- "forestgreen"
                                 } else {
                                   color <- "black"
                                 }
                                 list(color = color, fontWeight = "bold", borderRight = "1px solid #eee" )
                               }
                        ),
                        colDef(headerStyle = list(fontWeight = "bold", fontSize=20 ,fontFamily= "monospace"), minWidth=3000, html=TRUE, style = function(value) {
                        list(fontWeight = "bold", fontSize=20 ,fontFamily= "monospace") })
                      ),
                      tnames
                    )

                    position<-which(grepl("Protein", outtable[,1])==TRUE)
                    output$outtable<-renderReactable({reactable(outtable, pagination = FALSE, sortable=FALSE, rownames = FALSE, height= 375,
                                                                rowStyle = function(index) {
                                                                  if (index %in% position) {
                                                                    list(`border-top` = "thin solid")
                                                                  }
                                                                },
                                                                columns= columns



                                                                )})




                    output$message <- renderText("Data has been analyzed. Evidence and Histone similarity/antibodies tabs are active")


                    }else {output$dir <- renderText("You did not select a directory!")}

                }

               )
observeEvent(input$example,
             {
             showModal(modalDialog(
               title = "Example",
               HTML('<img src="www/schema.jpg" />'),
               easyClose = TRUE,
               footer = NULL
             ))}
             )

observeEvent(input$show_details,

             {
              vacant<-is.na(IP_table_w_hits_details[,8])
              if(vacant[as.numeric(input$show_details)]==FALSE){
              whichlevel<<-1
              output$description <- renderText(paste0("Hits for IPR# ", IP_table_w_hits[as.numeric(input$show_details),1], " (", IP_table_w_hits[as.numeric(input$show_details),3], ")"))
              specificdetails<-Interpro[which(Interpro[,12]==IP_table_w_hits[as.numeric(input$show_details),1]),]
              specificdetails<-specificdetails[,-2]
              specificdetails<-specificdetails[,-10]
              specificdetails<-specificdetails[,1:8]
              if(length(which (specificdetails[,8]=="-"))==0)
              {
              specificdetails<<-specificdetails<-specificdetails[order(decreasing=TRUE, as.numeric(specificdetails[,8])),]
              }
              uni<-unique(specificdetails[,1])


              who<-which(specificdetails[,1]==uni[1])
              for(k in 1:length(who))
              {

                if(k==1)
                {
                  plf<-paste("<br", specificdetails[who[k],2])
                  hmm1<-paste("<br>", specificdetails[who[k],3], "<br>")
                  hnum1<-paste("<br>", specificdetails[who[k],4], "<br>")
                  hd1<-paste("<br>", specificdetails[who[k],5], "<br>")
                  ms1<-paste("<br>", specificdetails[who[k],6], "<br>")
                  mst1<-paste("<br>", specificdetails[who[k],7], "<br>")
                  ev1<-paste("<br>", specificdetails[who[k],8], "<br>")
                  det1<-paste("<br>", specificdetails[who[k],9], "<br>")
                }else
                {
                  hmm2<-paste("<br>", specificdetails[who[k],3], "<br>")
                  hnum2<-paste("<br>", specificdetails[who[k],4], "<br>")
                  hd2<-paste("<br>", specificdetails[who[k],5], "<br>")
                  ev2<-paste("<br>", specificdetails[who[k],8], "<br>")
                  ms2<-paste("<br>", specificdetails[who[k],6], "<br>")
                  mst2<-paste("<br>", specificdetails[who[k],7], "<br>")



                  hmm1<-paste(hmm1,hmm2)
                  hnum1<-paste(hnum1,hnum2)
                  hd1<-paste(hd1,hd2)
                  ms1<-paste(ms1,ms2)
                  mst1<-paste(mst1, mst2)
                  ev1<-paste(ev1,ev2)

                }
              }

              hmmf<-hmm1
              hnumf<-hnum1
              hdf<-hd1
              msf<-ms1
              mstf<-mst1
              evf<-ev1

              if(length(uni)>1)
              {
              for (p in 2:length(uni))
              {

                who<-which(specificdetails==uni[p])
                for(k in 1:length(who))
                {

                  if(k==1)
                  {
                    pl<-paste("<br", specificdetails[who[k],2])
                    hmm1<-paste("<br>", specificdetails[who[k],3], "<br>")
                    hnum1<-paste("<br>", specificdetails[who[k],4], "<br>")
                    hd1<-paste("<br>", specificdetails[who[k],5], "<br>")
                    ms1<-paste("<br>", specificdetails[who[k],6], "<br>")
                    mst1<-paste("<br>", specificdetails[who[k],7], "<br>")
                    ev1<-paste("<br>", specificdetails[who[k],8], "<br>")

                  }else
                  {
                    hmm2<-paste("<br>", specificdetails[who[k],3], "<br>")
                    hnum2<-paste("<br>", specificdetails[who[k],4], "<br>")
                    hd2<-paste("<br>", specificdetails[who[k],5], "<br>")
                    ev2<-paste("<br>", specificdetails[who[k],8], "<br>")
                    ms2<-paste("<br>", specificdetails[who[k],6], "<br>")
                    mst2<-paste("<br>", specificdetails[who[k],7], "<br>")



                    hmm1<-paste(hmm1,hmm2)
                    hnum1<-paste(hnum1,hnum2)
                    hd1<-paste(hd1,hd2)
                    ms1<-paste(ms1,ms2)
                    mst1<-paste(mst1,mst2)
                    ev1<-paste(ev1,ev2)

                  }

                }
                plf<-c(plf,pl)
                hmmf<-c(hmmf,hmm1)
                hnumf<-c(hnumf,hnum1)
                hdf<-c(hdf,hd1)
                msf<-c(msf,ms1)
                mstf<-c(mstf,mst1)
                evf<-c(evf,ev1)







              }
              }
              specificdetails<-cbind(uni, plf, hmmf,hnumf,hdf,msf,mstf,evf)
              colnames(specificdetails)[1]<-c("Transcript number")
              colnames(specificdetails)[2]<-c("Protein length")
              colnames(specificdetails)[3]<-c("HMM Analysis")
              colnames(specificdetails)[4]<-c("HMM accession #")
              colnames(specificdetails)[5]<-c("HMM description")
              colnames(specificdetails)[6]<-c("Match start")
              colnames(specificdetails)[7]<-c("Match stop")
              colnames(specificdetails)[8]<-c("Score (e-value)")

              specificdetails<-cbind(specificdetails, details = NA)

              transcipthits<-unique(blastp[,1])

              for (y in 1:dim(specificdetails)[1])
              {
                if (any(transcipthits %in% specificdetails[y,1]))
                {
                specificdetails[y,9]<- as.character(tags$div(tags$button("Blast+ results for transcript")))}
              }

              specificdetails<<-specificdetails
              output$IPhits<-renderReactable({reactable(specificdetails, pagination = FALSE, searchable = TRUE, columns = list(

                "Transcript number" = colDef(
                  width = 125,
                  style = function(value) {
                    list(fontWeight = "bold")
                  }
                ),
                "Protein length" = colDef(
                  html=TRUE,
                  width = 125
                ),
                "HMM Analysis" = colDef(
                  html=TRUE,
                  width = 125
                ),
                "HMM accession #" = colDef(
                  html=TRUE,
                  width = 125
                ),
                "HMM description" = colDef(
                  html=TRUE,
                  minWidth=200
                ),
                "Match start" = colDef(
                  html=TRUE,
                  width = 125
                ),
                "Match stop" = colDef(
                  html=TRUE,
                  width = 125
                ),
                "Score (e-value)" = colDef(
                  html=TRUE,
                  width = 125
                ),
                  details = colDef(
                  name = "",
                  sortable = TRUE,
                  html=TRUE
                )),onClick = JS("function(rowInfo, column) {
    if (column.id !== 'details') {
      return
    }

    if (window.Shiny) {
      Shiny.setInputValue('show_details2', { index: rowInfo.index + 1 }, { priority: 'event' })
    }

  }"))})



            }

             }
             )

observeEvent(input$show_details2,

             {
               vacant<-is.na(specificdetails[,9])
               if(vacant[as.numeric(input$show_details2)]==FALSE){
               whichlevel<<-2
               output$description <- renderText(paste0("Blast results for transcript ID: ", specificdetails[as.numeric(input$show_details2),1]))
               blastp[,2]<-protname
               somestuff<-blastp[which(blastp[,1]==specificdetails[as.numeric(input$show_details2),1]),]
               somestuff<-somestuff[,-1]
               somestuff<-somestuff[order(as.numeric(somestuff[,11]), decreasing= TRUE), ]
               somestuff<-somestuff[,-11]


                   blast<-somestuff[,1]
                   listofstrings<-c()
                   for(j in 1:length(blast))
                   {
                     minimal<-unlist(strsplit(blast[j], " "))[c(-1)]
                     minimal<-minimal[-c(which(grepl("\\[", minimal)==TRUE):which(grepl("\\]", minimal)==TRUE))]
                     listofstrings[j]<-paste(minimal, collapse=" ")
                     rm(minimal)
                   }
                   freq<-as.data.frame(table(listofstrings))
                   colnames(freq)<-c("words","freq")
                   freq<-freq[order(freq[,2], decreasing = TRUE),]
                   if(dim(freq)[1]>10)
                   {
                     freq<-freq[1:10,]
                   }
                   showModal(modalDialog(tags$head(tags$style(HTML(".modal-dialog { width: 90vw; height: 90vw }"))),
                                         renderPlot(
                                           ggplot(freq, aes(x=words, y=freq)) +
                                             geom_segment( aes(x=words, xend=words, y=0, yend=freq), color="skyblue") +
                                             geom_point( color="darkblue", size=4, alpha=0.6) +
                                             theme_light() +
                                             coord_flip() +
                                             ylab("# of BLAST+ Hits") +
                                             xlab("Homologous Protein Identifier") +
                                             ggtitle(paste0("Top Homologous Protein Identifiers for transcript: ", specificdetails[as.numeric(input$show_details2),1]))+
                                             theme(
                                               panel.grid.major.y = element_blank(),
                                               panel.border = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text = element_text(size=14, face= "bold", colour="black"),
                                               axis.title.y  = element_text(size=14, face= "bold", colour="darkred", margin = margin(t = 0, r = 15, b = 0, l = 0)),
                                               axis.title.x  = element_text(size=14, face= "bold", colour="darkred"),
                                               plot.title = element_text(color="black", size=16, face="bold", hjust = 0.5),

                                             )







                                         )))



               colnames(somestuff)[1]<-c("Model Organism Sequence ID and Description")
               colnames(somestuff)[2]<-c("Percent Match")
               colnames(somestuff)[3]<-c("Alignment Length")
               colnames(somestuff)[4]<-c("# Mismatches")
               colnames(somestuff)[5]<-c("# gap openings")
               colnames(somestuff)[6]<-c("Alignment Start Transcript")
               colnames(somestuff)[7]<-c("Alignment Stop Transcript")
               colnames(somestuff)[8]<-c("Alignment Start Model Organism")
               colnames(somestuff)[9]<-c("Alignment Stop Model Organism")
               colnames(somestuff)[10]<-c("Score (e-value)")

               output$IPhits<-renderReactable({reactable(somestuff, pagination = FALSE, searchable = TRUE, columns = list(

                 "Model Organism Sequence ID and Description" = colDef(minWidth = 750,
                   style = function(value) {
                     list(fontWeight = "bold")
                   }
                 ),
                 "Percent Match"= colDef(width=150),
                 "Alignment Length" = colDef(width=150),
                 "# Mismatches" = colDef(width=150),
                 "# gap openings" = colDef(width=150),
                 "Alignment Start Transcript" = colDef(width=150),
                 "Alignment Stop Transcript" = colDef(width=150),
                 "Alignment Start Model Organism" = colDef(width=150),
                 "Alignment Stop Model Organism" = colDef(width=150),
                 "Score (e-value)" = colDef(width=150)
                 ))
               })



              } })







observeEvent(input$back,
             if (whichlevel==1)
             {
            output$description <- renderText("InterProScan HMM Hits")
             suppressWarnings(rm(whichlevel))
               description<-paste0("InterProScan HMM Hits")
               output$IPhits<-renderReactable({reactable(IP_table_w_hits_details, pagination = FALSE, searchable = TRUE,
                                                         columns = list(
                                                           Accession = colDef(width=90),
                                                           "HMM Protein Homology Hits" = colDef(width=120, na="0",
                                                                                                cell = icon_sets(icon_ref = "icontable", icon_color_ref ="colortable", IP_table_w_hits_details),
                                                                                                style = function(value) {
                                                                                                  if (is.na(value)) {
                                                                                                    color <- "lightgrey"
                                                                                                  } else
                                                                                                  {
                                                                                                    color <- "darkblue"
                                                                                                  }
                                                                                                  list(color = color, fontWeight = "bold")
                                                                                                }
                                                           ),
                                                           NAME = colDef(minWidth=200),
                                                           'HMM Database' = colDef(width=100),
                                                           Class = colDef(minWidth=135),
                                                           'Integrated Signatures' = colDef(minWidth=125),
                                                           'Go Terms' = colDef(minWidth=175),
                                                           icontable = colDef(show=FALSE),
                                                           colortable = colDef(show=FALSE),
                                                           details = colDef(
                                                             name = "",
                                                             sortable = TRUE,
                                                             html = TRUE)

                                                         ),
                                                         onClick = JS("function(rowInfo, column) {
    // Only handle click events on the 'details' column
    if (column.id !== 'details') {
      return
    }

    // Send the click event to Shiny, which will be available in input$show_details
    // Note that the row index starts at 0 in JavaScript, so we add 1
    if (window.Shiny) {
      Shiny.setInputValue('show_details', { index: rowInfo.index + 1 }, { priority: 'event' })
    }
  }")



               )})


             }
             else if (whichlevel==2)
             {
               whichlevel<<-1
               output$description <- renderText(paste0("Hits for IPR# ", IP_table_w_hits[as.numeric(input$show_details),1], " (", IP_table_w_hits[as.numeric(input$show_details),3], ")"))
               output$IPhits<-renderReactable({reactable(specificdetails, pagination = FALSE, searchable = TRUE, columns = list(

                 "Transcript number" = colDef(
                   width = 125,
                   style = function(value) {
                     list(fontWeight = "bold")
                   }
                 ),
                 "Protein length" = colDef(
                   html=TRUE,
                   width = 125
                 ),
                 "HMM Analysis" = colDef(
                   html=TRUE,
                   width = 125
                 ),
                 "HMM accession #" = colDef(
                   html=TRUE,
                   width = 125
                 ),
                 "HMM description" = colDef(
                   html=TRUE,
                   minWidth=200
                 ),
                 "Match start" = colDef(
                   html=TRUE,
                   width = 125
                 ),
                 "Match stop" = colDef(
                   html=TRUE,
                   width = 125
                 ),
                 "Score (e-value)" = colDef(
                   html=TRUE,
                   width = 125
                 ),
                 details = colDef(
                   name = "",
                   sortable = TRUE,
                   html=TRUE
                 )),onClick = JS("function(rowInfo, column) {
    // Only handle click events on the 'details' column
    if (column.id !== 'details') {
      return
    }

    // Send the click event to Shiny, which will be available in input$show_details2
    // Note that the row index starts at 0 in JavaScript, so we add 1
    if (window.Shiny) {
      Shiny.setInputValue('show_details2', { index: rowInfo.index + 1 }, { priority: 'event' })
    }
  }"))})



               }

             )
observeEvent(input$species,
             {
               sp<<-input$species
               outtable<-similarity(histone = his, species = sp)
               outtable<<-as.data.frame(outtable)



               tnames<-colnames(outtable)
               columns <- setNames(
                 list(
                   colDef(width=300,
                          sticky = "left",
                          html=TRUE,
                          align="right",
                          # Add a right border style to visually distinguish the sticky column

                          headerStyle = list(borderRight = "1px solid #eee",fontWeight = "bold", fontSize=20 ,fontFamily= "monospace"),

                          style = function(value) {
                            if (grepl("Deletion", value)) {
                              color <- "firebrick"
                            } else if (grepl("Substitution", value)) {
                              color <- "maroon"
                            } else if (grepl("Match", value)) {
                              color <- "forestgreen"
                            } else {
                              color <- "black"
                            }
                            list(color = color, fontWeight = "bold", borderRight = "1px solid #eee" )
                          }
                   ),
                   colDef(headerStyle = list(fontWeight = "bold", fontSize=20 ,fontFamily= "monospace"), minWidth=3000, html=TRUE, style = function(value) {
                     list(fontWeight = "bold", fontSize=20 ,fontFamily= "monospace") })
                 ),
                 tnames
               )

               position<-which(grepl("Protein", outtable[,1])==TRUE)
               output$outtable<-renderReactable({reactable(outtable, pagination = FALSE, sortable=FALSE, rownames = FALSE, height= 375,
                                                           rowStyle = function(index) {
                                                             if (index %in% position) {
                                                               list(`border-top` = "thin solid")
                                                             }
                                                           },
                                                           columns= columns



               )})
},ignoreInit = TRUE)
observeEvent(input$histone,
             {
               his<<-input$histone
               outtable<-similarity(histone = his, species = sp)
               outtable<<-as.data.frame(outtable)
               tnames<-colnames(outtable)
               columns <- setNames(
                 list(
                   colDef(width=300,
                          sticky = "left",
                          html=TRUE,
                          align="right",
                          # Add a right border style to visually distinguish the sticky column

                          headerStyle = list(borderRight = "1px solid #eee",fontWeight = "bold", fontSize=20 ,fontFamily= "monospace"),

                          style = function(value) {
                            if (grepl("Deletion", value)) {
                              color <- "firebrick"
                            } else if (grepl("Substitution", value)) {
                              color <- "maroon"
                            } else if (grepl("Match", value)) {
                              color <- "forestgreen"
                            } else {
                              color <- "black"
                            }
                            list(color = color, fontWeight = "bold", borderRight = "1px solid #eee" )
                          }
                   ),
                   colDef(headerStyle = list(fontWeight = "bold", fontSize=20 ,fontFamily= "monospace"), minWidth=3000, html=TRUE, style = function(value) {
                     list(fontWeight = "bold", fontSize=20 ,fontFamily= "monospace") })
                 ),
                 tnames
               )

               position<-which(grepl("Protein", outtable[,1])==TRUE)
               output$outtable<-renderReactable({reactable(outtable, pagination = FALSE, sortable=FALSE, rownames = FALSE, height= 375,
                                                           rowStyle = function(index) {
                                                             if (index %in% position) {
                                                               list(`border-top` = "thin solid")
                                                             }
                                                           },
                                                           columns= columns



               )})
},ignoreInit = TRUE)




}


ui <- dashboardPage(scrollToTop = TRUE, skin="green",
  #freshTheme = theme,
  header = dashboardHeader(title="PERCEPTIVE", controlbarIcon =icon("moon")),
  sidebar = dashboardSidebar(sidebarMenu(
    menuItem("Select Data", tabName = "DataInput", icon = icon("database")),
    menuItem("Predictions", tabName = "Predictions", icon = icon("lightbulb")),
    menuItem("Evidence", tabName = "Analysis", icon = icon("magnifying-glass-chart")),
    menuItem("Histone Similarity/Antibody", tabName = "Similarity", icon = icon("microscope"))

  )),
  body = dashboardBody(
    #use_theme(mytheme),
    tabItems(
      # First tab content
      tabItem(tabName = "DataInput",
              h5(HTML('<b>Select pathway to directory generated by PERCEPTIVE for your organism: </b>')),
              shinyDirButton('folder', 'Select Directory' , 'Select pathway to directory generated by PERCEPTIVE for your organism:', multiple = FALSE,
                             buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder")),
              h5(HTML("<b>Folder Selected: </b>")),
              verbatimTextOutput("dir", placeholder = TRUE),
              div(style = "display:inline-block; float:right",  actionButton('startanalysis', 'Start Analysis', icon("play"),
                                                             style="color: #fff; background-color: #28A745; border-color: #2e6da4")),
              verbatimTextOutput("message")
              ),

      # Second tab content

        tabItem(tabName = "Analysis",
              div( h5(HTML("<b>Analyzing data at path: </b>"), textOutput("source"))),
              h2(textOutput("description")),
              div(h5(HTML("<b>Legend</b>:"), div(icon("circle", style = "color:black;"), ":Histone", icon("pencil", style = "color:firebrick;"), ": Histone Acetyltransferase", icon("pencil", style = "color:darkorange;"), ":Histone Methyltransferase", icon("pencil", style = "color:lightblue;"), ":Other Histone Modifier"))),
              div(h5(div(icon("eraser", style = "color:darkgreen;"), ":Histone Demethylase", icon("eraser", style = "color:pink;"), ": Histone Deacetylase", icon("hammer", style = "color:blue;"), ":Chromatin Remodeler", icon("building", style = "color:purple;"), ":Chromatin Structural Maintenance", icon("dna", style = "color:black;"), ":DNA Methylation"))),
              div(style = "text-align: right",  actionButton('back', 'Back', icon("backward"),
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),


              div(actionButton("csvbutton", label = "Download current selection as CSV", icon("download"), style= "color: #fff; background-color: #337ab7; border-color: #2e6da4",
                               onclick = sprintf("Reactable.downloadDataCSV('%s', '%s')", "IPhits", filename = "data.csv"))),

              reactableOutput("IPhits")


      ),
      tabItem(tabName = "Similarity",


              h2("Histone Similarity/Antibody Compatibility"),

              div (radioButtons("histone", "Select a Histone",
                                c("H1" = "H1",
                                  "H2A" = "H2A",
                                  "H2B" = "H2B",
                                  "H3" = "H3",
                                  "H4" = "H4")), style="display:inline-block; vertical-align:top; "),
               div( radioButtons("species", "Select a Model Organism:",
                           c("Human" = "hh",
                             "Budding Yeast" = "sc")),style="display:inline-block; vertical-align:top;"),









              reactableOutput("outtable"),
              div(
                div(actionButton("csvbutton", label = "Download current selection as CSV", icon("download"), style= "color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                 onclick = sprintf("Reactable.downloadDataCSV('%s', '%s')", "IPhits", filename = "data.csv"))), style="text-align: right"),
              h6(HTML("<b>User Note: Most histone and histone PTM antibodies are peptide antibodies mounted against the N-terminal and C-terminal tails of human or budding yeast histones. If sequence homology is low, the sequence of the peptide antigen should be determined and confirmed in the novel organism protein sequence before purchasing/utilizing any commercial antibodies. The homology predicted above reflects the minimum Levenshtein edit distance for a novel organism peptide sequence. This may not predict PTM potential and adjacent modifiable residues could be modified (e.g., A novel organism has H3K35, while the canonical human/budding yeast residue is H3K36 and is me1/2/3 modified).</b> "))

    ),
    tabItem(tabName = "Predictions",


            h3("Predicted DNA and Histone Epigenetic Modifications"),
            div(h5(HTML("<b>Probability associated with modification is represented on a scale from Zero to Excellent, with corresponding percentages of criterion met as outlined here in this example for H3K4me3: </b>"))),
            div(style = "display:inline-block; float:right",  actionButton('example', 'Example', icon("lightbulb"),
                                                                           style="color: #fff; background-color: #28A745; border-color: #2e6da4")),
            reactableOutput("Predictions"),
    )
  )
),

 footer = dashboardFooter(left = "COPYRIGHT HOLDER: © 2024. Triad National Security, LLC. All rights reserved.", right = "LA-UR-24-21779"),
controlbar = dashboardControlbar(collapsed = TRUE, overlay= FALSE, skin="dark", skinSelector())


)




shinyApp(ui, server)
}
