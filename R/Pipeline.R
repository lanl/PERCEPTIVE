#' @export

Pipeline<-function()
{

  library(purrr)
  library(shiny)
  library(shinyFiles)
  library(DT)
  library(htmltools)
  library(fontawesome)
  library(fresh)
  library(shinydashboard)
  library(shinydashboardPlus)
  library(parallel)
  library(shinyscreenshot)






  server <- function(input, output, session){
    if (is.na(unlist(strsplit(system("whereis singularity", intern = TRUE), split=":"))[2]))
    {
      showModal(modalDialog(
        title = "WARNING",
        paste0("Singularity is a required dependency. Please install singularity. See installation guide, or visit singularity (apptainer) installation guide"),
        footer = NULL
      ))
    }
    setwd(find.package("PERCEPTIVEv3"))
    if(file.exists(paste0(find.package("PERCEPTIVEv3"),"/tmp")))
    {

    }else
    {dir.create(path = paste0(find.package("PERCEPTIVEv3"),"/tmp"))
      setwd(paste0(find.package("PERCEPTIVEv3"),"/tmp"))
      write.csv(as.data.frame(matrix(,5,1)), "Defaults.csv", row.names = FALSE)
    }


    if (Sys.info()["sysname"]!="Linux")
    {

    output$OS<-renderText(paste("This is a ", Sys.info()["sysname"], "operating system, which is incompatible with the PERCEPTIVE pipeline (for OS X users, Darwin is OS X). Please use a linux environment. If R has incorrectly identified your environment, please proceed below."))
    }else
    {output$OS<-renderText(Sys.info()["sysname"])}
    denovo<-"No"
    observeEvent(input$denovo, {
      if (input$denovo=="Yes")
      {output$needapath<- renderText({"<b>You are assembling a <i> de novo </i> genome. Please provide fastq/unaligned BAM file containing unassembled reads associated with your species of interest. Path to fastq/BAM:</b>"})
      }
      else if (input$denovo=="No")
      {output$needapath<- renderText({"<b>You are not assembling a <i> de novo </i> genome. Please provide fasta file of all chromosomes/scaffolds associated with your species of interest. Path to fasta:</b>"})
      }
    })
    volumes = getVolumes() ()
    observe({
      shinyFileChoose(input, 'fasta', roots=volumes)

    }

    )
    observe({
      shinyFileChoose(input, 'fastq', roots=volumes)

    }
    )
    observe({
      shinyFileChoose(input, 'fastq2', roots=volumes)

    }
    )
    observe({
      shinyFileChoose(input, 'fastqR1RNA', roots=volumes)

    }
    )

    observe({
      shinyFileChoose(input, 'fastqR2RNA', roots=volumes)

    }
    )

    observe({
      shinyDirChoose(input, 'endlocation', roots=volumes)

    }
    )

    observe({
      shinyFileChoose(input, 'orthodb', roots=volumes)

    }
    )

    observe({
      shinyFileChoose(input, 'braker', roots=volumes)

    }

    )

    observe({
      shinyFileChoose(input, 'perceptive', roots=volumes)

    }

    )

    observe({
      shinyFileChoose(input, 'genemark', roots=volumes)

    }

    )
    observe({
      shinyFileChoose(input, 'blastdb', roots=volumes)

    }
    )

    setwd(paste0(find.package("PERCEPTIVEv3"),"/tmp"))
defaults<<-read.csv("Defaults.csv", header = TRUE)
if(is.na(defaults[1,1]))
{
  updateTabItems(session, "inTabset", selected = "settings")
}


observeEvent(input$update,{
if(is.na(defaults[1,1]))
{
  brakerlocation<<- as.character(parseFilePaths(roots = volumes, input$braker)[4])
  perceptivelocation<<- as.character(parseFilePaths(roots = volumes, input$perceptive)[4])
  orthodblocation<<- as.character(parseFilePaths(roots = volumes, input$orthodb)[4])
  genemarklocation<<- as.character(parseFilePaths(roots = volumes, input$genemark)[4])
      if (brakerlocation=="character(0)" || perceptivelocation=="character(0)" || orthodblocation=="character(0)" || genemarklocation=="character(0)")
      {
        showModal(modalDialog(
          title = "WARNING",
          paste0("All selections not completed"),
          easyClose = TRUE,
          footer = NULL
        ))
      } else
  {
  system(paste("mv" , perceptivelocation, paste0(find.package("PERCEPTIVEv3"),"/tmp")))
  system(paste("mv" , brakerlocation, paste0(find.package("PERCEPTIVEv3"),"/tmp")))
  system(paste("mv" , orthodblocation, paste0(find.package("PERCEPTIVEv3"),"/tmp")))
  system(paste("mv" , genemarklocation, paste0(find.package("PERCEPTIVEv3"),"/tmp")))
  defaults[1,1]<<-1
  defaults[2:5,1]<<-paste0(find.package("PERCEPTIVEv3"),"/tmp")
  write.csv(defaults, "Defaults.csv", row.names = FALSE)
  updateTabItems(session, "inTabset", selected = "DataInput")

  }
}else
{
    brakerlocation<<- as.character(parseFilePaths(roots = volumes, input$braker)[4])
    perceptivelocation<<- as.character(parseFilePaths(roots = volumes, input$perceptive)[4])
    orthodblocation<<- as.character(parseFilePaths(roots = volumes, input$orthodb)[4])
    genemarklocation<<- as.character(parseFilePaths(roots = volumes, input$genemark)[4])
      if(brakerlocation!="character(0)")
      {
      system(paste("mv" , brakerlocation, paste0(find.package("PERCEPTIVEv3"),"/tmp")))
      }
      if(perceptivelocation!="character(0)")
      {
      system(paste("mv" , perceptivelocation, paste0(find.package("PERCEPTIVEv3"),"/tmp")))
      }
      if(brakerlocation!="character(0)")
      {
      system(paste("mv" , orthodblocation, paste0(find.package("PERCEPTIVEv3"),"/tmp")))
      }
      if(brakerlocation!="character(0)")
      {
      system(paste("mv" , genemarklocation, paste0(find.package("PERCEPTIVEv3"),"/tmp")))
      }
      defaults[2:5,1]<<-paste0(find.package("PERCEPTIVEv3"),"/tmp")
      write.csv(defaults, "Defaults.csv", row.names = FALSE)
}

})

    observeEvent(input$run,{
        speciesout<<-input$speciesname
        finallocation<<- as.character(parseDirPath(roots = volumes, input$endlocation))
        cladeout<<-input$clade
        nprocout<<-input$numcore
        evalout<<-input$eval
        novode<<-input$denovo
        RNAeval<<-input$RNA
        if(novode=="Yes" && RNAeval=="No")
        {
        datalocationout<<- as.character(parseFilePaths(roots = volumes, input$fastq)[4])
        short<<-input$length
            if(short=="Long")
            {
            technology<<-input$tech
            gsize<<-input$genomesize
            if (speciesout=="" || datalocationout=="" || cladeout=="" || nprocout=="" || evalout=="" || length(finallocation)==0 || technology=="" || gsize=="")
            {
              showModal(modalDialog(
                title = "WARNING",
                paste0("All selections not completed"),
                easyClose = TRUE,
                footer = NULL
              ))

            }else{
              if(technology=="Nanopore")
              {tnolog<<-"-nanopore-raw"
              }else if (technology=="Pacbio")
              {
              tnolog<<-"-pacbio"
              }else if (technology=="Pacbio HiFi")
              {tnolog<<-"-pacbio-hifi"}
              canuarguments<-input$passthrough
              if(canuarguments=="No")
              {
              fileConn<-file("Canu.sh")
              writeLines(
                c("mkdir tempforpipeline",
                  "cd tempforpipeline/",
                  paste("cp",noquote(datalocationout),  "input.fastq"),
                  paste0("singularity exec ../Perceptivev0.1.sif canu -d assembly -p assembled genomesize=",gsize, " ", tnolog, " input.fastq"),
                  "cp assembly/assembled.contig.fasta input.fasta"

                ), fileConn)
              close(fileConn)
              }else if (canuarguments=="Yes")
              {
                fileConn<-file("Canu.sh")
                writeLines(
                  c("mkdir tempforpipeline",
                    "cd tempforpipeline/",
                    paste("cp",noquote(datalocationout),  "input.fastq"),
                    paste0("singularity exec ../Perceptivev0.1.sif canu -d assembly -p assembled genomesize=",gsize, " ", tnolog, " input.fastq", canuargs),
                    "cp assembly/assembled.contig.fasta input.fasta"

                  ), fileConn)
                close(fileConn)
              }

              fileConn<-file("RepeatMasker.sh")
              writeLines(
                c("cd tempforpipeline/",
                  paste("singularity overlay create -s 10000 temp.img"),
                  paste("singularity exec --overlay temp.img ../Perceptivev0.1.sif RepeatMasker input.fasta --xsmall -species",cladeout, "-pa",nprocout)


                ), fileConn)
              close(fileConn)

              fileConn<-file("braker.sh")
              writeLines(

                c("cd tempforpipeline/",
                  paste("singularity exec ../braker3.sif braker.pl --genome=input.fasta.masked --threads", nprocout, "--prot_seq=../Eukaryota.fa"))
                , fileConn)
              close(fileConn)


              fileConn<-file("interpro.sh")
              writeLines(
                c(
                  "cd tempforpipeline/",
                  "mkdir interpro",
                  "cp braker/braker.aa interpro/",
                  "cd interpro/",
                  paste("sed", '"s/\\*//g"', "< braker.aa > braker.peptide"),
                  "rm braker.aa",
                  paste("singularity exec --overlay ../temp.img ../../Perceptivev0.1.sif interproscan.sh -i braker.peptide -cpu", nprocout, "-dp --goterms -appl AntiFam-7.0,CDD-3.20,Coils-2.2.1,FunFam-4.3.0,Gene3D-4.3.0,Hamap-2023_05,NCBIfam-15.0,PANTHER-19.0,Pfam-37.0,PIRSF-3.10,PIRSR-2023_05,PRINTS-42.0,ProSitePatterns-2023_05,ProSiteProfiles-2023_05,SFLD-4,SMART-9.0,SUPERFAMILY-1.75"),
                  "rm -r temp/",
                  "rm ../temp.img"

                ), fileConn)
              close(fileConn)

              blastexp<-input$uniqueblast
              if(blastexp=="No")
              {
                fileConn<-file("blasting.sh")
                writeLines(
                  c("cd tempforpipeline/interpro",
                    paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
                    "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
                    "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
                    "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
                    "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
                    paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db /dependencies/modelorgsprot/modelorgsprot -db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
                    paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
                  ), fileConn)
                close(fileConn)
              }else if (blastexp=="Yes")
              {
                blastnewdb<<- as.character(parseFilePaths(roots = volumes, input$blastdb)[4])
                fileConn<-file("blasting.sh")
                writeLines(
                  c("cd tempforpipeline/interpro",
                    paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
                    "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
                    "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
                    "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
                    "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
                    paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db", blastnewdb, "-db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
                    paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
                  ), fileConn)
                close(fileConn)
              }

              fileConn<-file("parsing.sh")
              writeLines(
                c("cd tempforpipeline/interpro",
                  "cut -f2 blastbpstricteval.table >blastpstrictevalgenes.txt",
                  "for line in $(cat blastpstrictevalgenes.txt); do grep -w $line ../../../temp.faa >> Speciesandprotein.txt; done",
                  paste("grep", '"IPR002119"', "braker.peptide.tsv | cut -f1 | sort -u > H2A.txt"),
                  paste("grep", '"IPR000558"', "braker.peptide.tsv | cut -f1 | sort -u > H2B.txt"),
                  paste("grep", '"IPR000164"', "braker.peptide.tsv | cut -f1 | sort -u > H3.txt"),
                  paste("grep", '"IPR001951"', "braker.peptide.tsv | cut -f1 | sort -u > H4.txt"),
                  paste("grep", '"IPR005819"', "braker.peptide.tsv | cut -f1 | sort -u > H1.txt"),
                  "grep -w -A 2 -f H1.txt  braker.peptide --no-group-separator > H1",
                  "grep -w -A 2 -f H3.txt  braker.peptide --no-group-separator > H3",
                  "grep -w -A 2 -f H4.txt  braker.peptide --no-group-separator > H4",
                  "grep -w -A 2 -f H2A.txt  braker.peptide --no-group-separator > H2A",
                  "grep -w -A 2 -f H2B.txt  braker.peptide --no-group-separator > H2B",
                  "rm H1.txt H2A.txt H2B.txt H3.txt H4.txt",
                  "mv braker.peptide.gff3 Interpro_annotation.gff3",
                  "mv braker.peptide.json Interpro_annotation.json",
                  "mv braker.peptide.tsv Interpro_annotation.tsv",
                  "mv braker.peptide.xml Interpro_annotation.xml",
                  "mkdir FilesforGUI",
                  "mv Speciesandprotein.txt FilesforGUI/",
                  "mv blastbpstricteval.table FilesforGUI/",
                  "mv InterproIdentified.txt FilesforGUI/",
                  "cp ../../../ListofINTERPRONUMBERS.csv FilesforGUI/",
                  "mv H1 H2A H2B H3 H4 FilesforGUI/",
                  "cp ../../../humanhistones.csv FilesforGUI/",
                  "cp ../../../Interprofunctions.csv FilesforGUI/",
                  "cp ../../../scerhistones.csv FilesforGUI/",
                  "cp ../../../associationTable.csv FilesforGUI/",
                  "mv FilesforGUI/ ../",
                  "cd ../"

                ), fileConn)
              close(fileConn)


              system("chmod +x *.sh")
              withProgress(message = "Annotation in progress, be patient", value=0, detail="0%: Running Canu",
                           {
                             system("./Canu.sh")
                             incProgress(0.05,detail = paste0("5%: Running RepeatMasker"))
                             system("./RepeatMasker.sh")
                             incProgress(0.15,detail = paste0("15%: Running BRAKER3 Pipeline"))
                             system("./braker.sh")
                             incProgress(0.55,detail = paste0("55%: Running Interproscan"))
                             system("./interpro.sh")
                             incProgress(0.80,detail = paste0("80%: Running BLAST+"))
                             system("./blasting.sh")
                             incProgress(0.90,detail = paste0("90%: Parsing Outputs"))
                             setwd("tempforpipeline/interpro")
                             blast<-read.table("nocommentedlines.txt", sep="\t", header=FALSE)
                             limited<-blast[which(blast[,11]<=evalout),]
                             write.table(limited, "blastbpstricteval.table", sep="\t", row.names=FALSE, col.names= FALSE, quote=FALSE)
                             setwd("../../")
                             system("./parsing.sh")
                             incProgress(0.99,detail = paste0("99%: Moving Files To Final Location"))
                             system(paste("mv tempforpipeline",paste0(finallocation,"/",speciesout)))
                             incProgress(1,detail = paste0("100%"))
                           })
            }

            }else(short=="short")
        {        if (speciesout=="" || datalocationout=="" || cladeout=="" || nprocout=="" || evalout=="" || length(finallocation)==0)
        {
          showModal(modalDialog(
            title = "WARNING",
            paste0("All selections not completed"),
            easyClose = TRUE,
            footer = NULL
          ))

        }else{
          vevletpaired<-input$paired
          if(velvetpaired=="No")
          {
          fileConn<-file("Velvet.sh")
          writeLines(
            c("mkdir tempforpipeline",
              "cd tempforpipeline/",
              "mkdir Velvet",
              "cd Velvet/",
              paste("cp",noquote(datalocationout),  "input.fastq"),
              paste("singularity exec ../../Perceptivev0.1.sif velveth assembly 31 -fastq input.fastq"),
              paste("singularity exec ../../Perceptivev0.1.sif velvetg assembly -exp_cov auto -cov_cutoff auto"),
              "cp assembly/contigs.fa ../input.fasta",
              "cd ../"

            ), fileConn)
          close(fileConn)
          }else if (velvetpaired=="Yes")
          {
            datalocationouttwo<<- as.character(parseFilePaths(roots = volumes, input$fastq2)[4])
            fileConn<-file("Velvet.sh")
            writeLines(
              c("mkdir tempforpipeline",
                "cd tempforpipeline/",
                "mkdir Velvet",
                "cd Velvet/",
                paste("cp",noquote(datalocationout),  "input.fastq"),
                paste("cp",noquote(datalocationouttwo),  "input2.fastq"),
                paste("singularity exec ../../Perceptivev0.1.sif velveth assembly 31 -shortPaired -fastq -separate input.fastq input2.fastq"),
                paste("singularity exec ../../Perceptivev0.1.sif velvetg assembly -exp_cov auto -cov_cutoff auto"),
                "cp assembly/contigs.fa ../input.fasta",
                "cd ../"

              ), fileConn)
            close(fileConn)
          }


          fileConn<-file("RepeatMasker.sh")
          writeLines(
            c("cd tempforpipeline/",
              paste("singularity overlay create -s 10000 temp.img"),
              paste("singularity exec --overlay temp.img ../Perceptivev0.1.sif RepeatMasker input.fasta --xsmall -species",cladeout, "-pa",nprocout)


            ), fileConn)
          close(fileConn)

          fileConn<-file("braker.sh")
          writeLines(

            c("cd tempforpipeline/",
              paste("singularity exec ../braker3.sif braker.pl --genome=input.fasta.masked --threads", nprocout, "--prot_seq=../Eukaryota.fa"))
            , fileConn)
          close(fileConn)


          fileConn<-file("interpro.sh")
          writeLines(
            c(
              "cd tempforpipeline/",
              "mkdir interpro",
              "cp braker/braker.aa interpro/",
              "cd interpro/",
              paste("sed", '"s/\\*//g"', "< braker.aa > braker.peptide"),
              "rm braker.aa",
              paste("singularity exec --overlay ../temp.img ../../Perceptivev0.1.sif interproscan.sh -i braker.peptide -cpu", nprocout, "-dp --goterms -appl AntiFam-7.0,CDD-3.20,Coils-2.2.1,FunFam-4.3.0,Gene3D-4.3.0,Hamap-2023_05,NCBIfam-15.0,PANTHER-19.0,Pfam-37.0,PIRSF-3.10,PIRSR-2023_05,PRINTS-42.0,ProSitePatterns-2023_05,ProSiteProfiles-2023_05,SFLD-4,SMART-9.0,SUPERFAMILY-1.75"),
              "rm -r temp/",
              "rm ../temp.img"

            ), fileConn)
          close(fileConn)

          blastexp<-input$uniqueblast
          if(blastexp=="No")
          {
            fileConn<-file("blasting.sh")
            writeLines(
              c("cd tempforpipeline/interpro",
                paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
                "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
                "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
                "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
                "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
                paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db /dependencies/modelorgsprot/modelorgsprot -db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
                paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
              ), fileConn)
            close(fileConn)
          }else if (blastexp=="Yes")
          {
            blastnewdb<<- as.character(parseFilePaths(roots = volumes, input$blastdb)[4])
            fileConn<-file("blasting.sh")
            writeLines(
              c("cd tempforpipeline/interpro",
                paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
                "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
                "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
                "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
                "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
                paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db", blastnewdb, "-db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
                paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
              ), fileConn)
            close(fileConn)
          }

          fileConn<-file("parsing.sh")
          writeLines(
            c("cd tempforpipeline/interpro",
              "cut -f2 blastbpstricteval.table >blastpstrictevalgenes.txt",
              "for line in $(cat blastpstrictevalgenes.txt); do grep -w $line ../../../temp.faa >> Speciesandprotein.txt; done",
              paste("grep", '"IPR002119"', "braker.peptide.tsv | cut -f1 | sort -u > H2A.txt"),
              paste("grep", '"IPR000558"', "braker.peptide.tsv | cut -f1 | sort -u > H2B.txt"),
              paste("grep", '"IPR000164"', "braker.peptide.tsv | cut -f1 | sort -u > H3.txt"),
              paste("grep", '"IPR001951"', "braker.peptide.tsv | cut -f1 | sort -u > H4.txt"),
              paste("grep", '"IPR005819"', "braker.peptide.tsv | cut -f1 | sort -u > H1.txt"),
              "grep -w -A 2 -f H1.txt  braker.peptide --no-group-separator > H1",
              "grep -w -A 2 -f H3.txt  braker.peptide --no-group-separator > H3",
              "grep -w -A 2 -f H4.txt  braker.peptide --no-group-separator > H4",
              "grep -w -A 2 -f H2A.txt  braker.peptide --no-group-separator > H2A",
              "grep -w -A 2 -f H2B.txt  braker.peptide --no-group-separator > H2B",
              "rm H1.txt H2A.txt H2B.txt H3.txt H4.txt",
              "mv braker.peptide.gff3 Interpro_annotation.gff3",
              "mv braker.peptide.json Interpro_annotation.json",
              "mv braker.peptide.tsv Interpro_annotation.tsv",
              "mv braker.peptide.xml Interpro_annotation.xml",
              "mkdir FilesforGUI",
              "mv Speciesandprotein.txt FilesforGUI/",
              "mv blastbpstricteval.table FilesforGUI/",
              "mv InterproIdentified.txt FilesforGUI/",
              "cp ../../../ListofINTERPRONUMBERS.csv FilesforGUI/",
              "mv H1 H2A H2B H3 H4 FilesforGUI/",
              "cp ../../../humanhistones.csv FilesforGUI/",
              "cp ../../../Interprofunctions.csv FilesforGUI/",
              "cp ../../../scerhistones.csv FilesforGUI/",
              "cp ../../../associationTable.csv FilesforGUI/",
              "mv FilesforGUI/ ../",
              "cd ../"

            ), fileConn)
          close(fileConn)


          system("chmod +x *.sh")
          withProgress(message = "Annotation in progress, be patient", value=0, detail="0%: Running Velvet",
                       {
                         system("./Velvet.sh")
                         incProgress(0.05,detail = paste0("5%: Running RepeatMasker"))
                         system("./RepeatMasker.sh")
                         incProgress(0.15,detail = paste0("15%: Running BRAKER3 Pipeline"))
                         system("./braker.sh")
                         incProgress(0.55,detail = paste0("55%: Running Interproscan"))
                         system("./interpro.sh")
                         incProgress(0.80,detail = paste0("80%: Running BLAST+"))
                         system("./blasting.sh")
                         incProgress(0.90,detail = paste0("90%: Parsing Outputs"))
                         setwd("tempforpipeline/interpro")
                         blast<-read.table("nocommentedlines.txt", sep="\t", header=FALSE)
                         limited<-blast[which(blast[,11]<=evalout),]
                         write.table(limited, "blastbpstricteval.table", sep="\t", row.names=FALSE, col.names= FALSE, quote=FALSE)
                         setwd("../../")
                         system("./parsing.sh")
                         incProgress(0.99,detail = paste0("99%: Moving Files To Final Location"))
                         system(paste("mv tempforpipeline",paste0(finallocation,"/",speciesout)))
                         incProgress(1,detail = paste0("100%"))
                       })
        }}
        }else if(novode=="Yes" && RNAeval=="Yes")
        {
          datalocationout<<- as.character(parseFilePaths(roots = volumes, input$fastq)[4])
          R1<<-as.character(parseFilePaths(roots = volumes, input$fastqR1RNA)[4])
          R2<<-as.character(parseFilePaths(roots = volumes, input$fastqR2RNA)[4])
          short<<-input$length
          if(short=="Long")
          {
            technology<<-input$tech
            gsize<<-input$genomesize
            if (speciesout=="" || datalocationout=="" || cladeout=="" || nprocout=="" || evalout=="" || length(finallocation)==0 || technology=="" || gsize=="")
            {
              showModal(modalDialog(
                title = "WARNING",
                paste0("All selections not completed"),
                easyClose = TRUE,
                footer = NULL
              ))

            }else{
              if(technology=="Nanopore")
              {tnolog<<-"-nanopore-raw"
              }else if (technology=="Pacbio")
              {
                tnolog<<-"-pacbio"
              }else if (technology=="Pacbio HiFi")
              {tnolog<<-"-pacbio-hifi"}
              fileConn<-file("Canu.sh")
              canuarguments<-input$passthrough
              if(canuarguments=="No")
              {
                fileConn<-file("Canu.sh")
                writeLines(
                  c("mkdir tempforpipeline",
                    "cd tempforpipeline/",
                    paste("cp",noquote(datalocationout),  "input.fastq"),
                    paste0("singularity exec ../Perceptivev0.1.sif canu -d assembly -p assembled genomesize=",gsize, " ", tnolog, " input.fastq"),
                    "cp assembly/assembled.contig.fasta input.fasta"

                  ), fileConn)
                close(fileConn)
              }else if (canuarguments=="Yes")
              {
                fileConn<-file("Canu.sh")
                writeLines(
                  c("mkdir tempforpipeline",
                    "cd tempforpipeline/",
                    paste("cp",noquote(datalocationout),  "input.fastq"),
                    paste0("singularity exec ../Perceptivev0.1.sif canu -d assembly -p assembled genomesize=",gsize, " ", tnolog, " input.fastq", canuargs),
                    "cp assembly/assembled.contig.fasta input.fasta"

                  ), fileConn)
                close(fileConn)
              }

              fileConn<-file("RepeatMasker.sh")
              writeLines(
                c("cd tempforpipeline/",
                  paste("singularity overlay create -s 10000 temp.img"),
                  paste("singularity exec --overlay temp.img ../Perceptivev0.1.sif RepeatMasker input.fasta --xsmall -species",cladeout, "-pa",nprocout)


                ), fileConn)
              close(fileConn)

              fileConn<-file("braker.sh")
              writeLines(

                c("cd tempforpipeline/",
                  paste("singularity exec ../braker3.sif braker.pl --genome=input.fasta.masked --threads", nprocout, "--rnaseq_sets_ids=ID1 --rnaseq_sets_dirs=../"))
                , fileConn)
              close(fileConn)


              fileConn<-file("interpro.sh")
              writeLines(
                c(
                  "cd tempforpipeline/",
                  "mkdir interpro",
                  "cp braker/braker.aa interpro/",
                  "cd interpro/",
                  paste("sed", '"s/\\*//g"', "< braker.aa > braker.peptide"),
                  "rm braker.aa",
                  paste("singularity exec --overlay ../temp.img ../../Perceptivev0.1.sif interproscan.sh -i braker.peptide -cpu", nprocout, "-dp --goterms -appl AntiFam-7.0,CDD-3.20,Coils-2.2.1,FunFam-4.3.0,Gene3D-4.3.0,Hamap-2023_05,NCBIfam-15.0,PANTHER-19.0,Pfam-37.0,PIRSF-3.10,PIRSR-2023_05,PRINTS-42.0,ProSitePatterns-2023_05,ProSiteProfiles-2023_05,SFLD-4,SMART-9.0,SUPERFAMILY-1.75"),
                  "rm -r temp/",
                  "rm ../temp.img"

                ), fileConn)
              close(fileConn)

              blastexp<-input$uniqueblast
              if(blastexp=="No")
              {
                fileConn<-file("blasting.sh")
                writeLines(
                  c("cd tempforpipeline/interpro",
                    paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
                    "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
                    "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
                    "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
                    "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
                    paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db /dependencies/modelorgsprot/modelorgsprot -db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
                    paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
                  ), fileConn)
                close(fileConn)
              }else if (blastexp=="Yes")
              {
                blastnewdb<<- as.character(parseFilePaths(roots = volumes, input$blastdb)[4])
                fileConn<-file("blasting.sh")
                writeLines(
                  c("cd tempforpipeline/interpro",
                    paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
                    "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
                    "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
                    "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
                    "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
                    paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db", blastnewdb, "-db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
                    paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
                  ), fileConn)
                close(fileConn)
              }

              fileConn<-file("parsing.sh")
              writeLines(
                c("cd tempforpipeline/interpro",
                  "cut -f2 blastbpstricteval.table >blastpstrictevalgenes.txt",
                  "for line in $(cat blastpstrictevalgenes.txt); do grep -w $line ../../../temp.faa >> Speciesandprotein.txt; done",
                  paste("grep", '"IPR002119"', "braker.peptide.tsv | cut -f1 | sort -u > H2A.txt"),
                  paste("grep", '"IPR000558"', "braker.peptide.tsv | cut -f1 | sort -u > H2B.txt"),
                  paste("grep", '"IPR000164"', "braker.peptide.tsv | cut -f1 | sort -u > H3.txt"),
                  paste("grep", '"IPR001951"', "braker.peptide.tsv | cut -f1 | sort -u > H4.txt"),
                  paste("grep", '"IPR005819"', "braker.peptide.tsv | cut -f1 | sort -u > H1.txt"),
                  "grep -w -A 2 -f H1.txt  braker.peptide --no-group-separator > H1",
                  "grep -w -A 2 -f H3.txt  braker.peptide --no-group-separator > H3",
                  "grep -w -A 2 -f H4.txt  braker.peptide --no-group-separator > H4",
                  "grep -w -A 2 -f H2A.txt  braker.peptide --no-group-separator > H2A",
                  "grep -w -A 2 -f H2B.txt  braker.peptide --no-group-separator > H2B",
                  "rm H1.txt H2A.txt H2B.txt H3.txt H4.txt",
                  "mv braker.peptide.gff3 Interpro_annotation.gff3",
                  "mv braker.peptide.json Interpro_annotation.json",
                  "mv braker.peptide.tsv Interpro_annotation.tsv",
                  "mv braker.peptide.xml Interpro_annotation.xml",
                  "mkdir FilesforGUI",
                  "mv Speciesandprotein.txt FilesforGUI/",
                  "mv blastbpstricteval.table FilesforGUI/",
                  "mv InterproIdentified.txt FilesforGUI/",
                  "cp ../../../ListofINTERPRONUMBERS.csv FilesforGUI/",
                  "mv H1 H2A H2B H3 H4 FilesforGUI/",
                  "cp ../../../humanhistones.csv FilesforGUI/",
                  "cp ../../../Interprofunctions.csv FilesforGUI/",
                  "cp ../../../scerhistones.csv FilesforGUI/",
                  "cp ../../../associationTable.csv FilesforGUI/",
                  "mv FilesforGUI/ ../",
                  "cd ../"

                ), fileConn)
              close(fileConn)


              system("chmod +x *.sh")
              withProgress(message = "Annotation in progress, be patient", value=0, detail="0%: Running Canu",
                           {
                             system("./Canu.sh")
                             incProgress(0.05,detail = paste0("5%: Running RepeatMasker"))
                             system("./RepeatMasker.sh")
                             incProgress(0.05,detail = paste0("5%: Running BRAKER3 Pipeline"))
                             if(R2=="")
                             {
                               system(paste("cp", R1, "ID1.fastq"))
                             }else
                             {
                               system(paste("cp", R1, "ID1_1.fastq"))
                               system(paste("cp", R2, "ID1_2.fastq"))
                             }
                             system("./braker.sh")
                             incProgress(0.55,detail = paste0("55%: Running Interproscan"))
                             system("./interpro.sh")
                             incProgress(0.80,detail = paste0("80%: Running BLAST+"))
                             system("./blasting.sh")
                             incProgress(0.90,detail = paste0("90%: Parsing Outputs"))
                             setwd("tempforpipeline/interpro")
                             blast<-read.table("nocommentedlines.txt", sep="\t", header=FALSE)
                             limited<-blast[which(blast[,11]<=evalout),]
                             write.table(limited, "blastbpstricteval.table", sep="\t", row.names=FALSE, col.names= FALSE, quote=FALSE)
                             setwd("../../")
                             system("./parsing.sh")
                             incProgress(0.99,detail = paste0("99%: Moving Files To Final Location"))
                             system(paste("mv tempforpipeline",paste0(finallocation,"/",speciesout)))
                             incProgress(1,detail = paste0("100%"))
                           })
            }

          }else(short=="short")
          {        if (speciesout=="" || datalocationout=="" || cladeout=="" || nprocout=="" || evalout=="" || length(finallocation)==0)
          {
            showModal(modalDialog(
              title = "WARNING",
              paste0("All selections not completed"),
              easyClose = TRUE,
              footer = NULL
            ))

          }else{

            vevletpaired<-input$paired
            if(velvetpaired=="No")
            {
              fileConn<-file("Velvet.sh")
              writeLines(
                c("mkdir tempforpipeline",
                  "cd tempforpipeline/",
                  "mkdir Velvet",
                  "cd Velvet/",
                  paste("cp",noquote(datalocationout),  "input.fastq"),
                  paste("singularity exec ../../Perceptivev0.1.sif velveth assembly 31 -fastq input.fastq"),
                  paste("singularity exec ../../Perceptivev0.1.sif velvetg assembly -exp_cov auto -cov_cutoff auto"),
                  "cp assembly/contigs.fa ../input.fasta",
                  "cd ../"

                ), fileConn)
              close(fileConn)
            }else if (velvetpaired=="Yes")
            {
              datalocationouttwo<<- as.character(parseFilePaths(roots = volumes, input$fastq2)[4])
              fileConn<-file("Velvet.sh")
              writeLines(
                c("mkdir tempforpipeline",
                  "cd tempforpipeline/",
                  "mkdir Velvet",
                  "cd Velvet/",
                  paste("cp",noquote(datalocationout),  "input.fastq"),
                  paste("cp",noquote(datalocationouttwo),  "input2.fastq"),
                  paste("singularity exec ../../Perceptivev0.1.sif velveth assembly 31 -shortPaired -fastq -separate input.fastq input2.fastq"),
                  paste("singularity exec ../../Perceptivev0.1.sif velvetg assembly -exp_cov auto -cov_cutoff auto"),
                  "cp assembly/contigs.fa ../input.fasta",
                  "cd ../"

                ), fileConn)
              close(fileConn)
            }

            fileConn<-file("RepeatMasker.sh")
            writeLines(
              c("cd tempforpipeline/",
                paste("singularity overlay create -s 10000 temp.img"),
                paste("singularity exec --overlay temp.img ../Perceptivev0.1.sif RepeatMasker input.fasta --xsmall -species",cladeout, "-pa",nprocout)


              ), fileConn)
            close(fileConn)

            fileConn<-file("braker.sh")
            writeLines(

              c("cd tempforpipeline/",
                paste("singularity exec ../braker3.sif braker.pl --genome=input.fasta.masked --threads", nprocout, "--rnaseq_sets_ids=ID1 --rnaseq_sets_dirs=../"))
              , fileConn)
            close(fileConn)


            fileConn<-file("interpro.sh")
            writeLines(
              c(
                "cd tempforpipeline/",
                "mkdir interpro",
                "cp braker/braker.aa interpro/",
                "cd interpro/",
                paste("sed", '"s/\\*//g"', "< braker.aa > braker.peptide"),
                "rm braker.aa",
                paste("singularity exec --overlay ../temp.img ../../Perceptivev0.1.sif interproscan.sh -i braker.peptide -cpu", nprocout, "-dp --goterms -appl AntiFam-7.0,CDD-3.20,Coils-2.2.1,FunFam-4.3.0,Gene3D-4.3.0,Hamap-2023_05,NCBIfam-15.0,PANTHER-19.0,Pfam-37.0,PIRSF-3.10,PIRSR-2023_05,PRINTS-42.0,ProSitePatterns-2023_05,ProSiteProfiles-2023_05,SFLD-4,SMART-9.0,SUPERFAMILY-1.75"),
                "rm -r temp/",
                "rm ../temp.img"

              ), fileConn)
            close(fileConn)

            blastexp<-input$uniqueblast
            if(blastexp=="No")
            {
              fileConn<-file("blasting.sh")
              writeLines(
                c("cd tempforpipeline/interpro",
                  paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
                  "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
                  "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
                  "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
                  "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
                  paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db /dependencies/modelorgsprot/modelorgsprot -db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
                  paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
                ), fileConn)
              close(fileConn)
            }else if (blastexp=="Yes")
            {
              blastnewdb<<- as.character(parseFilePaths(roots = volumes, input$blastdb)[4])
              fileConn<-file("blasting.sh")
              writeLines(
                c("cd tempforpipeline/interpro",
                  paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
                  "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
                  "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
                  "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
                  "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
                  paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db", blastnewdb, "-db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
                  paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
                ), fileConn)
              close(fileConn)
            }

            fileConn<-file("parsing.sh")
            writeLines(
              c("cd tempforpipeline/interpro",
                "cut -f2 blastbpstricteval.table >blastpstrictevalgenes.txt",
                "for line in $(cat blastpstrictevalgenes.txt); do grep -w $line ../../../temp.faa >> Speciesandprotein.txt; done",
                paste("grep", '"IPR002119"', "braker.peptide.tsv | cut -f1 | sort -u > H2A.txt"),
                paste("grep", '"IPR000558"', "braker.peptide.tsv | cut -f1 | sort -u > H2B.txt"),
                paste("grep", '"IPR000164"', "braker.peptide.tsv | cut -f1 | sort -u > H3.txt"),
                paste("grep", '"IPR001951"', "braker.peptide.tsv | cut -f1 | sort -u > H4.txt"),
                paste("grep", '"IPR005819"', "braker.peptide.tsv | cut -f1 | sort -u > H1.txt"),
                "grep -w -A 2 -f H1.txt  braker.peptide --no-group-separator > H1",
                "grep -w -A 2 -f H3.txt  braker.peptide --no-group-separator > H3",
                "grep -w -A 2 -f H4.txt  braker.peptide --no-group-separator > H4",
                "grep -w -A 2 -f H2A.txt  braker.peptide --no-group-separator > H2A",
                "grep -w -A 2 -f H2B.txt  braker.peptide --no-group-separator > H2B",
                "rm H1.txt H2A.txt H2B.txt H3.txt H4.txt",
                "mv braker.peptide.gff3 Interpro_annotation.gff3",
                "mv braker.peptide.json Interpro_annotation.json",
                "mv braker.peptide.tsv Interpro_annotation.tsv",
                "mv braker.peptide.xml Interpro_annotation.xml",
                "mkdir FilesforGUI",
                "mv Speciesandprotein.txt FilesforGUI/",
                "mv blastbpstricteval.table FilesforGUI/",
                "mv InterproIdentified.txt FilesforGUI/",
                "cp ../../../ListofINTERPRONUMBERS.csv FilesforGUI/",
                "mv H1 H2A H2B H3 H4 FilesforGUI/",
                "cp ../../../humanhistones.csv FilesforGUI/",
                "cp ../../../Interprofunctions.csv FilesforGUI/",
                "cp ../../../scerhistones.csv FilesforGUI/",
                "cp ../../../associationTable.csv FilesforGUI/",
                "mv FilesforGUI/ ../",
                "cd ../"

              ), fileConn)
            close(fileConn)


            system("chmod +x *.sh")
            withProgress(message = "Annotation in progress, be patient", value=0, detail="0%: Running Velvet",
                         {
                           system("./Velvet.sh")
                           incProgress(0.05,detail = paste0("5%: Running RepeatMasker"))
                           system("./RepeatMasker.sh")
                           incProgress(0.05,detail = paste0("5%: Running BRAKER3 Pipeline"))
                           if(R2=="")
                           {
                             system(paste("cp", R1, "ID1.fastq"))
                           }else
                           {
                             system(paste("cp", R1, "ID1_1.fastq"))
                             system(paste("cp", R2, "ID1_2.fastq"))
                           }
                           system("./braker.sh")
                           incProgress(0.15,detail = paste0("15%: Running BRAKER3 Pipeline"))
                           system("./braker.sh")
                           incProgress(0.55,detail = paste0("55%: Running Interproscan"))
                           system("./interpro.sh")
                           incProgress(0.80,detail = paste0("80%: Running BLAST+"))
                           system("./blasting.sh")
                           incProgress(0.90,detail = paste0("90%: Parsing Outputs"))
                           setwd("tempforpipeline/interpro")
                           blast<-read.table("nocommentedlines.txt", sep="\t", header=FALSE)
                           limited<-blast[which(blast[,11]<=evalout),]
                           write.table(limited, "blastbpstricteval.table", sep="\t", row.names=FALSE, col.names= FALSE, quote=FALSE)
                           setwd("../../")
                           system("./parsing.sh")
                           incProgress(0.99,detail = paste0("99%: Moving Files To Final Location"))
                           system(paste("mv tempforpipeline",paste0(finallocation,"/",speciesout)))
                           incProgress(1,detail = paste0("100%"))
                         })
          }}
        }else if (novode=="No" && RNAeval=="Yes")
        {
        datalocationout<<- as.character(parseFilePaths(roots = volumes, input$fasta)[4])
        R1<<-as.character(parseFilePaths(roots = volumes, input$fastqR1RNA)[4])
        R2<<-as.character(parseFilePaths(roots = volumes, input$fastqR2RNA)[4])
        if (speciesout=="" || datalocationout=="" || cladeout=="" || nprocout=="" || evalout=="" || length(finallocation)==0 || R1=="")
        {
          showModal(modalDialog(
            title = "WARNING",
            paste0("All selections not completed"),
            easyClose = TRUE,
            footer = NULL
          ))

        }else{
          fileConn<-file("RepeatMasker.sh")
          writeLines(
            c("mkdir tempforpipeline",
              "cd tempforpipeline/",
              paste("cp",noquote(datalocationout),  "input.fasta"),
              paste("singularity overlay create -s 10000 temp.img"),
              paste("singularity exec --overlay temp.img ../Perceptivev0.1.sif RepeatMasker input.fasta --xsmall -species",cladeout, "-pa",nprocout)


            ), fileConn)
          close(fileConn)



          fileConn<-file("braker.sh")
          writeLines(

            c("cd tempforpipeline/",
              paste("singularity exec ../braker3.sif braker.pl --genome=input.fasta.masked --threads", nprocout, "--rnaseq_sets_ids=ID1 --rnaseq_sets_dirs=../"))
            , fileConn)
          close(fileConn)


          fileConn<-file("interpro.sh")
          writeLines(
            c(
              "cd tempforpipeline/",
              "mkdir interpro",
              "cp braker/braker.aa interpro/",
              "cd interpro/",
              paste("sed", '"s/\\*//g"', "< braker.aa > braker.peptide"),
              "rm braker.aa",
              paste("singularity exec --overlay ../temp.img ../../Perceptivev0.1.sif interproscan.sh -i braker.peptide -cpu", nprocout, "-dp --goterms -appl AntiFam-7.0,CDD-3.20,Coils-2.2.1,FunFam-4.3.0,Gene3D-4.3.0,Hamap-2023_05,NCBIfam-15.0,PANTHER-19.0,Pfam-37.0,PIRSF-3.10,PIRSR-2023_05,PRINTS-42.0,ProSitePatterns-2023_05,ProSiteProfiles-2023_05,SFLD-4,SMART-9.0,SUPERFAMILY-1.75"),
              "rm -r temp/",
              "rm ../temp.img"

            ), fileConn)
          close(fileConn)

          blastexp<-input$uniqueblast
          if(blastexp=="No")
          {
          fileConn<-file("blasting.sh")
          writeLines(
            c("cd tempforpipeline/interpro",
              paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
              "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
              "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
              "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
              "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
              paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db /dependencies/modelorgsprot/modelorgsprot -db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
              paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
            ), fileConn)
          close(fileConn)
          }else if (blastexp=="Yes")
          {
            blastnewdb<<- as.character(parseFilePaths(roots = volumes, input$blastdb)[4])
            fileConn<-file("blasting.sh")
            writeLines(
              c("cd tempforpipeline/interpro",
                paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
                "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
                "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
                "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
                "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
                paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db", blastnewdb, "-db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
                paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
              ), fileConn)
            close(fileConn)
          }

          fileConn<-file("parsing.sh")
          writeLines(
            c("cd tempforpipeline/interpro",
              "cut -f2 blastbpstricteval.table >blastpstrictevalgenes.txt",
              "for line in $(cat blastpstrictevalgenes.txt); do grep -w $line ../../../temp.faa >> Speciesandprotein.txt; done",
              paste("grep", '"IPR002119"', "braker.peptide.tsv | cut -f1 | sort -u > H2A.txt"),
              paste("grep", '"IPR000558"', "braker.peptide.tsv | cut -f1 | sort -u > H2B.txt"),
              paste("grep", '"IPR000164"', "braker.peptide.tsv | cut -f1 | sort -u > H3.txt"),
              paste("grep", '"IPR001951"', "braker.peptide.tsv | cut -f1 | sort -u > H4.txt"),
              paste("grep", '"IPR005819"', "braker.peptide.tsv | cut -f1 | sort -u > H1.txt"),
              "grep -w -A 2 -f H1.txt  braker.peptide --no-group-separator > H1",
              "grep -w -A 2 -f H3.txt  braker.peptide --no-group-separator > H3",
              "grep -w -A 2 -f H4.txt  braker.peptide --no-group-separator > H4",
              "grep -w -A 2 -f H2A.txt  braker.peptide --no-group-separator > H2A",
              "grep -w -A 2 -f H2B.txt  braker.peptide --no-group-separator > H2B",
              "rm H1.txt H2A.txt H2B.txt H3.txt H4.txt",
              "mv braker.peptide.gff3 Interpro_annotation.gff3",
              "mv braker.peptide.json Interpro_annotation.json",
              "mv braker.peptide.tsv Interpro_annotation.tsv",
              "mv braker.peptide.xml Interpro_annotation.xml",
              "mkdir FilesforGUI",
              "mv Speciesandprotein.txt FilesforGUI/",
              "mv blastbpstricteval.table FilesforGUI/",
              "mv InterproIdentified.txt FilesforGUI/",
              "cp ../../../ListofINTERPRONUMBERS.csv FilesforGUI/",
              "mv H1 H2A H2B H3 H4 FilesforGUI/",
              "cp ../../../humanhistones.csv FilesforGUI/",
              "cp ../../../Interprofunctions.csv FilesforGUI/",
              "cp ../../../scerhistones.csv FilesforGUI/",
              "cp ../../../associationTable.csv FilesforGUI/",
              "mv FilesforGUI/ ../",
              "cd ../"

            ), fileConn)
          close(fileConn)


          system("chmod +x *.sh")
          withProgress(message = "Annotation in progress, be patient", value=0, detail="0%: Running RepeatMasker",
                       {
                         system("./RepeatMasker.sh")
                         incProgress(0.05,detail = paste0("5%: Running BRAKER3 Pipeline"))
                         if(R2=="")
                         {
                           system(paste("cp", R1, "ID1.fastq"))
                         }else
                         {
                           system(paste("cp", R1, "ID1_1.fastq"))
                           system(paste("cp", R2, "ID1_2.fastq"))
                         }
                         system("./braker.sh")
                         system(paste("rm *.fastq"))
                         incProgress(0.55,detail = paste0("55%: Running Interproscan"))
                         system("./interpro.sh")
                         incProgress(0.80,detail = paste0("80%: Running BLAST+"))
                         system("./blasting.sh")
                         incProgress(0.90,detail = paste0("90%: Parsing Outputs"))
                         setwd("tempforpipeline/interpro")
                         blast<-read.table("nocommentedlines.txt", sep="\t", header=FALSE)
                         limited<-blast[which(blast[,11]<=evalout),]
                         write.table(limited, "blastbpstricteval.table", sep="\t", row.names=FALSE, col.names= FALSE, quote=FALSE)
                         setwd("../../")
                         system("./parsing.sh")
                         incProgress(0.99,detail = paste0("99%: Moving Files To Final Location"))
                         system(paste("mv tempforpipeline",paste0(finallocation,"/",speciesout)))
                         incProgress(1,detail = paste0("100%"))
                       })
        }
        }else if (novode=="No" && RNAeval=="No")
        {
        datalocationout<<- as.character(parseFilePaths(roots = volumes, input$fasta)[4])
        if (speciesout=="" || datalocationout=="" || cladeout=="" || nprocout=="" || evalout=="" || length(finallocation)==0)
        {
          showModal(modalDialog(
            title = "WARNING",
            paste0("All selections not completed"),
            easyClose = TRUE,
            footer = NULL
          ))

        }else{
          fileConn<-file("RepeatMasker.sh")
          writeLines(
            c("mkdir tempforpipeline",
              "cd tempforpipeline/",
              paste("cp",noquote(datalocationout),  "input.fasta"),
              paste("singularity overlay create -s 10000 temp.img"),
              paste("singularity exec --overlay temp.img ../Perceptivev0.1.sif RepeatMasker input.fasta --xsmall -species",cladeout, "-pa",nprocout)


            ), fileConn)
          close(fileConn)

          fileConn<-file("braker.sh")
          writeLines(

            c("cd tempforpipeline/",
              paste("singularity exec ../braker3.sif braker.pl --genome=input.fasta.masked --threads", nprocout, "--prot_seq=../Eukaryota.fa"))
            , fileConn)
          close(fileConn)


          fileConn<-file("interpro.sh")
          writeLines(
            c(
              "cd tempforpipeline/",
              "mkdir interpro",
              "cp braker/braker.aa interpro/",
              "cd interpro/",
              paste("sed", '"s/\\*//g"', "< braker.aa > braker.peptide"),
              "rm braker.aa",
              paste("singularity exec --overlay ../temp.img ../../Perceptivev0.1.sif interproscan.sh -i braker.peptide -cpu", nprocout, "-dp --goterms -appl AntiFam-7.0,CDD-3.20,Coils-2.2.1,FunFam-4.3.0,Gene3D-4.3.0,Hamap-2023_05,NCBIfam-15.0,PANTHER-19.0,Pfam-37.0,PIRSF-3.10,PIRSR-2023_05,PRINTS-42.0,ProSitePatterns-2023_05,ProSiteProfiles-2023_05,SFLD-4,SMART-9.0,SUPERFAMILY-1.75"),
              "rm -r temp/",
              "rm ../temp.img"

            ), fileConn)
          close(fileConn)
          blastexp<-input$uniqueblast
          if(blastexp=="No")
          {
          fileConn<-file("blasting.sh")
          writeLines(
            c("cd tempforpipeline/interpro",
              paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
              "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
              "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
              "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
              "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
              paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db /dependencies/modelorgsprot/modelorgsprot -db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
              paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
            ), fileConn)
          close(fileConn)
          }else if (blastexp=="Yes")
          {
            blastnewdb<<- as.character(parseFilePaths(roots = volumes, input$blastdb)[4])
            fileConn<-file("blasting.sh")
            writeLines(
              c("cd tempforpipeline/interpro",
                paste("cat ../../../ListofINTERPRONUMBERS.csv | cut -d", as.character('","'), "-f1 | sort -u > justIPR"),
                "for line in $(cat justIPR); do grep -w $line braker.peptide.tsv >> InterproIdentified.txt; done",
                "cut -f1 InterproIdentified.txt> InterproIdentifiedjustgene.txt",
                "cat InterproIdentifiedjustgene.txt | uniq > InterproIdentifiedjustgeneuniq.txt",
                "grep -w -A 2 -f  InterproIdentifiedjustgeneuniq.txt braker.peptide --no-group-separator > peptides_interpro.faa",
                paste("
                singularity exec ../../Perceptivev0.1.sif blastp -query peptides_interpro.faa -db", blastnewdb, "-db_soft_mask 21 -outfmt 7 -out blastx_outfmt7_results.out -num_threads", nprocout),
                paste("grep -v", as.character('"#"'), "blastx_outfmt7_results.out > nocommentedlines.txt")
              ), fileConn)
            close(fileConn)
          }

          fileConn<-file("parsing.sh")
          writeLines(
            c("cd tempforpipeline/interpro",
              "cut -f2 blastbpstricteval.table >blastpstrictevalgenes.txt",
              "for line in $(cat blastpstrictevalgenes.txt); do grep -w $line ../../../temp.faa >> Speciesandprotein.txt; done",
              paste("grep", '"IPR002119"', "braker.peptide.tsv | cut -f1 | sort -u > H2A.txt"),
              paste("grep", '"IPR000558"', "braker.peptide.tsv | cut -f1 | sort -u > H2B.txt"),
              paste("grep", '"IPR000164"', "braker.peptide.tsv | cut -f1 | sort -u > H3.txt"),
              paste("grep", '"IPR001951"', "braker.peptide.tsv | cut -f1 | sort -u > H4.txt"),
              paste("grep", '"IPR005819"', "braker.peptide.tsv | cut -f1 | sort -u > H1.txt"),
              "grep -w -A 2 -f H1.txt  braker.peptide --no-group-separator > H1",
              "grep -w -A 2 -f H3.txt  braker.peptide --no-group-separator > H3",
              "grep -w -A 2 -f H4.txt  braker.peptide --no-group-separator > H4",
              "grep -w -A 2 -f H2A.txt  braker.peptide --no-group-separator > H2A",
              "grep -w -A 2 -f H2B.txt  braker.peptide --no-group-separator > H2B",
              "rm H1.txt H2A.txt H2B.txt H3.txt H4.txt",
              "mv braker.peptide.gff3 Interpro_annotation.gff3",
              "mv braker.peptide.json Interpro_annotation.json",
              "mv braker.peptide.tsv Interpro_annotation.tsv",
              "mv braker.peptide.xml Interpro_annotation.xml",
              "mkdir FilesforGUI",
              "mv Speciesandprotein.txt FilesforGUI/",
              "mv blastbpstricteval.table FilesforGUI/",
              "mv InterproIdentified.txt FilesforGUI/",
              "cp ../../../ListofINTERPRONUMBERS.csv FilesforGUI/",
              "mv H1 H2A H2B H3 H4 FilesforGUI/",
              "cp ../../../humanhistones.csv FilesforGUI/",
              "cp ../../../Interprofunctions.csv FilesforGUI/",
              "cp ../../../scerhistones.csv FilesforGUI/",
              "cp ../../../associationTable.csv FilesforGUI/",
              "mv FilesforGUI/ ../",
              "cd ../"

            ), fileConn)
          close(fileConn)


          system("chmod +x *.sh")
          withProgress(message = "Annotation in progress, be patient", value=0, detail="0%: Running RepeatMasker",
                       {
                         system("./RepeatMasker.sh")
                         incProgress(0.05,detail = paste0("5%: Running BRAKER3 Pipeline"))
                         system("./braker.sh")
                         incProgress(0.55,detail = paste0("55%: Running Interproscan"))
                         system("./interpro.sh")
                         incProgress(0.80,detail = paste0("80%: Running BLAST+"))
                         system("./blasting.sh")
                         incProgress(0.90,detail = paste0("90%: Parsing Outputs"))
                         setwd("tempforpipeline/interpro")
                         blast<-read.table("nocommentedlines.txt", sep="\t", header=FALSE)
                         limited<-blast[which(blast[,11]<=evalout),]
                         write.table(limited, "blastbpstricteval.table", sep="\t", row.names=FALSE, col.names= FALSE, quote=FALSE)
                         setwd("../../")
                         system("./parsing.sh")
                         incProgress(0.99,detail = paste0("99%: Moving Files To Final Location"))
                         system(paste("mv tempforpipeline",paste0(finallocation,"/",speciesout)))
                         incProgress(1,detail = paste0("100%"))
                       })
        }
        }






    })


  }



  ui <- dashboardPage(scrollToTop = TRUE, skin="green",
                      #freshTheme = theme,
                      header = dashboardHeader(title="PERCEPTIVE Helper", controlbarIcon =icon("moon")),
                      sidebar = dashboardSidebar(sidebarMenu(id = "inTabset",
                        menuItem("Select Data", tabName = "DataInput", icon = icon("database")),
                        menuItem("Environment Settings", tabName = "settings", icon = icon("gear"))
                      )),
                      tags$head(
                        tags$style(
                          HTML(".shiny-notification {
           height: 200px;
           width: 800px;
           position:fixed;
           top: calc(50% - 50px);
           left: calc(50% - 400px);
           font-size: 200%;
           text-align: center;
           }
           "
                          ,".shiny-notification-close {display: none}")
                        )
                      ),
                      body = dashboardBody(

                        #use_theme(mytheme),
                        tabItems(
                          # First tab content
                          tabItem(tabName = "DataInput",
                                  div( h5(HTML("<b>Operating System Detected: </b>"), textOutput("OS"))),
                                  div( h5(HTML("<b>Select location to save output files: </b>"))),
                                  suppressWarnings (shinyDirButton('endlocation', 'Select Folder' , 'Select location to save files:', multiple = FALSE,
                                                   buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder"))),
                                  div( h5(HTML("<b>Select computation resources: </b>"))),
                                  div( h5(HTML(paste0("Select # cores to use, recommended is 90% of available resources. <b>", detectCores(), ": cores detected. </b>")))),
                                  textInput("numcore", "Number of cores:", value = "", width = NULL, placeholder = NULL),
                                  div( h5(HTML("<b>What is your species name? Please only use alphanumeric characters and no spaces. For example H_sapiens, Drosophila_melanogaster, or CR1998 are appropriate inputs. </b>"))),
                                  textInput("speciesname", "Species name:", value = "", width = NULL, placeholder = NULL),
                                  div( h5(HTML("<b>Do you want to perform <i> de novo </i> genome assembly (beta)? </b>"))),
                                  radioButtons(
                                    "denovo", "Select an option", c("Yes", "No"), inline = TRUE, selected="No"
                                  ),
                                  htmlOutput("needapath", placeholder = TRUE),
                                  conditionalPanel(condition="input.denovo == 'No'", suppressWarnings (shinyFilesButton('fasta', 'Select FASTA' , 'Select pathway to fasta:', multiple = FALSE,
                                                 buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder")))),
                                  conditionalPanel(condition="input.denovo == 'Yes'",

                                                   div( h5(HTML("<b>Are your raw reads long (ONT etc.) or short (Illumina etc.)? </b>"))),
                                                   radioButtons(
                                                     "length", "Select an option", c("Long", "Short"), inline = TRUE
                                                   ),
                                                   conditionalPanel(condition = "input.length =='Short'",
                                                                    div(h5(HTML("Is your input paired ended short reads?"))),
                                                                    radioButtons(
                                                                      "paired", "Select an option", c("No", "Yes"), inline = TRUE, selected="No"
                                                                    ),
                                                                    conditionalPanel(condition = "input.paired =='No'",
                                                                                     div(h5(HTML("<b>Expected file is unzipped fastq from Illumina sequencer. PERCEPTIVE will not use gzipped fastq.</b>"))),
                                                                                     suppressWarnings (shinyFilesButton('fastq', 'Select FASTQ' , 'Select pathway to fastq:', multiple = FALSE,
                                                                                                                        buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder")))
                                                                    ),
                                                                    conditionalPanel(condition = "input.paired =='Yes'",
                                                                                     div( h5(HTML("<b>Expected files is unzipped fastq from Illumina sequencer. PERCEPTIVE will not use gzipped fastq.</b>"))),
                                                                                     suppressWarnings (shinyFilesButton('fastq', 'Select FASTQ R1' , 'Select pathway to fastq R1:', multiple = FALSE,
                                                                                                                        buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder"))),
                                                                                     suppressWarnings (shinyFilesButton('fastq2', 'Select FASTQ R2' , 'Select pathway to fastq R2:', multiple = FALSE,
                                                                                                                        buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder"))),
                                                                    ),
                                                                    ),
                                                   conditionalPanel(condition = "input.length =='Long'",
                                                                    div( h5(HTML("<b>Expected file is unzipped fastq from Nanopore or Pacbio sequencer. PERCEPTIVE will not use gzipped fastq.  </b> Canu expects Nanopore/Pacbio reads to be untrimmed and uncorrected (however corrected or trimmed reads are fine). With respect to HiFi reads, Canu expects corrected and trimmed reads and may not provide expected results with untrimmed/uncorrected reads."))),
                                                                    suppressWarnings (shinyFilesButton('fastq', 'Select FASTQ' , 'Select pathway to fastq:', multiple = FALSE,
                                                                                                       buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder"))),
                                                                    radioButtons(
                                                                      "tech", "Select a long read technology", c("Nanopore", "Pacbio", "Pacbio HiFi"), inline = TRUE
                                                                    ),
                                                                    div(h5(HTML("<b>What is the estimated size of your genome in bases? Use jellyfish or GenomeSource for estimates, or the size of a closely related genome. </b>"))),
                                                                    textInput("genomesize", "Genome Size:", value = "", width = NULL, placeholder = NULL),
                                                                    div(h5(HTML(paste("<b>Do you want to pass additional arguments to Canu other than genome size (beta, not recommended)? </b>Please read the docs", tags$a(href="https://canu.readthedocs.io/en/latest/parameter-reference.html", "here."))))),
                                                                    radioButtons(
                                                                      "passthrough", "Select an option", c("No", "Yes"), inline = TRUE, selected="No"
                                                                    ),
                                                                    conditionalPanel(condition = "input.passthrough =='No'",),
                                                                    conditionalPanel(condition = "input.passthrough =='Yes'",
                                                                                     div(h5(HTML("<b>Do not pass 'canu', -p <assembly-prefix>, -d <assembly-directory>, genomeSize=<number>[g|m|k], [-pacbio|-nanopore|-pacbio-hifi], or path to fastq. </b>"))),
                                                                                     textInput("canuargs", "Additional arguments:", value = "", width = NULL, placeholder = NULL)
                                                                    ),
                                                   ),
                                                   ),

                                  div( h5(HTML("<b>Do you want to use unaligened short read total RNA-seq data to enhance gene predictions?</b>"))),
                                  radioButtons(
                                    "RNA", "Select an option", c("Yes", "No"), inline = TRUE, selected="No"
                                  ),
                                  conditionalPanel(condition = "input.RNA =='No'",),
                                  conditionalPanel(condition="input.RNA == 'Yes'",
                                                   div( h5(HTML("<b>If RNA-seq sequencing was paired please select a path for both Read 1 and Read 2 (R1/R2) otherwise only select a path for R1. PERCEPTIVE will not use gzipped fastq!</b>"))),
                                                                                   suppressWarnings (shinyFilesButton('fastqR1RNA', 'Select pathway to unaligned RNA FASTQ Read 1' , 'Select pathway to unaligned RNA fastq Read 1:', multiple = FALSE,
                                                                                   buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder"))),
                                                                                   suppressWarnings (shinyFilesButton('fastqR2RNA', 'Select pathway to unaligned RNA FASTQ Read 2' , 'Select pathway to unaligned RNA fastq Read 2:', multiple = FALSE,
                                                                                   buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder")))),

                                  div( h5(HTML("<b>Name a well studied species that is within the same clade or that is related to this organism. </b> </br>For example, arabadopsis, cerevisae, chlamydomonas, human, or mouse. RepeatMasker will use this input to mask repeat sequences which might cause downstream spurious annotation. </b>"))),
                                  textInput("clade", "Clade:", value = "", width = NULL, placeholder = NULL),

                                  div( h6(HTML(" </br>Note: Interspursed repeats are mostly transposable elements in different states of erosion.
Thus, depending on evolutionary divergence since the origination of a transposable element, interspersed repeats are generally conserved in a clade of species.

</br> In principal, all unique clade names occurring in this database (http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html) can be used. Additional rules and examples: all inputs longer than one word must be in quotes. Input can range in specificity, for example <i> 'chlamydomonas reinhardtii'</i>, chimpanzee, fungi, canidae, mammals, etc. Picking a more specific option is prefered. "))),

                                  div( h5(HTML("</br> Please select an e-value cutoff for BLAST results. For example, 0.05 will result in a false discovery rate of 5%. Lower values will yield more stringent results, but 0.05 is recommended to start."))),
                                  textInput("eval", "e-value:", value = "0.05", width = NULL, placeholder = NULL),
                                  div( h5(HTML("Do you want to override PERCEPTIVE defaults and blast against your own non-model blast database? Please generate this database using instructions on github."))),
                                  radioButtons(
                                    "uniqueblast", "Select an option", c("No", "Yes"), inline = TRUE, selected="No"
                                  ),
                                  conditionalPanel(condition = "input.uniqueblast =='No'",),
                                  conditionalPanel(condition = "input.uniqueblast =='Yes'",
                                                   suppressWarnings (shinyFilesButton('blastdb', 'Select pathway to blastDB' , 'Select pathway to blastDB:', multiple = FALSE,
                                                                                      buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder"))),
                                  ),
                                  div(style = "display:inline-block; float:right; padding:10px",  actionButton("run", "Run PERCEPTIVE Pipeline", icon("gear"))),
                                  div(style = "display:inline-block; float:right", h6(HTML(" </br>Note: On a machine with 36 cores and 64GB of memory, PERCEPTIVE runs for ~24hrs for a 300MB genome.</b>"))),





                          ),



                          tabItem(tabName = "settings",
                                  div( h5(HTML("<b>Please identify the path to the singularity container (Perceptivev0.1.sif) downloaded from gdrive: </b>"))),
                                  suppressWarnings (shinyFilesButton('perceptive', 'Select Perceptivev0.1.sif' , 'Select pathway Perceptivev0.1.sif:', multiple = FALSE,
                                                   buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder"))),
                                  div( h5(HTML("<b>Please identify the path to the singularity container (braker3.sif) downloaded from gdrive: </b>"))),
                                  suppressWarnings (shinyFilesButton('braker', 'Select braker3.sif' , 'Select pathway to braker3.sif:', multiple = FALSE,
                                                   buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder"))),
                                  div( h5(HTML("<b>Please identify the path to the Eukaryota.fa OrthoDB database downloaded from gdrive: </b>"))),
                                  suppressWarnings (shinyFilesButton('orthodb', 'Select Eukaryota.fa' , 'Select pathway to Eukaryota.fa:', multiple = FALSE,
                                                   buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder"))),
                                div(h5(HTML(paste("<b> Please identify the path to the license for GeneMark which can be obtained", tags$a(href="http://topaz.gatech.edu/GeneMark/license_download.cgi", "here."), " License file should be unzipped using command gzip -d</b>")))),
                                suppressWarnings (shinyFilesButton('genemark', 'Select GeneMark License' , 'Select pathway to gmeskey:', multiple = FALSE,
                                                 buttonType = "default", class = NULL, style="color: #fff; background-color: #337ab7; border-color: #2e6da4", icon("folder"))),
                                div( h6(HTML("GeneMark Licenses expire. If PERCEPTIVE fails to complete, verify that your license is up to date. "))),
                                div(style = "display:inline-block; float:right; padding:10px",  actionButton("update", "Set Environmental Variables", icon("handshake"))),



                          )
                        )
                      ),

                      footer = dashboardFooter(left = "COPYRIGHT HOLDER: © 2024. Triad National Security, LLC. All rights reserved.", right = "LA-UR-24-21779"),
                      controlbar = dashboardControlbar(collapsed = TRUE, overlay= FALSE, skin="dark", skinSelector())


  )




  shinyApp(ui, server)
}
