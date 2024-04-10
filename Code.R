######....................................................................
#.........................................................................
#Autor: Tomás Texis
#Título: Barras_apiladas_microbiota
#Fecha: 6/Noviembre/23
#Objetivo: Ocuparemos barras apiladas para describir la composición de la
#microbiota intestinal de los pacientes
#.........................................................................
######....................................................................

#.....................................................................................
#.....................................Paquetes........................................
#.........................................................

library(ggthemes)
library(tidyverse)
library(ggplot2)
library(psych)
library(ggpubr) #para la funcion stat_compare_means
library(ggtext) #titulos en italicas
library(scales) #escalas
library(forcats) #Change order within factors ggplot
library(RColorBrewer)
library(ggdist)
getwd()
setwd("C://Users//hp i5 8va//Documents//farmacogenomica//Graficos")



#Leemos el archivos Tabla de asvs unida de los 3 lotes
Taxones_merged_corrected_w_batch <-
  read.table("C:\\Users\\hp i5 8va\\Documents\\farmacogenomica\\Genus_merged_batches_ConQuR_correction (1).txt", 
             header = TRUE) 

# Reorganizamos datos para plot a lo largo de las columnas 2 a 228 en dos 
# columnas: "Bacteria" (que contiene los nombres de las columnas originales) y 
#"Freq" (que contiene los valores)
#Reordenamos para el complot de barras
#Reordenamos para el plot de barras
Taxones_merged_corrected_w_batch_gather <-
  gather(Taxones_merged_corrected_w_batch,
         2:228, key = Bacteria, value = Freq)

#.....................................................................................
#####Separar taxones#####
#.....................................................................................
#separamos en base a la columna "bacteria"
Taxones_merged_corrected_w_batch_split <-
  separate(Taxones_merged_corrected_w_batch_gather,
           Bacteria,
           c(NA,
             "Kingdom",
             "Phylum",
             "Clase",
             "Order",
             "Family",
             "Genus",
             "Species"),
           sep = "[dpcofgs]__",
           fill = "right")
#le quitamos el punto a el texto de la columan
Taxones_merged_corrected_w_batch_split <-
  lapply(Taxones_merged_corrected_w_batch_split,
         gsub,
         pattern='\\.', replacement='')%>%
  as.data.frame()
#Con esta linea de codigo nos percatamso que existen Firmicutes________ y Firmicutes
#como identificadores diferneste lo cual no es correcto
table(Taxones_merged_corrected_w_batch_split$Phylum)
table(Taxones_merged_corrected_w_batch_split$Genus)

#####Homogeinizar taxones#####
#.....................................................................................
#Se cambia cualquier especie en formato "uncultured" por NA
#igual unidentified, gut_metagenome, human_gut, metagenome,
#wallaby_gut
Taxones_merged_corrected_w_batch_homogenizado <-
  Taxones_merged_corrected_w_batch_split

Taxones_merged_corrected_w_batch_homogenizado$Phylum <-
  gsub(Taxones_merged_corrected_w_batch_homogenizado$Phylum,
       pattern = "Firmicutes ",
       replacement = "Firmicutes")

Taxones_merged_corrected_w_batch_homogenizado$Phylum <-
  gsub(Taxones_merged_corrected_w_batch_homogenizado$Phylum,
       pattern = "Firmicutes________",
       replacement = "Firmicutes")

Taxones_merged_corrected_w_batch_homogenizado$Genus <-
  gsub(Taxones_merged_corrected_w_batch_homogenizado$Genus,
       pattern = "Bacteroides ",
       replacement = "Bacteroides")

Taxones_merged_corrected_w_batch_homogenizado$Genus <-
  gsub(Taxones_merged_corrected_w_batch_homogenizado$Genus,
       pattern = "uncultured ",
       replacement = "uncultured")

Taxones_merged_corrected_w_batch_homogenizado$Genus <-
  gsub(Taxones_merged_corrected_w_batch_homogenizado$Genus,
       pattern = "uncultured",
       replacement = "Uncultured")

Taxones_merged_corrected_w_batch_homogenizado[is.na(Taxones_merged_corrected_w_batch_homogenizado)] <- "No Asignado"


#.....................................................................................
#####Descriptiva taxones#####
#.....................................................................................
Taxones_merged_corrected_w_batch_homogenizado$Freq <-
  as.integer(Taxones_merged_corrected_w_batch_homogenizado$Freq)

#Debido a que ya tenemos reportada la Frecuencia, solo es cuestion
#de sumarla de acuerdo a su correspondeinte taxon
Descriptiva_Genera <-
  aggregate(Taxones_merged_corrected_w_batch_homogenizado$Freq,
            by=list(Genera=Taxones_merged_corrected_w_batch_homogenizado$Genus),
            FUN=sum)

Descriptiva_Genera <-
  mutate(Descriptiva_Genera,
         Freq = x,
         Freq_percent = x/sum(Descriptiva_Genera$x)*100) %>%
  select(-x)

#Lo mismo sucede con Phylum
Descriptiva_Phylum <-
  aggregate(Taxones_merged_corrected_w_batch_homogenizado$Freq,
            by=list(Phylum=Taxones_merged_corrected_w_batch_homogenizado$Phylum),
            FUN=sum)

Descriptiva_Phylum <-
  mutate(Descriptiva_Phylum,
         Freq = x,
         Freq_percent = x/sum(Descriptiva_Phylum$x)*100) %>%
  select(-x)

#Especie
Descriptiva_Especie <-
  aggregate(Taxones_merged_corrected_w_batch_homogenizado$Freq,
            by=list(Especie=Taxones_merged_corrected_w_batch_homogenizado$Species),
            FUN=sum)

Descriptiva_Especie <-
  mutate(Descriptiva_Especie,
         Freq = x,
         Freq_percent = x/sum(Descriptiva_Especie$x)*100) %>%
  select(-x)

#.....................................................................................
#####Plot phylum only#####
#.....................................................................................
Descriptiva_phylum_plot <-
  mutate(Descriptiva_Phylum,
         Base = "SILVA",
         Phylum_top = ifelse(Freq_percent > 0.5, Phylum, "Otros"),
         Phylum_top_format = str_replace(Phylum_top,
                                         "(.*acter.*a)", "*\\1*"),
         Phylum_top_format = str_replace(Phylum_top_format,
                                         "(.*tes)", "*\\1*"),
         Phylum_top_format = str_replace(Phylum_top_format,
                                         "(.*biota)", "*\\1*"),
         Phylum_top_format_order = fct_reorder(Phylum_top_format,
                                               Freq_percent,
                                               .desc = T))

ggplot(Descriptiva_phylum_plot,
       aes(x = Base,
           y = Freq_percent,
           fill = Phylum_top_format_order)) +
  geom_col(alpha = 0.9) +
  labs(y = "Frecuencia relativa, %") +
  scale_x_discrete(labels=c("SILVA" = "Phylum")) +
  scale_fill_manual(values = c(brewer.pal(3, "Set2"),
                               "#A6D854",
                               "#9F2B68",
                               "grey")) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.key.size = unit(10, "pt"))

ggsave(filename = "Descriptiva_phylum_merged_batches_corrected_all.png",
       plot = last_plot(),
       width = 3,
       height = 5,
       device='png',
       dpi=600)
#.....................................................................................
#####Plot genus only#####
#.....................................................................................
Descriptiva_genera_plot <-
  mutate(Descriptiva_Genera,
         Base = "SILVA",
         Genera_top = ifelse(Freq_percent > 3.5, Genera, "Otros"),
         Genera_top = ifelse(Genera == "No Asignado",
                             "No Asignado", Genera_top),
         Genera_top_format = str_replace(Genera_top,
                                         "(.*es.*)", "*\\1*"),
         Genera_top_format = str_replace(Genera_top_format,
                                         "(.*bacter.*)", "*\\1*"),
         Genera_top_format = str_replace(Genera_top_format,
                                         "(.*Esch.*)", "*\\1*"),
         Genera_top_format = str_replace(Genera_top_format,
                                         "(.*_.*)", "*\\1*"),
         Genera_top_format_order = fct_reorder(Genera_top_format,
                                               Freq_percent,
                                               .desc = T))

Descriptiva_genera_plot$Genera_top_format_order <-
  factor(Descriptiva_genera_plot$Genera_top_format_order,
         levels = c("*Bacteroides*","*Prevotella_9*",
                    "*Faecalibacterium*", "*Lachnospiraceae_NK4A136_group*",
                    "*Alistipes*","*EscherichiaShigella*",
                    "Otros","No Asignado"))

ggplot(Descriptiva_genera_plot, 
       aes(x = Base,
           y = Freq_percent,
           fill = Genera_top_format_order)) +
  geom_col(alpha = 0.9) +
  labs(y = "Frecuencia relativa, %") +
  scale_x_discrete(labels=c("SILVA" = "Genero")) +
  theme_classic() +
  scale_fill_manual(values = c("#66C2A5",
                               "#FC8D62",
                               "#6a83bb",
                               "#E78AC3",
                               "#A6D854",
                               "red",
                               "grey",
                               "black")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.key.size = unit(10, "pt"))

ggsave(filename = "Descriptiva_genera_merged_batches_corrected_ALL.png",
       plot = last_plot(),
       width = 4,
       height = 5,
       device='png',
       dpi=600)

#.....................................................................................
#####Plot individuos por Phylum#####
#.....................................................................................
#Arreglar datos para plot
Descriptiva_phylum_plot_ID<-
  Taxones_merged_corrected_w_batch_homogenizado %>%
  select(ID,
         Phylum,
         Freq) %>%
  group_by(ID,
           Phylum) %>%
  summarise(Freq = sum(Freq)) %>%
  mutate(Freq_rel = Freq/sum(Freq)) %>%
  ungroup() %>%
  mutate(Phylum_top = ifelse(Phylum == "Firmicutes", "*Firmicutes*",
                             ifelse(Phylum == "Bacteroidota", "*Bacteroidota*",
                                    ifelse(Phylum == "Proteobacteria", "*Proteobacteria*",
                                           ifelse(Phylum == "Verrucomicrobiota", "*Verrucomicrobiota*",
                                                  ifelse(Phylum == "Actinobacteriota", "*Actinobacteriota*",
                                                         "Others"))))))

Tablas_phylum_ID_Freq_rel <-
  select(Descriptiva_phylum_plot_ID,
         ID,
         Freq_rel,
         Phylum_top) %>%
  group_by(ID, Phylum_top) %>%
  summarize(Freq_rel = round(sum(Freq_rel),2)) %>% #Aqui agrupamos las bacterias "others"
  ungroup() %>%
  spread(key = Phylum_top,
         value = Freq_rel)
#Guardamos la base en una tabla de txt esto fue por que lo solicito la Dra 
write.table(Tablas_phylum_ID_Freq_rel,
            file = "Tablas_phylum_ID_Freq_rel.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

#Crear una base que ordene segun la freq rel de Firmicutes
Descriptiva_phylum_plot_ID_order <-
  filter(Descriptiva_phylum_plot_ID,
         Phylum == "Firmicutes") %>%
  arrange(desc(Freq_rel)) %>%
  ungroup() %>%
  mutate(order = 1:nrow(.)) %>%
  select(ID, order)

#Unimos las bases que contienen el ID ordenado y los datos
Descriptiva_phylum_plot_ID_full <-
  full_join(x = Descriptiva_phylum_plot_ID_order,
            y = Descriptiva_phylum_plot_ID,
            by = "ID")

#Ordenar factores de phylum
Descriptiva_phylum_plot_ID_full$Phylum_top <-
  factor(Descriptiva_phylum_plot_ID_full$Phylum_top,
         levels = c("*Firmicutes*", "*Bacteroidota*",
                    "*Proteobacteria*", "*Verrucomicrobiota*",
                    "*Actinobacteriota*", "Others"))

#Plot
Descripcion_Phylum_all_sample <-
  ggplot(Descriptiva_phylum_plot_ID_full,
       aes(x = as.factor(order),
           y = Freq_rel,
           fill = Phylum_top)) +
  geom_col(width = 0.9) +
  labs(y = "Relative abundance",
       x = "Samples") +
  theme_classic() +
  scale_fill_manual(values = c(brewer.pal(3, "Dark2"),
                               "#A6D854",
                               "#9F2B68",
                               "grey")) +
  scale_x_discrete(breaks = Descriptiva_phylum_plot_ID_full$order) +
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom")


ggsave(filename = "Abundanciadephylum.png",
       plot = last_plot(),
       width = 10,
       height = 6,
       device='png',
       dpi=600)


#.....................................................................................
#####Plot de genero por individuo#####
#.....................................................................................
#Arreglar datos para plot

Descriptiva_Genus_plot_ID<-
  Taxones_merged_corrected_w_batch_homogenizado %>%
  select(ID,
         Genus,
         Freq) %>%
  group_by(ID,
           Genus) %>%
  summarise(Freq = sum(Freq)) %>%
  mutate(Freq_rel = Freq/sum(Freq)) %>%
  ungroup() %>%
  mutate(Genus_top = ifelse(Genus == "Bacteroides", "*Bacteroides*",
                             ifelse(Genus == "Prevotella_9", "*Prevotella_9*",
                                    ifelse(Genus == "Faecalibacterium", "*Faecalibacterium*",
                                           ifelse(Genus == "Lachnospiraceae_NK4A136_group", "*Lachnospiraceae_NK4A136_group*",
                                                  ifelse(Genus == "Alistipes", "*Alistipes*",
                                                         ifelse(Genus == "EscherichiaShigella", "*EscherichiaShigella*",
                                                                ifelse(Genus == "Eubacterium_xylanophilum_group", "*Eubacterium_xylanophilum_group*",
                                                                          #ifelse(Genus == "UCG002", "*UCG002*",
                                                                       ifelse(Genus == "No Asignado",
                                                                       "No Asignado",
                                                         "Others")))))))))


##### con estos vimos que generso tenia un %5 maypr por muestra y en cuantas muestras aparecian 
### para saber cuales utilizamos para poner en el plot

Generosmayoracinco <- Descriptiva_Genus_plot_ID %>%
  group_by(Genus)%>%
  filter(Freq_rel > 5)
table(Generosmayoracinco$Genus)

Tablas_Genus_ID_Freq_rel <-
  select(Descriptiva_Genus_plot_ID,
         ID,
         Freq_rel,
         Genus_top) %>%
  group_by(ID, Genus_top) %>%
  summarize(Freq_rel = round(sum(Freq_rel),2)) %>% #Aqui agrupamos las bacterias "others"
  ungroup() %>%
  spread(key = Genus_top,
         value = Freq_rel)

#Crear una base que ordene segun la freq rel de Bacteroides"
Descriptiva_Genus_plot_ID_order <-
  filter(Descriptiva_Genus_plot_ID,
         Genus == "Bacteroides") %>%
  arrange(desc(Freq_rel)) %>%
  ungroup() %>%
  mutate(order = 1:nrow(.)) %>%
  select(ID, order)

#Unimos las bases que contienen el ID ordenado y los datos
Descriptiva_Genus_plot_ID_full <-
  full_join(x = Descriptiva_Genus_plot_ID_order,
            y = Descriptiva_Genus_plot_ID,
            by = "ID")%>%
  mutate(ID = as.double(ID)) #volvemos la columa ID numerica

#Ordenar factores de phylum
Descriptiva_Genus_plot_ID_full$Genus_top <-
  factor(Descriptiva_Genus_plot_ID_full$Genus_top,
         levels = c("*Bacteroides*","*Prevotella_9*",
                    "*Faecalibacterium*", "*Lachnospiraceae_NK4A136_group*",
                    "*Alistipes*", "*EscherichiaShigella*","*UCG002*",
                    "*Eubacterium_xylanophilum_group*","Others","No Asignado"))

sigloxxl <- vroom::vroom("C:/Users/hp i5 8va/Dropbox/Farmacogenomica microbioma/Microbiota adultos mayore/Base_clinica_para_microbioma_257samples.txt")%>%
  filter(ID %in% c(20035, 25733,20035,569, 174, 26004,20077,236,820,458,494,372,515,
                   25485,2800,21158,20708,22610,
                   435,26042,795,21431,871,888,871,49))
siglo21 <- sigloxxl%>%
  select(ID,EDAD,PS,PD)%>%
  mutate( rangosedad =ifelse( EDAD <= 69, "Adultos Mayores 60-69 Años",
                              ifelse( EDAD >= 70 & EDAD <= 79, "Adultos Mayores 70-79 Años",
                                      ifelse(EDAD >= 80, "Adultos Mayores ≥80 Años",NA_character_)))) %>%
  mutate(ID = as.character(ID))

union <-
  inner_join(siglo21, Descriptiva_Genus_plot_ID, by = "ID")%>%
  group_by(ID) %>% 
  filter(!(any(Genus == "Bacteroides" & Freq_rel > 0.05)))%>%
  filter(!(any(Genus == "Prevotella_9" & Freq_rel > 0.05)))%>%
  filter(!(ID %in% c(21431, 871, 888, 871, 49))) %>%
  group_by(ID)
  #filter(ID %in% c(20035, 25733,20035,569, 174, 26004,20077,236,458,515,
   #                 25485,2800,21158,20708,
    #                435,26042,795,21431, 26004,558))%>%
filter(!(any(ID= 20035, 25733,20035,569, 174, 26004,20077,236,820,458,494,372,515,
             25485,2800,21158,20708,22610,
             435,26042,795,21431,871,888,871,49)))%>%

    select(ID,
         Genus,
         Freq) %>%
  group_by(Genus) %>%
  summarise(Freq = sum(Freq)) %>%
  mutate(Freq_rel = Freq/sum(Freq)) %>%
  ungroup()

### Union de las 2 bases 
#Descriptiva_Genus_plot_ID_full <- left_join(Descriptiva_Genus_plot_ID_full, siglo21, by = "ID") %>%
 # filter(rangosedad == "≥80")

#Plot
Descripcion_genero_all_sample <- 
ggplot(Descriptiva_Genus_plot_ID_full,
       aes(x = as.factor(order),
           y = Freq_rel,
           fill = Genus_top)) +
  geom_col(width = 0.9) +
  labs(y = "Relative abundance, %",
       x = "Samples (N=241)") +
  theme_classic() +
  scale_fill_manual(values = c( brewer.pal(5,"Set2"),
                                "red",
                                "#9F2B68",
                                "grey",
                                "black")) +
  scale_x_discrete(breaks = Descriptiva_Genus_plot_ID_full$order) +
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.key.size = unit(15, "pt"),
        legend.position = "bottom")


ggsave(filename = "AbundanciagencontodosgenersoN.png",
       plot = Descripcion_genero_all_sample,
       width = 10,
       height = 6,
       device='png',
       dpi=600)
