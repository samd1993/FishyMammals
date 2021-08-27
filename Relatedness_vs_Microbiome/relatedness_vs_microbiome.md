
## Generate fish tree–&gt;make distance matrix of it–&gt;then melt matrix for linear model

``` r
library(fishtree)
library(reshape2)
genomic_fish <- c("Epibulus insidiator","Naso lituratus","Chlorurus spilurus","Ctenochaetus striatus","Cephalopholis argus","Cephalopholis urodeta","Epinephelus merra","Halichoeres trimaculatus","Dascyllus flavicaudus","Aulostomus chinensis","Abudefduf sexfasciatus","Chaetodon ornatissimus","Chaetodon reticulatus","Chaetodon auriga","Rhinecanthus aculeatus","Zebrasoma scopas","Sargocentron spiniferum","Myripristis violacea", "Stegastes nigricans","Myripristis berndti")

tree <- fishtree::fishtree_phylogeny(specie=genomic_fish)

library(ape)
write.tree(tree, file="fish-tree.nwk")

relatedness_matrix_fish <- ape::vcv(tree, corr = TRUE)
nrow(relatedness_matrix_fish)
```

    ## [1] 20

``` r
rmf=melt(relatedness_matrix_fish)
rmf
```

    ##                         Var1                     Var2     value
    ## 1      Dascyllus_flavicaudus    Dascyllus_flavicaudus 1.0000000
    ## 2     Abudefduf_sexfasciatus    Dascyllus_flavicaudus 0.6523906
    ## 3        Stegastes_nigricans    Dascyllus_flavicaudus 0.6242110
    ## 4       Aulostomus_chinensis    Dascyllus_flavicaudus 0.1176839
    ## 5          Epinephelus_merra    Dascyllus_flavicaudus 0.1046071
    ## 6      Cephalopholis_urodeta    Dascyllus_flavicaudus 0.1046071
    ## 7        Cephalopholis_argus    Dascyllus_flavicaudus 0.1046071
    ## 8           Chaetodon_auriga    Dascyllus_flavicaudus 0.1046071
    ## 9      Chaetodon_reticulatus    Dascyllus_flavicaudus 0.1046071
    ## 10    Chaetodon_ornatissimus    Dascyllus_flavicaudus 0.1046071
    ## 11     Ctenochaetus_striatus    Dascyllus_flavicaudus 0.1046071
    ## 12          Zebrasoma_scopas    Dascyllus_flavicaudus 0.1046071
    ## 13            Naso_lituratus    Dascyllus_flavicaudus 0.1046071
    ## 14    Rhinecanthus_aculeatus    Dascyllus_flavicaudus 0.1046071
    ## 15  Halichoeres_trimaculatus    Dascyllus_flavicaudus 0.1046071
    ## 16        Chlorurus_spilurus    Dascyllus_flavicaudus 0.1046071
    ## 17       Epibulus_insidiator    Dascyllus_flavicaudus 0.1046071
    ## 18   Sargocentron_spiniferum    Dascyllus_flavicaudus 0.0000000
    ## 19       Myripristis_berndti    Dascyllus_flavicaudus 0.0000000
    ## 20      Myripristis_violacea    Dascyllus_flavicaudus 0.0000000
    ## 21     Dascyllus_flavicaudus   Abudefduf_sexfasciatus 0.6523906
    ## 22    Abudefduf_sexfasciatus   Abudefduf_sexfasciatus 1.0000000
    ## 23       Stegastes_nigricans   Abudefduf_sexfasciatus 0.6242110
    ## 24      Aulostomus_chinensis   Abudefduf_sexfasciatus 0.1176839
    ## 25         Epinephelus_merra   Abudefduf_sexfasciatus 0.1046071
    ## 26     Cephalopholis_urodeta   Abudefduf_sexfasciatus 0.1046071
    ## 27       Cephalopholis_argus   Abudefduf_sexfasciatus 0.1046071
    ## 28          Chaetodon_auriga   Abudefduf_sexfasciatus 0.1046071
    ## 29     Chaetodon_reticulatus   Abudefduf_sexfasciatus 0.1046071
    ## 30    Chaetodon_ornatissimus   Abudefduf_sexfasciatus 0.1046071
    ## 31     Ctenochaetus_striatus   Abudefduf_sexfasciatus 0.1046071
    ## 32          Zebrasoma_scopas   Abudefduf_sexfasciatus 0.1046071
    ## 33            Naso_lituratus   Abudefduf_sexfasciatus 0.1046071
    ## 34    Rhinecanthus_aculeatus   Abudefduf_sexfasciatus 0.1046071
    ## 35  Halichoeres_trimaculatus   Abudefduf_sexfasciatus 0.1046071
    ## 36        Chlorurus_spilurus   Abudefduf_sexfasciatus 0.1046071
    ## 37       Epibulus_insidiator   Abudefduf_sexfasciatus 0.1046071
    ## 38   Sargocentron_spiniferum   Abudefduf_sexfasciatus 0.0000000
    ## 39       Myripristis_berndti   Abudefduf_sexfasciatus 0.0000000
    ## 40      Myripristis_violacea   Abudefduf_sexfasciatus 0.0000000
    ## 41     Dascyllus_flavicaudus      Stegastes_nigricans 0.6242110
    ## 42    Abudefduf_sexfasciatus      Stegastes_nigricans 0.6242110
    ## 43       Stegastes_nigricans      Stegastes_nigricans 1.0000000
    ## 44      Aulostomus_chinensis      Stegastes_nigricans 0.1176839
    ## 45         Epinephelus_merra      Stegastes_nigricans 0.1046071
    ## 46     Cephalopholis_urodeta      Stegastes_nigricans 0.1046071
    ## 47       Cephalopholis_argus      Stegastes_nigricans 0.1046071
    ## 48          Chaetodon_auriga      Stegastes_nigricans 0.1046071
    ## 49     Chaetodon_reticulatus      Stegastes_nigricans 0.1046071
    ## 50    Chaetodon_ornatissimus      Stegastes_nigricans 0.1046071
    ## 51     Ctenochaetus_striatus      Stegastes_nigricans 0.1046071
    ## 52          Zebrasoma_scopas      Stegastes_nigricans 0.1046071
    ## 53            Naso_lituratus      Stegastes_nigricans 0.1046071
    ## 54    Rhinecanthus_aculeatus      Stegastes_nigricans 0.1046071
    ## 55  Halichoeres_trimaculatus      Stegastes_nigricans 0.1046071
    ## 56        Chlorurus_spilurus      Stegastes_nigricans 0.1046071
    ## 57       Epibulus_insidiator      Stegastes_nigricans 0.1046071
    ## 58   Sargocentron_spiniferum      Stegastes_nigricans 0.0000000
    ## 59       Myripristis_berndti      Stegastes_nigricans 0.0000000
    ## 60      Myripristis_violacea      Stegastes_nigricans 0.0000000
    ## 61     Dascyllus_flavicaudus     Aulostomus_chinensis 0.1176839
    ## 62    Abudefduf_sexfasciatus     Aulostomus_chinensis 0.1176839
    ## 63       Stegastes_nigricans     Aulostomus_chinensis 0.1176839
    ## 64      Aulostomus_chinensis     Aulostomus_chinensis 1.0000000
    ## 65         Epinephelus_merra     Aulostomus_chinensis 0.1046071
    ## 66     Cephalopholis_urodeta     Aulostomus_chinensis 0.1046071
    ## 67       Cephalopholis_argus     Aulostomus_chinensis 0.1046071
    ## 68          Chaetodon_auriga     Aulostomus_chinensis 0.1046071
    ## 69     Chaetodon_reticulatus     Aulostomus_chinensis 0.1046071
    ## 70    Chaetodon_ornatissimus     Aulostomus_chinensis 0.1046071
    ## 71     Ctenochaetus_striatus     Aulostomus_chinensis 0.1046071
    ## 72          Zebrasoma_scopas     Aulostomus_chinensis 0.1046071
    ## 73            Naso_lituratus     Aulostomus_chinensis 0.1046071
    ## 74    Rhinecanthus_aculeatus     Aulostomus_chinensis 0.1046071
    ## 75  Halichoeres_trimaculatus     Aulostomus_chinensis 0.1046071
    ## 76        Chlorurus_spilurus     Aulostomus_chinensis 0.1046071
    ## 77       Epibulus_insidiator     Aulostomus_chinensis 0.1046071
    ## 78   Sargocentron_spiniferum     Aulostomus_chinensis 0.0000000
    ## 79       Myripristis_berndti     Aulostomus_chinensis 0.0000000
    ## 80      Myripristis_violacea     Aulostomus_chinensis 0.0000000
    ## 81     Dascyllus_flavicaudus        Epinephelus_merra 0.1046071
    ## 82    Abudefduf_sexfasciatus        Epinephelus_merra 0.1046071
    ## 83       Stegastes_nigricans        Epinephelus_merra 0.1046071
    ## 84      Aulostomus_chinensis        Epinephelus_merra 0.1046071
    ## 85         Epinephelus_merra        Epinephelus_merra 1.0000000
    ## 86     Cephalopholis_urodeta        Epinephelus_merra 0.6638028
    ## 87       Cephalopholis_argus        Epinephelus_merra 0.6638028
    ## 88          Chaetodon_auriga        Epinephelus_merra 0.2569510
    ## 89     Chaetodon_reticulatus        Epinephelus_merra 0.2569510
    ## 90    Chaetodon_ornatissimus        Epinephelus_merra 0.2569510
    ## 91     Ctenochaetus_striatus        Epinephelus_merra 0.2046493
    ## 92          Zebrasoma_scopas        Epinephelus_merra 0.2046493
    ## 93            Naso_lituratus        Epinephelus_merra 0.2046493
    ## 94    Rhinecanthus_aculeatus        Epinephelus_merra 0.1802150
    ## 95  Halichoeres_trimaculatus        Epinephelus_merra 0.1296265
    ## 96        Chlorurus_spilurus        Epinephelus_merra 0.1296265
    ## 97       Epibulus_insidiator        Epinephelus_merra 0.1296265
    ## 98   Sargocentron_spiniferum        Epinephelus_merra 0.0000000
    ## 99       Myripristis_berndti        Epinephelus_merra 0.0000000
    ## 100     Myripristis_violacea        Epinephelus_merra 0.0000000
    ## 101    Dascyllus_flavicaudus    Cephalopholis_urodeta 0.1046071
    ## 102   Abudefduf_sexfasciatus    Cephalopholis_urodeta 0.1046071
    ## 103      Stegastes_nigricans    Cephalopholis_urodeta 0.1046071
    ## 104     Aulostomus_chinensis    Cephalopholis_urodeta 0.1046071
    ## 105        Epinephelus_merra    Cephalopholis_urodeta 0.6638028
    ## 106    Cephalopholis_urodeta    Cephalopholis_urodeta 1.0000000
    ## 107      Cephalopholis_argus    Cephalopholis_urodeta 0.7251554
    ## 108         Chaetodon_auriga    Cephalopholis_urodeta 0.2569510
    ## 109    Chaetodon_reticulatus    Cephalopholis_urodeta 0.2569510
    ## 110   Chaetodon_ornatissimus    Cephalopholis_urodeta 0.2569510
    ## 111    Ctenochaetus_striatus    Cephalopholis_urodeta 0.2046493
    ## 112         Zebrasoma_scopas    Cephalopholis_urodeta 0.2046493
    ## 113           Naso_lituratus    Cephalopholis_urodeta 0.2046493
    ## 114   Rhinecanthus_aculeatus    Cephalopholis_urodeta 0.1802150
    ## 115 Halichoeres_trimaculatus    Cephalopholis_urodeta 0.1296265
    ## 116       Chlorurus_spilurus    Cephalopholis_urodeta 0.1296265
    ## 117      Epibulus_insidiator    Cephalopholis_urodeta 0.1296265
    ## 118  Sargocentron_spiniferum    Cephalopholis_urodeta 0.0000000
    ## 119      Myripristis_berndti    Cephalopholis_urodeta 0.0000000
    ## 120     Myripristis_violacea    Cephalopholis_urodeta 0.0000000
    ## 121    Dascyllus_flavicaudus      Cephalopholis_argus 0.1046071
    ## 122   Abudefduf_sexfasciatus      Cephalopholis_argus 0.1046071
    ## 123      Stegastes_nigricans      Cephalopholis_argus 0.1046071
    ## 124     Aulostomus_chinensis      Cephalopholis_argus 0.1046071
    ## 125        Epinephelus_merra      Cephalopholis_argus 0.6638028
    ## 126    Cephalopholis_urodeta      Cephalopholis_argus 0.7251554
    ## 127      Cephalopholis_argus      Cephalopholis_argus 1.0000000
    ## 128         Chaetodon_auriga      Cephalopholis_argus 0.2569510
    ## 129    Chaetodon_reticulatus      Cephalopholis_argus 0.2569510
    ## 130   Chaetodon_ornatissimus      Cephalopholis_argus 0.2569510
    ## 131    Ctenochaetus_striatus      Cephalopholis_argus 0.2046493
    ## 132         Zebrasoma_scopas      Cephalopholis_argus 0.2046493
    ## 133           Naso_lituratus      Cephalopholis_argus 0.2046493
    ## 134   Rhinecanthus_aculeatus      Cephalopholis_argus 0.1802150
    ## 135 Halichoeres_trimaculatus      Cephalopholis_argus 0.1296265
    ## 136       Chlorurus_spilurus      Cephalopholis_argus 0.1296265
    ## 137      Epibulus_insidiator      Cephalopholis_argus 0.1296265
    ## 138  Sargocentron_spiniferum      Cephalopholis_argus 0.0000000
    ## 139      Myripristis_berndti      Cephalopholis_argus 0.0000000
    ## 140     Myripristis_violacea      Cephalopholis_argus 0.0000000
    ## 141    Dascyllus_flavicaudus         Chaetodon_auriga 0.1046071
    ## 142   Abudefduf_sexfasciatus         Chaetodon_auriga 0.1046071
    ## 143      Stegastes_nigricans         Chaetodon_auriga 0.1046071
    ## 144     Aulostomus_chinensis         Chaetodon_auriga 0.1046071
    ## 145        Epinephelus_merra         Chaetodon_auriga 0.2569510
    ## 146    Cephalopholis_urodeta         Chaetodon_auriga 0.2569510
    ## 147      Cephalopholis_argus         Chaetodon_auriga 0.2569510
    ## 148         Chaetodon_auriga         Chaetodon_auriga 1.0000000
    ## 149    Chaetodon_reticulatus         Chaetodon_auriga 0.8472185
    ## 150   Chaetodon_ornatissimus         Chaetodon_auriga 0.8472185
    ## 151    Ctenochaetus_striatus         Chaetodon_auriga 0.2046493
    ## 152         Zebrasoma_scopas         Chaetodon_auriga 0.2046493
    ## 153           Naso_lituratus         Chaetodon_auriga 0.2046493
    ## 154   Rhinecanthus_aculeatus         Chaetodon_auriga 0.1802150
    ## 155 Halichoeres_trimaculatus         Chaetodon_auriga 0.1296265
    ## 156       Chlorurus_spilurus         Chaetodon_auriga 0.1296265
    ## 157      Epibulus_insidiator         Chaetodon_auriga 0.1296265
    ## 158  Sargocentron_spiniferum         Chaetodon_auriga 0.0000000
    ## 159      Myripristis_berndti         Chaetodon_auriga 0.0000000
    ## 160     Myripristis_violacea         Chaetodon_auriga 0.0000000
    ## 161    Dascyllus_flavicaudus    Chaetodon_reticulatus 0.1046071
    ## 162   Abudefduf_sexfasciatus    Chaetodon_reticulatus 0.1046071
    ## 163      Stegastes_nigricans    Chaetodon_reticulatus 0.1046071
    ## 164     Aulostomus_chinensis    Chaetodon_reticulatus 0.1046071
    ## 165        Epinephelus_merra    Chaetodon_reticulatus 0.2569510
    ## 166    Cephalopholis_urodeta    Chaetodon_reticulatus 0.2569510
    ## 167      Cephalopholis_argus    Chaetodon_reticulatus 0.2569510
    ## 168         Chaetodon_auriga    Chaetodon_reticulatus 0.8472185
    ## 169    Chaetodon_reticulatus    Chaetodon_reticulatus 1.0000000
    ## 170   Chaetodon_ornatissimus    Chaetodon_reticulatus 0.9668092
    ## 171    Ctenochaetus_striatus    Chaetodon_reticulatus 0.2046493
    ## 172         Zebrasoma_scopas    Chaetodon_reticulatus 0.2046493
    ## 173           Naso_lituratus    Chaetodon_reticulatus 0.2046493
    ## 174   Rhinecanthus_aculeatus    Chaetodon_reticulatus 0.1802150
    ## 175 Halichoeres_trimaculatus    Chaetodon_reticulatus 0.1296265
    ## 176       Chlorurus_spilurus    Chaetodon_reticulatus 0.1296265
    ## 177      Epibulus_insidiator    Chaetodon_reticulatus 0.1296265
    ## 178  Sargocentron_spiniferum    Chaetodon_reticulatus 0.0000000
    ## 179      Myripristis_berndti    Chaetodon_reticulatus 0.0000000
    ## 180     Myripristis_violacea    Chaetodon_reticulatus 0.0000000
    ## 181    Dascyllus_flavicaudus   Chaetodon_ornatissimus 0.1046071
    ## 182   Abudefduf_sexfasciatus   Chaetodon_ornatissimus 0.1046071
    ## 183      Stegastes_nigricans   Chaetodon_ornatissimus 0.1046071
    ## 184     Aulostomus_chinensis   Chaetodon_ornatissimus 0.1046071
    ## 185        Epinephelus_merra   Chaetodon_ornatissimus 0.2569510
    ## 186    Cephalopholis_urodeta   Chaetodon_ornatissimus 0.2569510
    ## 187      Cephalopholis_argus   Chaetodon_ornatissimus 0.2569510
    ## 188         Chaetodon_auriga   Chaetodon_ornatissimus 0.8472185
    ## 189    Chaetodon_reticulatus   Chaetodon_ornatissimus 0.9668092
    ## 190   Chaetodon_ornatissimus   Chaetodon_ornatissimus 1.0000000
    ## 191    Ctenochaetus_striatus   Chaetodon_ornatissimus 0.2046493
    ## 192         Zebrasoma_scopas   Chaetodon_ornatissimus 0.2046493
    ## 193           Naso_lituratus   Chaetodon_ornatissimus 0.2046493
    ## 194   Rhinecanthus_aculeatus   Chaetodon_ornatissimus 0.1802150
    ## 195 Halichoeres_trimaculatus   Chaetodon_ornatissimus 0.1296265
    ## 196       Chlorurus_spilurus   Chaetodon_ornatissimus 0.1296265
    ## 197      Epibulus_insidiator   Chaetodon_ornatissimus 0.1296265
    ## 198  Sargocentron_spiniferum   Chaetodon_ornatissimus 0.0000000
    ## 199      Myripristis_berndti   Chaetodon_ornatissimus 0.0000000
    ## 200     Myripristis_violacea   Chaetodon_ornatissimus 0.0000000
    ## 201    Dascyllus_flavicaudus    Ctenochaetus_striatus 0.1046071
    ## 202   Abudefduf_sexfasciatus    Ctenochaetus_striatus 0.1046071
    ## 203      Stegastes_nigricans    Ctenochaetus_striatus 0.1046071
    ## 204     Aulostomus_chinensis    Ctenochaetus_striatus 0.1046071
    ## 205        Epinephelus_merra    Ctenochaetus_striatus 0.2046493
    ## 206    Cephalopholis_urodeta    Ctenochaetus_striatus 0.2046493
    ## 207      Cephalopholis_argus    Ctenochaetus_striatus 0.2046493
    ## 208         Chaetodon_auriga    Ctenochaetus_striatus 0.2046493
    ## 209    Chaetodon_reticulatus    Ctenochaetus_striatus 0.2046493
    ## 210   Chaetodon_ornatissimus    Ctenochaetus_striatus 0.2046493
    ## 211    Ctenochaetus_striatus    Ctenochaetus_striatus 1.0000000
    ## 212         Zebrasoma_scopas    Ctenochaetus_striatus 0.6998904
    ## 213           Naso_lituratus    Ctenochaetus_striatus 0.5551201
    ## 214   Rhinecanthus_aculeatus    Ctenochaetus_striatus 0.1802150
    ## 215 Halichoeres_trimaculatus    Ctenochaetus_striatus 0.1296265
    ## 216       Chlorurus_spilurus    Ctenochaetus_striatus 0.1296265
    ## 217      Epibulus_insidiator    Ctenochaetus_striatus 0.1296265
    ## 218  Sargocentron_spiniferum    Ctenochaetus_striatus 0.0000000
    ## 219      Myripristis_berndti    Ctenochaetus_striatus 0.0000000
    ## 220     Myripristis_violacea    Ctenochaetus_striatus 0.0000000
    ## 221    Dascyllus_flavicaudus         Zebrasoma_scopas 0.1046071
    ## 222   Abudefduf_sexfasciatus         Zebrasoma_scopas 0.1046071
    ## 223      Stegastes_nigricans         Zebrasoma_scopas 0.1046071
    ## 224     Aulostomus_chinensis         Zebrasoma_scopas 0.1046071
    ## 225        Epinephelus_merra         Zebrasoma_scopas 0.2046493
    ## 226    Cephalopholis_urodeta         Zebrasoma_scopas 0.2046493
    ## 227      Cephalopholis_argus         Zebrasoma_scopas 0.2046493
    ## 228         Chaetodon_auriga         Zebrasoma_scopas 0.2046493
    ## 229    Chaetodon_reticulatus         Zebrasoma_scopas 0.2046493
    ## 230   Chaetodon_ornatissimus         Zebrasoma_scopas 0.2046493
    ## 231    Ctenochaetus_striatus         Zebrasoma_scopas 0.6998904
    ## 232         Zebrasoma_scopas         Zebrasoma_scopas 1.0000000
    ## 233           Naso_lituratus         Zebrasoma_scopas 0.5551201
    ## 234   Rhinecanthus_aculeatus         Zebrasoma_scopas 0.1802150
    ## 235 Halichoeres_trimaculatus         Zebrasoma_scopas 0.1296265
    ## 236       Chlorurus_spilurus         Zebrasoma_scopas 0.1296265
    ## 237      Epibulus_insidiator         Zebrasoma_scopas 0.1296265
    ## 238  Sargocentron_spiniferum         Zebrasoma_scopas 0.0000000
    ## 239      Myripristis_berndti         Zebrasoma_scopas 0.0000000
    ## 240     Myripristis_violacea         Zebrasoma_scopas 0.0000000
    ## 241    Dascyllus_flavicaudus           Naso_lituratus 0.1046071
    ## 242   Abudefduf_sexfasciatus           Naso_lituratus 0.1046071
    ## 243      Stegastes_nigricans           Naso_lituratus 0.1046071
    ## 244     Aulostomus_chinensis           Naso_lituratus 0.1046071
    ## 245        Epinephelus_merra           Naso_lituratus 0.2046493
    ## 246    Cephalopholis_urodeta           Naso_lituratus 0.2046493
    ## 247      Cephalopholis_argus           Naso_lituratus 0.2046493
    ## 248         Chaetodon_auriga           Naso_lituratus 0.2046493
    ## 249    Chaetodon_reticulatus           Naso_lituratus 0.2046493
    ## 250   Chaetodon_ornatissimus           Naso_lituratus 0.2046493
    ## 251    Ctenochaetus_striatus           Naso_lituratus 0.5551201
    ## 252         Zebrasoma_scopas           Naso_lituratus 0.5551201
    ## 253           Naso_lituratus           Naso_lituratus 1.0000000
    ## 254   Rhinecanthus_aculeatus           Naso_lituratus 0.1802150
    ## 255 Halichoeres_trimaculatus           Naso_lituratus 0.1296265
    ## 256       Chlorurus_spilurus           Naso_lituratus 0.1296265
    ## 257      Epibulus_insidiator           Naso_lituratus 0.1296265
    ## 258  Sargocentron_spiniferum           Naso_lituratus 0.0000000
    ## 259      Myripristis_berndti           Naso_lituratus 0.0000000
    ## 260     Myripristis_violacea           Naso_lituratus 0.0000000
    ## 261    Dascyllus_flavicaudus   Rhinecanthus_aculeatus 0.1046071
    ## 262   Abudefduf_sexfasciatus   Rhinecanthus_aculeatus 0.1046071
    ## 263      Stegastes_nigricans   Rhinecanthus_aculeatus 0.1046071
    ## 264     Aulostomus_chinensis   Rhinecanthus_aculeatus 0.1046071
    ## 265        Epinephelus_merra   Rhinecanthus_aculeatus 0.1802150
    ## 266    Cephalopholis_urodeta   Rhinecanthus_aculeatus 0.1802150
    ## 267      Cephalopholis_argus   Rhinecanthus_aculeatus 0.1802150
    ## 268         Chaetodon_auriga   Rhinecanthus_aculeatus 0.1802150
    ## 269    Chaetodon_reticulatus   Rhinecanthus_aculeatus 0.1802150
    ## 270   Chaetodon_ornatissimus   Rhinecanthus_aculeatus 0.1802150
    ## 271    Ctenochaetus_striatus   Rhinecanthus_aculeatus 0.1802150
    ## 272         Zebrasoma_scopas   Rhinecanthus_aculeatus 0.1802150
    ## 273           Naso_lituratus   Rhinecanthus_aculeatus 0.1802150
    ## 274   Rhinecanthus_aculeatus   Rhinecanthus_aculeatus 1.0000000
    ## 275 Halichoeres_trimaculatus   Rhinecanthus_aculeatus 0.1296265
    ## 276       Chlorurus_spilurus   Rhinecanthus_aculeatus 0.1296265
    ## 277      Epibulus_insidiator   Rhinecanthus_aculeatus 0.1296265
    ## 278  Sargocentron_spiniferum   Rhinecanthus_aculeatus 0.0000000
    ## 279      Myripristis_berndti   Rhinecanthus_aculeatus 0.0000000
    ## 280     Myripristis_violacea   Rhinecanthus_aculeatus 0.0000000
    ## 281    Dascyllus_flavicaudus Halichoeres_trimaculatus 0.1046071
    ## 282   Abudefduf_sexfasciatus Halichoeres_trimaculatus 0.1046071
    ## 283      Stegastes_nigricans Halichoeres_trimaculatus 0.1046071
    ## 284     Aulostomus_chinensis Halichoeres_trimaculatus 0.1046071
    ## 285        Epinephelus_merra Halichoeres_trimaculatus 0.1296265
    ## 286    Cephalopholis_urodeta Halichoeres_trimaculatus 0.1296265
    ## 287      Cephalopholis_argus Halichoeres_trimaculatus 0.1296265
    ## 288         Chaetodon_auriga Halichoeres_trimaculatus 0.1296265
    ## 289    Chaetodon_reticulatus Halichoeres_trimaculatus 0.1296265
    ## 290   Chaetodon_ornatissimus Halichoeres_trimaculatus 0.1296265
    ## 291    Ctenochaetus_striatus Halichoeres_trimaculatus 0.1296265
    ## 292         Zebrasoma_scopas Halichoeres_trimaculatus 0.1296265
    ## 293           Naso_lituratus Halichoeres_trimaculatus 0.1296265
    ## 294   Rhinecanthus_aculeatus Halichoeres_trimaculatus 0.1296265
    ## 295 Halichoeres_trimaculatus Halichoeres_trimaculatus 1.0000000
    ## 296       Chlorurus_spilurus Halichoeres_trimaculatus 0.4350973
    ## 297      Epibulus_insidiator Halichoeres_trimaculatus 0.4350973
    ## 298  Sargocentron_spiniferum Halichoeres_trimaculatus 0.0000000
    ## 299      Myripristis_berndti Halichoeres_trimaculatus 0.0000000
    ## 300     Myripristis_violacea Halichoeres_trimaculatus 0.0000000
    ## 301    Dascyllus_flavicaudus       Chlorurus_spilurus 0.1046071
    ## 302   Abudefduf_sexfasciatus       Chlorurus_spilurus 0.1046071
    ## 303      Stegastes_nigricans       Chlorurus_spilurus 0.1046071
    ## 304     Aulostomus_chinensis       Chlorurus_spilurus 0.1046071
    ## 305        Epinephelus_merra       Chlorurus_spilurus 0.1296265
    ## 306    Cephalopholis_urodeta       Chlorurus_spilurus 0.1296265
    ## 307      Cephalopholis_argus       Chlorurus_spilurus 0.1296265
    ## 308         Chaetodon_auriga       Chlorurus_spilurus 0.1296265
    ## 309    Chaetodon_reticulatus       Chlorurus_spilurus 0.1296265
    ## 310   Chaetodon_ornatissimus       Chlorurus_spilurus 0.1296265
    ## 311    Ctenochaetus_striatus       Chlorurus_spilurus 0.1296265
    ## 312         Zebrasoma_scopas       Chlorurus_spilurus 0.1296265
    ## 313           Naso_lituratus       Chlorurus_spilurus 0.1296265
    ## 314   Rhinecanthus_aculeatus       Chlorurus_spilurus 0.1296265
    ## 315 Halichoeres_trimaculatus       Chlorurus_spilurus 0.4350973
    ## 316       Chlorurus_spilurus       Chlorurus_spilurus 1.0000000
    ## 317      Epibulus_insidiator       Chlorurus_spilurus 0.5139098
    ## 318  Sargocentron_spiniferum       Chlorurus_spilurus 0.0000000
    ## 319      Myripristis_berndti       Chlorurus_spilurus 0.0000000
    ## 320     Myripristis_violacea       Chlorurus_spilurus 0.0000000
    ## 321    Dascyllus_flavicaudus      Epibulus_insidiator 0.1046071
    ## 322   Abudefduf_sexfasciatus      Epibulus_insidiator 0.1046071
    ## 323      Stegastes_nigricans      Epibulus_insidiator 0.1046071
    ## 324     Aulostomus_chinensis      Epibulus_insidiator 0.1046071
    ## 325        Epinephelus_merra      Epibulus_insidiator 0.1296265
    ## 326    Cephalopholis_urodeta      Epibulus_insidiator 0.1296265
    ## 327      Cephalopholis_argus      Epibulus_insidiator 0.1296265
    ## 328         Chaetodon_auriga      Epibulus_insidiator 0.1296265
    ## 329    Chaetodon_reticulatus      Epibulus_insidiator 0.1296265
    ## 330   Chaetodon_ornatissimus      Epibulus_insidiator 0.1296265
    ## 331    Ctenochaetus_striatus      Epibulus_insidiator 0.1296265
    ## 332         Zebrasoma_scopas      Epibulus_insidiator 0.1296265
    ## 333           Naso_lituratus      Epibulus_insidiator 0.1296265
    ## 334   Rhinecanthus_aculeatus      Epibulus_insidiator 0.1296265
    ## 335 Halichoeres_trimaculatus      Epibulus_insidiator 0.4350973
    ## 336       Chlorurus_spilurus      Epibulus_insidiator 0.5139098
    ## 337      Epibulus_insidiator      Epibulus_insidiator 1.0000000
    ## 338  Sargocentron_spiniferum      Epibulus_insidiator 0.0000000
    ## 339      Myripristis_berndti      Epibulus_insidiator 0.0000000
    ## 340     Myripristis_violacea      Epibulus_insidiator 0.0000000
    ## 341    Dascyllus_flavicaudus  Sargocentron_spiniferum 0.0000000
    ## 342   Abudefduf_sexfasciatus  Sargocentron_spiniferum 0.0000000
    ## 343      Stegastes_nigricans  Sargocentron_spiniferum 0.0000000
    ## 344     Aulostomus_chinensis  Sargocentron_spiniferum 0.0000000
    ## 345        Epinephelus_merra  Sargocentron_spiniferum 0.0000000
    ## 346    Cephalopholis_urodeta  Sargocentron_spiniferum 0.0000000
    ## 347      Cephalopholis_argus  Sargocentron_spiniferum 0.0000000
    ## 348         Chaetodon_auriga  Sargocentron_spiniferum 0.0000000
    ## 349    Chaetodon_reticulatus  Sargocentron_spiniferum 0.0000000
    ## 350   Chaetodon_ornatissimus  Sargocentron_spiniferum 0.0000000
    ## 351    Ctenochaetus_striatus  Sargocentron_spiniferum 0.0000000
    ## 352         Zebrasoma_scopas  Sargocentron_spiniferum 0.0000000
    ## 353           Naso_lituratus  Sargocentron_spiniferum 0.0000000
    ## 354   Rhinecanthus_aculeatus  Sargocentron_spiniferum 0.0000000
    ## 355 Halichoeres_trimaculatus  Sargocentron_spiniferum 0.0000000
    ## 356       Chlorurus_spilurus  Sargocentron_spiniferum 0.0000000
    ## 357      Epibulus_insidiator  Sargocentron_spiniferum 0.0000000
    ## 358  Sargocentron_spiniferum  Sargocentron_spiniferum 1.0000000
    ## 359      Myripristis_berndti  Sargocentron_spiniferum 0.5678036
    ## 360     Myripristis_violacea  Sargocentron_spiniferum 0.5678036
    ## 361    Dascyllus_flavicaudus      Myripristis_berndti 0.0000000
    ## 362   Abudefduf_sexfasciatus      Myripristis_berndti 0.0000000
    ## 363      Stegastes_nigricans      Myripristis_berndti 0.0000000
    ## 364     Aulostomus_chinensis      Myripristis_berndti 0.0000000
    ## 365        Epinephelus_merra      Myripristis_berndti 0.0000000
    ## 366    Cephalopholis_urodeta      Myripristis_berndti 0.0000000
    ## 367      Cephalopholis_argus      Myripristis_berndti 0.0000000
    ## 368         Chaetodon_auriga      Myripristis_berndti 0.0000000
    ## 369    Chaetodon_reticulatus      Myripristis_berndti 0.0000000
    ## 370   Chaetodon_ornatissimus      Myripristis_berndti 0.0000000
    ## 371    Ctenochaetus_striatus      Myripristis_berndti 0.0000000
    ## 372         Zebrasoma_scopas      Myripristis_berndti 0.0000000
    ## 373           Naso_lituratus      Myripristis_berndti 0.0000000
    ## 374   Rhinecanthus_aculeatus      Myripristis_berndti 0.0000000
    ## 375 Halichoeres_trimaculatus      Myripristis_berndti 0.0000000
    ## 376       Chlorurus_spilurus      Myripristis_berndti 0.0000000
    ## 377      Epibulus_insidiator      Myripristis_berndti 0.0000000
    ## 378  Sargocentron_spiniferum      Myripristis_berndti 0.5678036
    ## 379      Myripristis_berndti      Myripristis_berndti 1.0000000
    ## 380     Myripristis_violacea      Myripristis_berndti 0.9513783
    ## 381    Dascyllus_flavicaudus     Myripristis_violacea 0.0000000
    ## 382   Abudefduf_sexfasciatus     Myripristis_violacea 0.0000000
    ## 383      Stegastes_nigricans     Myripristis_violacea 0.0000000
    ## 384     Aulostomus_chinensis     Myripristis_violacea 0.0000000
    ## 385        Epinephelus_merra     Myripristis_violacea 0.0000000
    ## 386    Cephalopholis_urodeta     Myripristis_violacea 0.0000000
    ## 387      Cephalopholis_argus     Myripristis_violacea 0.0000000
    ## 388         Chaetodon_auriga     Myripristis_violacea 0.0000000
    ## 389    Chaetodon_reticulatus     Myripristis_violacea 0.0000000
    ## 390   Chaetodon_ornatissimus     Myripristis_violacea 0.0000000
    ## 391    Ctenochaetus_striatus     Myripristis_violacea 0.0000000
    ## 392         Zebrasoma_scopas     Myripristis_violacea 0.0000000
    ## 393           Naso_lituratus     Myripristis_violacea 0.0000000
    ## 394   Rhinecanthus_aculeatus     Myripristis_violacea 0.0000000
    ## 395 Halichoeres_trimaculatus     Myripristis_violacea 0.0000000
    ## 396       Chlorurus_spilurus     Myripristis_violacea 0.0000000
    ## 397      Epibulus_insidiator     Myripristis_violacea 0.0000000
    ## 398  Sargocentron_spiniferum     Myripristis_violacea 0.5678036
    ## 399      Myripristis_berndti     Myripristis_violacea 0.9513783
    ## 400     Myripristis_violacea     Myripristis_violacea 1.0000000

``` r
write.csv(rmf,file="relatedness_matrix_fish_melted.csv")
```

### In between these steps I basically generate a microbiome distance matrix in python and join the 2 matrices together to line up each species with relatedness value and microbiome distance value which become X and Y..this is the resulting data file with 100k points so you may need to subset a small test file to play with first

``` r
setwd("C:/Users/samde/OneDrive - UCLA IT Services/Fish Project/ISLAND COMPILED PROJECT/R/Relatedness")
```

``` r
rm=read.csv('relatedness_vs_microbiome_full_fish.csv', header=TRUE, sep=',')
head(rm)
```

    ##                 SpeciesA               SpeciesB
    ## 1 Chaetodon_ornatissimus Chaetodon_ornatissimus
    ## 2    Cephalopholis_argus    Cephalopholis_argus
    ## 3       Zebrasoma_scopas       Zebrasoma_scopas
    ## 4 Abudefduf_sexfasciatus Abudefduf_sexfasciatus
    ## 5   Myripristis_violacea   Myripristis_violacea
    ## 6         Naso_lituratus         Naso_lituratus
    ##                                         concat relatedness similarity
    ## 1 Chaetodon_ornatissimusChaetodon_ornatissimus           1          1
    ## 2       Cephalopholis_argusCephalopholis_argus           1          1
    ## 3             Zebrasoma_scopasZebrasoma_scopas           1          1
    ## 4 Abudefduf_sexfasciatusAbudefduf_sexfasciatus           1          1
    ## 5     Myripristis_violaceaMyripristis_violacea           1          1
    ## 6                 Naso_lituratusNaso_lituratus           1          1
    ##         dietA       dietB               dietAandB
    ## 1 corallivore corallivore corallivore:corallivore
    ## 2   carnivore   carnivore     carnivore:carnivore
    ## 3   herbivore   herbivore     herbivore:herbivore
    ## 4 planktivore planktivore planktivore:planktivore
    ## 5   carnivore   carnivore     carnivore:carnivore
    ## 6   herbivore   herbivore     herbivore:herbivore

``` r
library(ggplot2)
```

``` r
ggplot(rm, aes(x=relatedness, y=similarity, color=dietA)) +
  geom_point()+
  geom_jitter(width=0.015)
```

![](relatedness_vs_microbiome_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggplot(rm, aes(x=relatedness, y=similarity, color=dietB)) +
  geom_point()+
  geom_jitter(width=0.015)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)
```

    ## `geom_smooth()` using formula 'y ~ x'

![](relatedness_vs_microbiome_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
\#\#\# So the two above plots I basiclly colored it by only dietA or
diet B separately but as you can see they look identicaly which tells me
that color coding by every diet comparison is not necessary. See below
(or downlowd the html file I linked). As you can see it still is
identical.

### SO what lines tell us is the herbivory and corallivory exhibit the most phylosymbiosis (i.e. herbivore species have more similar gut microbiomes within their own species and close relatives compared to carnivores)

``` r
ggplot(rm, aes(x=relatedness, y=similarity, color=dietAandB)) +
  geom_point()+
  geom_jitter(width=0.015)
```

![](relatedness_vs_microbiome_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
