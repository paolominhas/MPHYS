# The breakdown of the data files here:

- KrakowScatter.root - The processed data from the Krakow simulation, 1 million events - final data
- HIBEAMScatter.root - The processed data from the HIBEAM simulation, 1 million events - final data















- hibeam_pions.root - The raw data from HIBEAM detector 100000 events, 1.5 GB

```cpp
Readout from processing hibeam_pions.root with the analysis script:

The input file contains source data.
The input file contains target data!
The input file contains TPC data.
The input file contains SEC data.
Begin analysis.
0 events read.
10000 events read.
20000 events read.
30000 events read.
40000 events read.
50000 events read.
60000 events read.
70000 events read.
80000 events read.
90000 events read.

******************************************************************************
*Tree    :hibeam    : Processed HIBEAM simulation output                     *
*Entries :   100000 : Total =      4196322517 bytes  File  Size =  635203118 *
*        :          : Tree compression factor =   6.53                       *
******************************************************************************
*Br    0 :PrimaryTrackID : vector<int>                                       *
*Entries :   100000 : Total  Size=    1803176 bytes  File Size  =     150092 *
*Baskets :       23 : Basket Size=     112128 bytes  Compression=  11.44     *
*............................................................................*
*Br    1 :PrimaryPDG : vector<int>                                           *
*Entries :   100000 : Total  Size=    1803060 bytes  File Size  =     149973 *
*Baskets :       23 : Basket Size=     112128 bytes  Compression=  11.45     *
*............................................................................*
*Br    2 :PrimaryEkin : vector<double>                                       *
*Entries :   100000 : Total  Size=    2203089 bytes  File Size  =     953218 *
*Baskets :       23 : Basket Size=     128000 bytes  Compression=   2.20     *
*............................................................................*
*Br    3 :PrimaryTime : vector<double>                                       *
*Entries :   100000 : Total  Size=    2203089 bytes  File Size  =     156304 *
*Baskets :       23 : Basket Size=     128000 bytes  Compression=  13.42     *
*............................................................................*
*Br    4 :PrimaryPosX : vector<double>                                       *
*Entries :   100000 : Total  Size=    2203089 bytes  File Size  =    1007282 *
*Baskets :       23 : Basket Size=     128000 bytes  Compression=   2.08     *
*............................................................................*
*Br    5 :PrimaryPosY : vector<double>                                       *
*Entries :   100000 : Total  Size=    2203089 bytes  File Size  =    1006834 *
*Baskets :       23 : Basket Size=     128000 bytes  Compression=   2.08     *
*............................................................................*
*Br    6 :PrimaryPosZ : vector<double>                                       *
*Entries :   100000 : Total  Size=    2203089 bytes  File Size  =     998187 *
*Baskets :       23 : Basket Size=     128000 bytes  Compression=   2.10     *
*............................................................................*
*Br    7 :PrimaryMomX : vector<double>                                       *
*Entries :   100000 : Total  Size=    2203089 bytes  File Size  =    1003346 *
*Baskets :       23 : Basket Size=     128000 bytes  Compression=   2.09     *
*............................................................................*
*Br    8 :PrimaryMomY : vector<double>                                       *
*Entries :   100000 : Total  Size=    2203089 bytes  File Size  =    1003640 *
*Baskets :       23 : Basket Size=     128000 bytes  Compression=   2.09     *
*............................................................................*
*Br    9 :PrimaryMomZ : vector<double>                                       *
*Entries :   100000 : Total  Size=    2203089 bytes  File Size  =    1003633 *
*Baskets :       23 : Basket Size=     128000 bytes  Compression=   2.09     *
*............................................................................*
*Br   10 :PrimaryWeight : vector<double>                                     *
*Entries :   100000 : Total  Size=    2203147 bytes  File Size  =     156622 *
*Baskets :       23 : Basket Size=     128000 bytes  Compression=  13.40     *
*............................................................................*
*Br   11 :target_TrackID : vector<int>                                       *
*Entries :   100000 : Total  Size=    1562833 bytes  File Size  =     284825 *
*Baskets :       22 : Basket Size=     102400 bytes  Compression=   5.23     *
*............................................................................*
*Br   12 :target_LayerID : vector<int>                                       *
*Entries :   100000 : Total  Size=    1562833 bytes  File Size  =     274853 *
*Baskets :       22 : Basket Size=     102400 bytes  Compression=   5.41     *
*............................................................................*
*Br   13 :target_LayerID1 : vector<int>                                      *
*Entries :   100000 : Total  Size=    1562861 bytes  File Size  =     273694 *
*Baskets :       22 : Basket Size=     102400 bytes  Compression=   5.44     *
*............................................................................*
*Br   14 :target_LayerID2 : vector<int>                                      *
*Entries :   100000 : Total  Size=    1562861 bytes  File Size  =     273694 *
*Baskets :       22 : Basket Size=     102400 bytes  Compression=   5.44     *
*............................................................................*
*Br   15 :target_PDG : vector<int>                                           *
*Entries :   100000 : Total  Size=    1562721 bytes  File Size  =     284475 *
*Baskets :       22 : Basket Size=     102400 bytes  Compression=   5.23     *
*............................................................................*
*Br   16 :target_EDep : vector<double>                                       *
*Entries :   100000 : Total  Size=    1722633 bytes  File Size  =     585799 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   2.80     *
*............................................................................*
*Br   17 :target_Time : vector<double>                                       *
*Entries :   100000 : Total  Size=    1722633 bytes  File Size  =     586076 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   2.80     *
*............................................................................*
*Br   18 :target_TrackLength : vector<double>                                *
*Entries :   100000 : Total  Size=    1722836 bytes  File Size  =     585561 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   2.80     *
*............................................................................*
*Br   19 :target_Position_X : vector<double>                                 *
*Entries :   100000 : Total  Size=    1722807 bytes  File Size  =     602218 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   2.72     *
*............................................................................*
*Br   20 :target_Position_Y : vector<double>                                 *
*Entries :   100000 : Total  Size=    1722807 bytes  File Size  =     602113 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   2.73     *
*............................................................................*
*Br   21 :target_Position_Z : vector<double>                                 *
*Entries :   100000 : Total  Size=    1722807 bytes  File Size  =     308674 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   5.32     *
*............................................................................*
*Br   22 :target_GlobalPosition_X : vector<double>                           *
*Entries :   100000 : Total  Size=    1722981 bytes  File Size  =     602542 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   2.72     *
*............................................................................*
*Br   23 :target_GlobalPosition_Y : vector<double>                           *
*Entries :   100000 : Total  Size=    1722981 bytes  File Size  =     602392 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   2.72     *
*............................................................................*
*Br   24 :target_GlobalPosition_Z : vector<double>                           *
*Entries :   100000 : Total  Size=    1722981 bytes  File Size  =     299732 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   5.47     *
*............................................................................*
*Br   25 :target_Momentum_X : vector<double>                                 *
*Entries :   100000 : Total  Size=    1722807 bytes  File Size  =     597221 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   2.75     *
*............................................................................*
*Br   26 :target_Momentum_Y : vector<double>                                 *
*Entries :   100000 : Total  Size=    1722807 bytes  File Size  =     597232 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   2.75     *
*............................................................................*
*Br   27 :target_Momentum_Z : vector<double>                                 *
*Entries :   100000 : Total  Size=    1722807 bytes  File Size  =     585891 *
*Baskets :       23 : Basket Size=     109056 bytes  Compression=   2.80     *
*............................................................................*
*Branch  :SEC                                                                *
*Entries :   100000 : BranchElement (see below)                              *
*............................................................................*
*Br   28 :TObject   : BASE                                                   *
*Entries :   100000 : Total  Size=    1405889 bytes  File Size  =     153182 *
*Baskets :       50 : Basket Size=      52736 bytes  Compression=   9.00     *
*............................................................................*
*Br   29 :nHits     : Int_t                                                  *
*Entries :   100000 : Total  Size=     402629 bytes  File Size  =     131963 *
*Baskets :       22 : Basket Size=      51200 bytes  Compression=   2.90     *
*............................................................................*
*Br   30 :ringNo    : vector<int>                                            *
*Entries :   100000 : Total  Size=   13486320 bytes  File Size  =    3637014 *
*Baskets :      163 : Basket Size=     175616 bytes  Compression=   3.67     *
*............................................................................*
*Br   31 :crystalNo : vector<int>                                            *
*Entries :   100000 : Total  Size=   13486827 bytes  File Size  =    4337592 *
*Baskets :      163 : Basket Size=     175616 bytes  Compression=   3.08     *
*............................................................................*
*Br   32 :Edep      : vector<double>                                         *
*Entries :   100000 : Total  Size=   25564248 bytes  File Size  =   21305595 *
*Baskets :      257 : Basket Size=     298496 bytes  Compression=   1.20     *
*............................................................................*
*Br   33 :Time      : vector<double>                                         *
*Entries :   100000 : Total  Size=   25564248 bytes  File Size  =   23876864 *
*Baskets :      257 : Basket Size=     298496 bytes  Compression=   1.07     *
*............................................................................*
*Br   34 :Pos_X     : vector<double>                                         *
*Entries :   100000 : Total  Size=   25564617 bytes  File Size  =   24071706 *
*Baskets :      258 : Basket Size=     299008 bytes  Compression=   1.06     *
*............................................................................*
*Br   35 :Pos_Y     : vector<double>                                         *
*Entries :   100000 : Total  Size=   25564617 bytes  File Size  =   24028618 *
*Baskets :      258 : Basket Size=     299008 bytes  Compression=   1.06     *
*............................................................................*
*Br   36 :Pos_Z     : vector<double>                                         *
*Entries :   100000 : Total  Size=   25564617 bytes  File Size  =   22461573 *
*Baskets :      258 : Basket Size=     299008 bytes  Compression=   1.14     *
*............................................................................*
*Branch  :TPC                                                                *
*Entries :   100000 : BranchElement (see below)                              *
*............................................................................*
*Br   37 :TObject   : BASE                                                   *
*Entries :   100000 : Total  Size=    1405889 bytes  File Size  =     153182 *
*Baskets :       50 : Basket Size=      52736 bytes  Compression=   9.00     *
*............................................................................*
*Br   38 :nHits     : Int_t                                                  *
*Entries :   100000 : Total  Size=     402629 bytes  File Size  =      94263 *
*Baskets :       22 : Basket Size=      51200 bytes  Compression=   4.06     *
*............................................................................*
*Br   39 :Edep      : vector<double>                                         *
*Entries :   100000 : Total  Size=    5596824 bytes  File Size  =    4509627 *
*Baskets :      105 : Basket Size=      95232 bytes  Compression=   1.22     *
*............................................................................*
*Br   40 :Time      : vector<double>                                         *
*Entries :   100000 : Total  Size=    1405721 bytes  File Size  =     152095 *
*Baskets :       50 : Basket Size=      52736 bytes  Compression=   9.07     *
*............................................................................*
*Br   41 :Pos_X     : vector<double>                                         *
*Entries :   100000 : Total  Size=    1405777 bytes  File Size  =     152664 *
*Baskets :       50 : Basket Size=      52736 bytes  Compression=   9.03     *
*............................................................................*
*Br   42 :Pos_Y     : vector<double>                                         *
*Entries :   100000 : Total  Size=    1405777 bytes  File Size  =     152664 *
*Baskets :       50 : Basket Size=      52736 bytes  Compression=   9.03     *
*............................................................................*
*Br   43 :Pos_Z     : vector<double>                                         *
*Entries :   100000 : Total  Size=    1405777 bytes  File Size  =     152664 *
*Baskets :       50 : Basket Size=      52736 bytes  Compression=   9.03     *
*............................................................................*
*Br   44 :padRow    : vector<int>                                            *
*Entries :   100000 : Total  Size= 1323747171 bytes  File Size  =   60088740 *
*Baskets :     3532 : Basket Size=   25600000 bytes  Compression=  21.78     *
*............................................................................*
*Br   45 :padColumn : vector<int>                                            *
*Entries :   100000 : Total  Size= 1323757785 bytes  File Size  =  108877350 *
*Baskets :     3532 : Basket Size=   25600000 bytes  Compression=  12.02     *
*............................................................................*
*Br   46 :timestamp : vector<int>                                            *
*Entries :   100000 : Total  Size= 1323757785 bytes  File Size  =  319654073 *
*Baskets :     3532 : Basket Size=   25600000 bytes  Compression=   4.09     *
*............................................................................*
*Br   47 :nEl       : vector<int>                                            *
*Entries :   100000 : Total  Size=    3501761 bytes  File Size  =    1561271 *
*Baskets :       83 : Basket Size=      73728 bytes  Compression=   2.22     *
*............................................................................*
Finished. 100000 events read in total.
```