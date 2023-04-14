# openFoamUtils

## turbulenceModels
### SARC
Модель турбулентности Спаларта-Аллмараса с поправкой на кривизну линий тока:
Shur, M.L. Turbulence Modeling in Rotating and Curved Channels: Assessing the Spalart-Shur Correction / M.L.  Shur,  M.K.  Streles,  A.K.  Travin,  P.R.  Spalart  //  AIAA  Journal. – 2000. – V. 38. – No 5. –  P. 784 – 792. 

### SARCM
Модель турбулентности Спаларта-Аллмараса с поправкой на кривизну линий тока (модифицированная):
Qiang,  Z.  A  new  simpler  rotation/curvature  correction method for Spalart–Allmaras turbulence model / Qiang Zhang, Yong Yang // Chinese Journal of Aeronautics. – 2013. – V. 26. –  Issue 2. – 2013. – P. 326 – 333. 

## kOmegaSSTRC
Модель турбулентности k-omega-sst с поправкой на кривизну линий тока:
Smirnov, P. Sensitization of the SST turbulence model to rotation and curvature by applying the Spalart – Shur correction  term  /  P.  Smirnov,  F.  Menter  //  Proc.  ASME  Turbo  Expo 2008: Power for Land, Sea and Air, 2008. 

## kOmegaSSTRCHellsten
Модель турбулентности k-omega-sst с поправкой на кривизну линий тока в представлении Hellsten
Hellsten,  A.  Some  Improvements  in  Menter's  k-omega SST Turbulence Model / A. Hellsten // AIAA-98-2554, 1998.

## Utilities
### aspectRatio
Утилита для создания векторного поля с отношением сторон ячеек для последующего направленного измельчения

### deltasCalc
Утилита для расчета дельт сетки

### heatCalc
Утилита для расчета теплового потока от поверхности по полю градиентов температуры и температуропроводности

### yPlusRAS
Расчет y+ поверхностей (в OpenFOAM используется y* вместо y+)

## Solvers
### prandtlSolver
Решатель для первой в мире модели турбулентности. Практической ценности нет - только в образовательных целях
