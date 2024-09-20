class Primer {
  constructor(seq, template) {
    // Tables of nearest-neighbor thermodynamics for DNA bases

    // Hybridization match table of dH and dS
    this.dH_mm = [
      8888.0, 8888.0, 8888.0, 4700.0, 8888.0, 8888.0, // AA/AA AA/AC AA/AG AA/AT AA/AN AA/A- 
      8888.0, 8888.0, 8888.0, 7600.0, 8888.0, 8888.0, // AA/CA AA/CC AA/CG AA/CT AA/CN AA/C- 
      8888.0, 8888.0, 8888.0, 3000.0, 8888.0, 8888.0, // AA/GA AA/GC AA/GG AA/GT AA/GN AA/G- 
      1200.0, 2300.0, -600.0, -7900, 2300.0, 2300.0, // AA/TA AA/TC AA/TG AA/TT AA/TN AA/T- 
      8888.0, 8888.0, 8888.0, 7600.0, 8888.0, 8888.0, // AA/NA AA/NC AA/NG AA/NT AA/NN AA/N- 
      8888.0, 8888.0, 8888.0, 7600.0, 8888.0, 8888.0, // AA/-A AA/-C AA/-G AA/-T AA/-N AA/-- 
      8888.0, 8888.0, -2900.0, 8888.0, 8888.0, 8888.0, // AC/AA AC/AC AC/AG AC/AT AC/AN AC/A- 
      8888.0, 8888.0, -700.0, 8888.0, 8888.0, 8888.0, // AC/CA AC/CC AC/CG AC/CT AC/CN AC/C- 
      8888.0, 8888.0, 500.0, 8888.0, 8888.0, 8888.0, // AC/GA AC/GC AC/GG AC/GT AC/GN AC/G- 
      5300.0, 0.0, -8400, 700.0, 5300.0, 5300.0, // AC/TA AC/TC AC/TG AC/TT AC/TN AC/T- 
      8888.0, 8888.0, 500.0, 8888.0, 8888.0, 8888.0, // AC/NA AC/NC AC/NG AC/NT AC/NN AC/N- 
      8888.0, 8888.0, 500.0, 8888.0, 8888.0, 8888.0, // AC/-A AC/-C AC/-G AC/-T AC/-N AC/-- 
      8888.0, -900.0, 8888.0, 8888.0, 8888.0, 8888.0, // AG/AA AG/AC AG/AG AG/AT AG/AN AG/A- 
      8888.0, 600.0, 8888.0, 8888.0, 8888.0, 8888.0, // AG/CA AG/CC AG/CG AG/CT AG/CN AG/C- 
      8888.0, -4000.0, 8888.0, 8888.0, 8888.0, 8888.0, // AG/GA AG/GC AG/GG AG/GT AG/GN AG/G- 
      -700.0, -7800, -3100.0, 1000.0, 1000.0, 1000.0, // AG/TA AG/TC AG/TG AG/TT AG/TN AG/T- 
      8888.0, 600.0, 8888.0, 8888.0, 8888.0, 8888.0, // AG/NA AG/NC AG/NG AG/NT AG/NN AG/N- 
      8888.0, 600.0, 8888.0, 8888.0, 8888.0, 8888.0, // AG/-A AG/-C AG/-G AG/-T AG/-N AG/-- 
      1200.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // AT/AA AT/AC AT/AG AT/AT AT/AN AT/A- 
      5300.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // AT/CA AT/CC AT/CG AT/CT AT/CN AT/C- 
      -700.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // AT/GA AT/GC AT/GG AT/GT AT/GN AT/G- 
      -7200, -1200.0, -4000.0, -2700.0, 0, 0, // AT/TA AT/TC AT/TG AT/TT AT/TN AT/T- 
      5300.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // AT/NA AT/NC AT/NG AT/NT AT/NN AT/N- 
      5300.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // AT/-A AT/-C AT/-G AT/-T AT/-N AT/-- 
      1200.0, -900.0, -2900.0, 4700.0, 8888.0, 8888.0, // AN/AA AN/AC AN/AG AN/AT AN/AN AN/A- 
      5300.0, 600.0, -700.0, 7600.0, 8888.0, 8888.0, // AN/CA AN/CC AN/CG AN/CT AN/CN AN/C- 
      -700.0, -4000.0, 500.0, 3000.0, 8888.0, 8888.0, // AN/GA AN/GC AN/GG AN/GT AN/GN AN/G- 
      -7200, -7800, -8400, -7900, 0, 0, // AN/TA AN/TC AN/TG AN/TT AN/TN AN/T- 
      5300.0, 600.0, 500.0, 7600.0, 8888.0, 8888.0, // AN/NA AN/NC AN/NG AN/NT AN/NN AN/N- 
      5300.0, 600.0, 500.0, 7600.0, 8888.0, 8888.0, // AN/-A AN/-C AN/-G AN/-T AN/-N AN/-- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // A-/AA A-/AC A-/AG A-/AT A-/AN A-/A- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // A-/CA A-/CC A-/CG A-/CT A-/CN A-/C- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // A-/GA A-/GC A-/GG A-/GT A-/GN A-/G- 
      5300.0, 2300.0, -600.0, 1000.0, 5300.0, 5300.0, // A-/TA A-/TC A-/TG A-/TT A-/TN A-/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // A-/NA A-/NC A-/NG A-/NT A-/NN A-/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // A-/-A A-/-C A-/-G A-/-T A-/-N A-/-- 
      8888.0, 8888.0, 8888.0, 3400.0, 8888.0, 8888.0, // CA/AA CA/AC CA/AG CA/AT CA/AN CA/A- 
      8888.0, 8888.0, 8888.0, 6100.0, 8888.0, 8888.0, // CA/CA CA/CC CA/CG CA/CT CA/CN CA/C- 
      -900.0, 1900.0, -700.0, -8500, 1900.0, 1900.0, // CA/GA CA/GC CA/GG CA/GT CA/GN CA/G- 
      8888.0, 8888.0, 8888.0, 1000.0, 8888.0, 8888.0, // CA/TA CA/TC CA/TG CA/TT CA/TN CA/T- 
      8888.0, 8888.0, 8888.0, 6100.0, 8888.0, 8888.0, // CA/NA CA/NC CA/NG CA/NT CA/NN CA/N- 
      8888.0, 8888.0, 8888.0, 6100.0, 8888.0, 8888.0, // CA/-A CA/-C CA/-G CA/-T CA/-N CA/-- 
      8888.0, 8888.0, 5200.0, 8888.0, 8888.0, 8888.0, // CC/AA CC/AC CC/AG CC/AT CC/AN CC/A- 
      8888.0, 8888.0, 3600.0, 8888.0, 8888.0, 8888.0, // CC/CA CC/CC CC/CG CC/CT CC/CN CC/C- 
      600.0, -1500.0, -8000, -800.0, 600.0, 600.0, // CC/GA CC/GC CC/GG CC/GT CC/GN CC/G- 
      8888.0, 8888.0, 5200.0, 8888.0, 8888.0, 8888.0, // CC/TA CC/TC CC/TG CC/TT CC/TN CC/T- 
      8888.0, 8888.0, 5200.0, 8888.0, 8888.0, 8888.0, // CC/NA CC/NC CC/NG CC/NT CC/NN CC/N- 
      8888.0, 8888.0, 5200.0, 8888.0, 8888.0, 8888.0, // CC/-A CC/-C CC/-G CC/-T CC/-N CC/-- 
      8888.0, 1900.0, 8888.0, 8888.0, 8888.0, 8888.0, // CG/AA CG/AC CG/AG CG/AT CG/AN CG/A- 
      8888.0, -1500.0, 8888.0, 8888.0, 8888.0, 8888.0, // CG/CA CG/CC CG/CG CG/CT CG/CN CG/C- 
      -4000.0, -10600, -4900.0, -4100.0, 0, 0, // CG/GA CG/GC CG/GG CG/GT CG/GN CG/G- 
      8888.0, -1500.0, 8888.0, 8888.0, 8888.0, 8888.0, // CG/TA CG/TC CG/TG CG/TT CG/TN CG/T- 
      8888.0, 1900.0, 8888.0, 8888.0, 8888.0, 8888.0, // CG/NA CG/NC CG/NG CG/NT CG/NN CG/N- 
      8888.0, 1900.0, 8888.0, 8888.0, 8888.0, 8888.0, // CG/-A CG/-C CG/-G CG/-T CG/-N CG/-- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // CT/AA CT/AC CT/AG CT/AT CT/AN CT/A- 
      0.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // CT/CA CT/CC CT/CG CT/CT CT/CN CT/C- 
      -7800, -1500.0, -2800.0, -5000.0, 0, 0, // CT/GA CT/GC CT/GG CT/GT CT/GN CT/G- 
      -1200.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // CT/TA CT/TC CT/TG CT/TT CT/TN CT/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // CT/NA CT/NC CT/NG CT/NT CT/NN CT/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // CT/-A CT/-C CT/-G CT/-T CT/-N CT/-- 
      8888.0, 1900.0, 5200.0, 3400.0, 8888.0, 8888.0, // CN/AA CN/AC CN/AG CN/AT CN/AN CN/A- 
      0.0, -1500.0, 3600.0, 6100.0, 8888.0, 8888.0, // CN/CA CN/CC CN/CG CN/CT CN/CN CN/C- 
      -7800, -10600, -8000, -8500, 0, 0, // CN/GA CN/GC CN/GG CN/GT CN/GN CN/G- 
      -1200.0, -1500.0, 5200.0, 1000.0, 8888.0, 8888.0, // CN/TA CN/TC CN/TG CN/TT CN/TN CN/T- 
      8888.0, 1900.0, 5200.0, 6100.0, 8888.0, 8888.0, // CN/NA CN/NC CN/NG CN/NT CN/NN CN/N- 
      8888.0, 1900.0, 5200.0, 6100.0, 8888.0, 8888.0, // CN/-A CN/-C CN/-G CN/-T CN/-N CN/-- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // C-/AA C-/AC C-/AG C-/AT C-/AN C-/A- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // C-/CA C-/CC C-/CG C-/CT C-/CN C-/C- 
      600.0, 1900.0, -700.0, -800.0, 1900.0, 1900.0, // C-/GA C-/GC C-/GG C-/GT C-/GN C-/G- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // C-/TA C-/TC C-/TG C-/TT C-/TN C-/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // C-/NA C-/NC C-/NG C-/NT C-/NN C-/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // C-/-A C-/-C C-/-G C-/-T C-/-N C-/-- 
      8888.0, 8888.0, 8888.0, 700.0, 8888.0, 8888.0, // GA/AA GA/AC GA/AG GA/AT GA/AN GA/A- 
      -2900.0, 5200.0, -600.0, -8200, 5200.0, 5200.0, // GA/CA GA/CC GA/CG GA/CT GA/CN GA/C- 
      8888.0, 8888.0, 8888.0, 1600.0, 8888.0, 8888.0, // GA/GA GA/GC GA/GG GA/GT GA/GN GA/G- 
      8888.0, 8888.0, 8888.0, -1300.0, 8888.0, 8888.0, // GA/TA GA/TC GA/TG GA/TT GA/TN GA/T- 
      8888.0, 8888.0, 8888.0, 1600.0, 8888.0, 8888.0, // GA/NA GA/NC GA/NG GA/NT GA/NN GA/N- 
      8888.0, 8888.0, 8888.0, 1600.0, 8888.0, 8888.0, // GA/-A GA/-C GA/-G GA/-T GA/-N GA/-- 
      8888.0, 8888.0, -600.0, 8888.0, 8888.0, 8888.0, // GC/AA GC/AC GC/AG GC/AT GC/AN GC/A- 
      -700.0, 3600.0, -9800, 2300.0, 3600.0, 3600.0, // GC/CA GC/CC GC/CG GC/CT GC/CN GC/C- 
      8888.0, 8888.0, -6000.0, 8888.0, 8888.0, 8888.0, // GC/GA GC/GC GC/GG GC/GT GC/GN GC/G- 
      8888.0, 8888.0, -4400.0, 8888.0, 8888.0, 8888.0, // GC/TA GC/TC GC/TG GC/TT GC/TN GC/T- 
      8888.0, 8888.0, 0, 8888.0, 8888.0, 8888.0, // GC/NA GC/NC GC/NG GC/NT GC/NN GC/N- 
      8888.0, 8888.0, 0, 8888.0, 8888.0, 8888.0, // GC/-A GC/-C GC/-G GC/-T GC/-N GC/-- 
      8888.0, -700.0, 8888.0, 8888.0, 8888.0, 8888.0, // GG/AA GG/AC GG/AG GG/AT GG/AN GG/A- 
      500.0, -8000, -6000.0, 3300.0, 3300.0, 3300.0, // GG/CA GG/CC GG/CG GG/CT GG/CN GG/C- 
      8888.0, -4900.0, 8888.0, 8888.0, 8888.0, 8888.0, // GG/GA GG/GC GG/GG GG/GT GG/GN GG/G- 
      8888.0, -2800.0, 8888.0, 5800.0, 8888.0, 8888.0, // GG/TA GG/TC GG/TG GG/TT GG/TN GG/T- 
      8888.0, 0, 8888.0, 8888.0, 8888.0, 8888.0, // GG/NA GG/NC GG/NG GG/NT GG/NN GG/N- 
      8888.0, 0, 8888.0, 8888.0, 8888.0, 8888.0, // GG/-A GG/-C GG/-G GG/-T GG/-N GG/-- 
      -600.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // GT/AA GT/AC GT/AG GT/AT GT/AN GT/A- 
      -8400, 5200.0, -4400.0, -2200.0, 5200.0, 5200.0, // GT/CA GT/CC GT/CG GT/CT GT/CN GT/C- 
      -3100.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // GT/GA GT/GC GT/GG GT/GT GT/GN GT/G- 
      -4000.0, 8888.0, 4100.0, 8888.0, 8888.0, 8888.0, // GT/TA GT/TC GT/TG GT/TT GT/TN GT/T- 
      0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // GT/NA GT/NC GT/NG GT/NT GT/NN GT/N- 
      0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // GT/-A GT/-C GT/-G GT/-T GT/-N GT/-- 
      -600.0, -700.0, -600.0, 700.0, 8888.0, 8888.0, // GN/AA GN/AC GN/AG GN/AT GN/AN GN/A- 
      -8400, -8000, -9800, -8200, 3300.0, 3300.0, // GN/CA GN/CC GN/CG GN/CT GN/CN GN/C- 
      -3100.0, -4900.0, -6000.0, 1600.0, 8888.0, 8888.0, // GN/GA GN/GC GN/GG GN/GT GN/GN GN/G- 
      -4000.0, -2800.0, -4400.0, -1300.0, 8888.0, 8888.0, // GN/TA GN/TC GN/TG GN/TT GN/TN GN/T- 
      0, 0, 0, 1600.0, 8888.0, 8888.0, // GN/NA GN/NC GN/NG GN/NT GN/NN GN/N- 
      0, 0, 0, 1600.0, 8888.0, 8888.0, // GN/-A GN/-C GN/-G GN/-T GN/-N GN/-- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // G-/AA G-/AC G-/AG G-/AT G-/AN G-/A- 
      500.0, 5200.0, -600.0, 3300.0, 5200.0, 5200.0, // G-/CA G-/CC G-/CG G-/CT G-/CN G-/C- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // G-/GA G-/GC G-/GG G-/GT G-/GN G-/G- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // G-/TA G-/TC G-/TG G-/TT G-/TN G-/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // G-/NA G-/NC G-/NG G-/NT G-/NN G-/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // G-/-A G-/-C G-/-G G-/-T G-/-N G-/-- 
      4700.0, 3400.0, 700.0, -7200, 4700.0, 4700.0, // TA/AA TA/AC TA/AG TA/AT TA/AN TA/A- 
      8888.0, 8888.0, 8888.0, 1200.0, 8888.0, 8888.0, // TA/CA TA/CC TA/CG TA/CT TA/CN TA/C- 
      8888.0, 8888.0, 8888.0, -100.0, 8888.0, 8888.0, // TA/GA TA/GC TA/GG TA/GT TA/GN TA/G- 
      8888.0, 8888.0, -100.0, 200.0, 8888.0, 8888.0, // TA/TA TA/TC TA/TG TA/TT TA/TN TA/T- 
      8888.0, 8888.0, 8888.0, 1200.0, 8888.0, 8888.0, // TA/NA TA/NC TA/NG TA/NT TA/NN TA/N- 
      8888.0, 8888.0, 8888.0, 1200.0, 8888.0, 8888.0, // TA/-A TA/-C TA/-G TA/-T TA/-N TA/-- 
      7600.0, 6100.0, -8200, 1200.0, 7600.0, 7600.0, // TC/AA TC/AC TC/AG TC/AT TC/AN TC/A- 
      8888.0, 8888.0, 2300.0, 8888.0, 8888.0, 8888.0, // TC/CA TC/CC TC/CG TC/CT TC/CN TC/C- 
      8888.0, 8888.0, 3300.0, 8888.0, 8888.0, 8888.0, // TC/GA TC/GC TC/GG TC/GT TC/GN TC/G- 
      8888.0, 8888.0, -2200.0, 8888.0, 8888.0, 8888.0, // TC/TA TC/TC TC/TG TC/TT TC/TN TC/T- 
      8888.0, 8888.0, 3300.0, 8888.0, 8888.0, 8888.0, // TC/NA TC/NC TC/NG TC/NT TC/NN TC/N- 
      8888.0, 8888.0, 3300.0, 8888.0, 8888.0, 8888.0, // TC/-A TC/-C TC/-G TC/-T TC/-N TC/-- 
      3000.0, -8500, 1600.0, -100.0, 3000.0, 3000.0, // TG/AA TG/AC TG/AG TG/AT TG/AN TG/A- 
      8888.0, -800.0, 8888.0, 8888.0, 8888.0, 8888.0, // TG/CA TG/CC TG/CG TG/CT TG/CN TG/C- 
      8888.0, -4100.0, 8888.0, -1400.0, 8888.0, 8888.0, // TG/GA TG/GC TG/GG TG/GT TG/GN TG/G- 
      8888.0, -5000.0, 8888.0, 8888.0, 8888.0, 8888.0, // TG/TA TG/TC TG/TG TG/TT TG/TN TG/T- 
      8888.0, 0, 8888.0, 8888.0, 8888.0, 8888.0, // TG/NA TG/NC TG/NG TG/NT TG/NN TG/N- 
      8888.0, 0, 8888.0, 8888.0, 8888.0, 8888.0, // TG/-A TG/-C TG/-G TG/-T TG/-N TG/-- 
      -7900, 1000.0, -1300.0, 200.0, 1000.0, 1000.0, // TT/AA TT/AC TT/AG TT/AT TT/AN TT/A- 
      700.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // TT/CA TT/CC TT/CG TT/CT TT/CN TT/C- 
      1000.0, 8888.0, 5800.0, 8888.0, 8888.0, 8888.0, // TT/GA TT/GC TT/GG TT/GT TT/GN TT/G- 
      -2700.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // TT/TA TT/TC TT/TG TT/TT TT/TN TT/T- 
      1000.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // TT/NA TT/NC TT/NG TT/NT TT/NN TT/N- 
      1000.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // TT/-A TT/-C TT/-G TT/-T TT/-N TT/-- 
      -7900, -8500, -8200, -7200, 1000.0, 1000.0, // TN/AA TN/AC TN/AG TN/AT TN/AN TN/A- 
      700.0, -800.0, 2300.0, 1200.0, 8888.0, 8888.0, // TN/CA TN/CC TN/CG TN/CT TN/CN TN/C- 
      1000.0, -4100.0, 3300.0, -1400.0, 8888.0, 8888.0, // TN/GA TN/GC TN/GG TN/GT TN/GN TN/G- 
      -2700.0, -5000.0, -2200.0, 200.0, 8888.0, 8888.0, // TN/TA TN/TC TN/TG TN/TT TN/TN TN/T- 
      1000.0, 0, 3300.0, 1200.0, 8888.0, 8888.0, // TN/NA TN/NC TN/NG TN/NT TN/NN TN/N- 
      1000.0, 0, 3300.0, 1200.0, 8888.0, 8888.0, // TN/-A TN/-C TN/-G TN/-T TN/-N TN/-- 
      7600.0, 6100.0, 1600.0, 1200.0, 7600.0, 7600.0, // T-/AA T-/AC T-/AG T-/AT T-/AN T-/A- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // T-/CA T-/CC T-/CG T-/CT T-/CN T-/C- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // T-/GA T-/GC T-/GG T-/GT T-/GN T-/G- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // T-/TA T-/TC T-/TG T-/TT T-/TN T-/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // T-/NA T-/NC T-/NG T-/NT T-/NN T-/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // T-/-A T-/-C T-/-G T-/-T T-/-N T-/-- 
      4700.0, 3400.0, 700.0, -7200, 4700.0, 4700.0, // NA/AA NA/AC NA/AG NA/AT NA/AN NA/A- 
      -2900.0, 5200.0, -600.0, -8200, 5200.0, 5200.0, // NA/CA NA/CC NA/CG NA/CT NA/CN NA/C- 
      -900.0, 1900.0, -700.0, -8500, 1900.0, 1900.0, // NA/GA NA/GC NA/GG NA/GT NA/GN NA/G- 
      1200.0, 2300.0, -600.0, -7900, 2300.0, 2300.0, // NA/TA NA/TC NA/TG NA/TT NA/TN NA/T- 
      8888.0, 8888.0, 8888.0, 1200.0, 8888.0, 8888.0, // NA/NA NA/NC NA/NG NA/NT NA/NN NA/N- 
      8888.0, 8888.0, 8888.0, 1200.0, 8888.0, 8888.0, // NA/-A NA/-C NA/-G NA/-T NA/-N NA/-- 
      7600.0, 6100.0, -8200, 1200.0, 7600.0, 7600.0, // NC/AA NC/AC NC/AG NC/AT NC/AN NC/A- 
      -700.0, 3600.0, -9800, 2300.0, 3600.0, 3600.0, // NC/CA NC/CC NC/CG NC/CT NC/CN NC/C- 
      600.0, -1500.0, -8000, -800.0, 600.0, 600.0, // NC/GA NC/GC NC/GG NC/GT NC/GN NC/G- 
      5300.0, 0.0, -8400, 700.0, 5300.0, 5300.0, // NC/TA NC/TC NC/TG NC/TT NC/TN NC/T- 
      8888.0, 8888.0, 0, 8888.0, 8888.0, 8888.0, // NC/NA NC/NC NC/NG NC/NT NC/NN NC/N- 
      8888.0, 8888.0, 0, 8888.0, 8888.0, 8888.0, // NC/-A NC/-C NC/-G NC/-T NC/-N NC/-- 
      3000.0, -8500, 1600.0, -100.0, 3000.0, 3000.0, // NG/AA NG/AC NG/AG NG/AT NG/AN NG/A- 
      500.0, -8000, -6000.0, 3300.0, 3300.0, 3300.0, // NG/CA NG/CC NG/CG NG/CT NG/CN NG/C- 
      -4000.0, -10600, -4900.0, -4100.0, 0, 0, // NG/GA NG/GC NG/GG NG/GT NG/GN NG/G- 
      -700.0, -7800, -3100.0, 1000.0, 1000.0, 1000.0, // NG/TA NG/TC NG/TG NG/TT NG/TN NG/T- 
      8888.0, 0, 8888.0, 8888.0, 8888.0, 8888.0, // NG/NA NG/NC NG/NG NG/NT NG/NN NG/N- 
      8888.0, 0, 8888.0, 8888.0, 8888.0, 8888.0, // NG/-A NG/-C NG/-G NG/-T NG/-N NG/-- 
      -7900, 1000.0, -1300.0, 200.0, 1000.0, 1000.0, // NT/AA NT/AC NT/AG NT/AT NT/AN NT/A- 
      -8400, 5200.0, -4400.0, -2200.0, 5200.0, 5200.0, // NT/CA NT/CC NT/CG NT/CT NT/CN NT/C- 
      -7800, -1500.0, -2800.0, -5000.0, 0, 0, // NT/GA NT/GC NT/GG NT/GT NT/GN NT/G- 
      -7200, -1200.0, -4000.0, -2700.0, 0, 0, // NT/TA NT/TC NT/TG NT/TT NT/TN NT/T- 
      0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // NT/NA NT/NC NT/NG NT/NT NT/NN NT/N- 
      0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // NT/-A NT/-C NT/-G NT/-T NT/-N NT/-- 
      -7900, -8500, -8200, -7200, 1000.0, 1000.0, // NN/AA NN/AC NN/AG NN/AT NN/AN NN/A- 
      -8400, -8000, -9800, -8200, 3300.0, 3300.0, // NN/CA NN/CC NN/CG NN/CT NN/CN NN/C- 
      -7800, -10600, -8000, -8500, 0, 0, // NN/GA NN/GC NN/GG NN/GT NN/GN NN/G- 
      -7200, -7800, -8400, -7900, 0, 0, // NN/TA NN/TC NN/TG NN/TT NN/TN NN/T- 
      0, 0, 0, 1200.0, 8888.0, 8888.0, // NN/NA NN/NC NN/NG NN/NT NN/NN NN/N- 
      0, 0, 0, 1200.0, 8888.0, 8888.0, // NN/-A NN/-C NN/-G NN/-T NN/-N NN/-- 
      7600.0, 6100.0, 1600.0, 1200.0, 7600.0, 7600.0, // N-/AA N-/AC N-/AG N-/AT N-/AN N-/A- 
      500.0, 5200.0, -600.0, 3300.0, 5200.0, 5200.0, // N-/CA N-/CC N-/CG N-/CT N-/CN N-/C- 
      600.0, 1900.0, -700.0, -800.0, 1900.0, 1900.0, // N-/GA N-/GC N-/GG N-/GT N-/GN N-/G- 
      5300.0, 2300.0, -600.0, 1000.0, 5300.0, 5300.0, // N-/TA N-/TC N-/TG N-/TT N-/TN N-/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // N-/NA N-/NC N-/NG N-/NT N-/NN N-/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // N-/-A N-/-C N-/-G N-/-T N-/-N N-/-- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -A/AA -A/AC -A/AG -A/AT -A/AN -A/A- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -A/CA -A/CC -A/CG -A/CT -A/CN -A/C- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -A/GA -A/GC -A/GG -A/GT -A/GN -A/G- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -A/TA -A/TC -A/TG -A/TT -A/TN -A/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -A/NA -A/NC -A/NG -A/NT -A/NN -A/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -A/-A -A/-C -A/-G -A/-T -A/-N -A/-- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -C/AA -C/AC -C/AG -C/AT -C/AN -C/A- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -C/CA -C/CC -C/CG -C/CT -C/CN -C/C- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -C/GA -C/GC -C/GG -C/GT -C/GN -C/G- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -C/TA -C/TC -C/TG -C/TT -C/TN -C/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -C/NA -C/NC -C/NG -C/NT -C/NN -C/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -C/-A -C/-C -C/-G -C/-T -C/-N -C/-- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -G/AA -G/AC -G/AG -G/AT -G/AN -G/A- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -G/CA -G/CC -G/CG -G/CT -G/CN -G/C- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -G/GA -G/GC -G/GG -G/GT -G/GN -G/G- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -G/TA -G/TC -G/TG -G/TT -G/TN -G/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -G/NA -G/NC -G/NG -G/NT -G/NN -G/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -G/-A -G/-C -G/-G -G/-T -G/-N -G/-- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -T/AA -T/AC -T/AG -T/AT -T/AN -T/A- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -T/CA -T/CC -T/CG -T/CT -T/CN -T/C- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -T/GA -T/GC -T/GG -T/GT -T/GN -T/G- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -T/TA -T/TC -T/TG -T/TT -T/TN -T/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -T/NA -T/NC -T/NG -T/NT -T/NN -T/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -T/-A -T/-C -T/-G -T/-T -T/-N -T/-- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -N/AA -N/AC -N/AG -N/AT -N/AN -N/A- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -N/CA -N/CC -N/CG -N/CT -N/CN -N/C- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -N/GA -N/GC -N/GG -N/GT -N/GN -N/G- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -N/TA -N/TC -N/TG -N/TT -N/TN -N/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -N/NA -N/NC -N/NG -N/NT -N/NN -N/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // -N/-A -N/-C -N/-G -N/-T -N/-N -N/-- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // --/AA --/AC --/AG --/AT --/AN --/A- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // --/CA --/CC --/CG --/CT --/CN --/C- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // --/GA --/GC --/GG --/GT --/GN --/G- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // --/TA --/TC --/TG --/TT --/TN --/T- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // --/NA --/NC --/NG --/NT --/NN --/N- 
      8888.0, 8888.0, 8888.0, 8888.0, 8888.0, 8888.0, // --/-A --/-C --/-G --/-T --/-N --/-- 
    ]


    this.dS_mm = [
      40.0, 40.0, 40.0, 12.9, 40.0, 40.0, // AA/AA AA/AC AA/AG AA/AT AA/AN AA/A- 
      40.0, 40.0, 40.0, 20.2, 40.0, 40.0, // AA/CA AA/CC AA/CG AA/CT AA/CN AA/C- 
      40.0, 40.0, 40.0, 7.4, 40.0, 40.0, // AA/GA AA/GC AA/GG AA/GT AA/GN AA/G- 
      1.7, 4.6, -2.3, -22.2, 4.6, 4.6, // AA/TA AA/TC AA/TG AA/TT AA/TN AA/T- 
      40.0, 40.0, 40.0, 20.2, 40.0, 40.0, // AA/NA AA/NC AA/NG AA/NT AA/NN AA/N- 
      40.0, 40.0, 40.0, 20.2, 40.0, 40.0, // AA/-A AA/-C AA/-G AA/-T AA/-N AA/-- 
      40.0, 40.0, -9.8, 40.0, 40.0, 40.0, // AC/AA AC/AC AC/AG AC/AT AC/AN AC/A- 
      40.0, 40.0, -3.8, 40.0, 40.0, 40.0, // AC/CA AC/CC AC/CG AC/CT AC/CN AC/C- 
      40.0, 40.0, 3.2, 40.0, 40.0, 40.0, // AC/GA AC/GC AC/GG AC/GT AC/GN AC/G- 
      14.6, -4.4, -22.4, 0.2, 14.6, 14.6, // AC/TA AC/TC AC/TG AC/TT AC/TN AC/T- 
      40.0, 40.0, 3.2, 40.0, 40.0, 40.0, // AC/NA AC/NC AC/NG AC/NT AC/NN AC/N- 
      40.0, 40.0, 3.2, 40.0, 40.0, 40.0, // AC/-A AC/-C AC/-G AC/-T AC/-N AC/-- 
      40.0, -4.2, 40.0, 40.0, 40.0, 40.0, // AG/AA AG/AC AG/AG AG/AT AG/AN AG/A- 
      40.0, -0.6, 40.0, 40.0, 40.0, 40.0, // AG/CA AG/CC AG/CG AG/CT AG/CN AG/C- 
      40.0, -13.2, 40.0, 40.0, 40.0, 40.0, // AG/GA AG/GC AG/GG AG/GT AG/GN AG/G- 
      -2.3, -21.0, -9.5, 0.9, 0.9, 0.9, // AG/TA AG/TC AG/TG AG/TT AG/TN AG/T- 
      40.0, 0, 40.0, 40.0, 40.0, 40.0, // AG/NA AG/NC AG/NG AG/NT AG/NN AG/N- 
      40.0, 0, 40.0, 40.0, 40.0, 40.0, // AG/-A AG/-C AG/-G AG/-T AG/-N AG/-- 
      1.7, 40.0, 40.0, 40.0, 40.0, 40.0, // AT/AA AT/AC AT/AG AT/AT AT/AN AT/A- 
      14.6, 40.0, 40.0, 40.0, 40.0, 40.0, // AT/CA AT/CC AT/CG AT/CT AT/CN AT/C- 
      -2.3, 40.0, 40.0, 40.0, 40.0, 40.0, // AT/GA AT/GC AT/GG AT/GT AT/GN AT/G- 
      -20.4, -6.2, -8.3, -10.8, 0, 0, // AT/TA AT/TC AT/TG AT/TT AT/TN AT/T- 
      14.6, 40.0, 40.0, 40.0, 40.0, 40.0, // AT/NA AT/NC AT/NG AT/NT AT/NN AT/N- 
      14.6, 40.0, 40.0, 40.0, 40.0, 40.0, // AT/-A AT/-C AT/-G AT/-T AT/-N AT/-- 
      1.7, -4.2, -9.8, 12.9, 40.0, 40.0, // AN/AA AN/AC AN/AG AN/AT AN/AN AN/A- 
      14.6, -0.6, -3.8, 20.2, 40.0, 40.0, // AN/CA AN/CC AN/CG AN/CT AN/CN AN/C- 
      -2.3, -13.2, 3.2, 7.4, 40.0, 40.0, // AN/GA AN/GC AN/GG AN/GT AN/GN AN/G- 
      -20.4, -21.0, -22.4, -22.2, 0, 0, // AN/TA AN/TC AN/TG AN/TT AN/TN AN/T- 
      14.6, 0, 3.2, 20.2, 40.0, 40.0, // AN/NA AN/NC AN/NG AN/NT AN/NN AN/N- 
      14.6, 0, 3.2, 20.2, 40.0, 40.0, // AN/-A AN/-C AN/-G AN/-T AN/-N AN/-- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // A-/AA A-/AC A-/AG A-/AT A-/AN A-/A- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // A-/CA A-/CC A-/CG A-/CT A-/CN A-/C- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // A-/GA A-/GC A-/GG A-/GT A-/GN A-/G- 
      14.6, 4.6, -2.3, 0.9, 14.6, 14.6, // A-/TA A-/TC A-/TG A-/TT A-/TN A-/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // A-/NA A-/NC A-/NG A-/NT A-/NN A-/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // A-/-A A-/-C A-/-G A-/-T A-/-N A-/-- 
      40.0, 40.0, 40.0, 8.0, 40.0, 40.0, // CA/AA CA/AC CA/AG CA/AT CA/AN CA/A- 
      40.0, 40.0, 40.0, 16.4, 40.0, 40.0, // CA/CA CA/CC CA/CG CA/CT CA/CN CA/C- 
      -4.2, 3.7, -2.3, -22.7, 3.7, 3.7, // CA/GA CA/GC CA/GG CA/GT CA/GN CA/G- 
      40.0, 40.0, 40.0, 0.7, 40.0, 40.0, // CA/TA CA/TC CA/TG CA/TT CA/TN CA/T- 
      40.0, 40.0, 40.0, 16.4, 40.0, 40.0, // CA/NA CA/NC CA/NG CA/NT CA/NN CA/N- 
      40.0, 40.0, 40.0, 16.4, 40.0, 40.0, // CA/-A CA/-C CA/-G CA/-T CA/-N CA/-- 
      40.0, 40.0, 14.2, 40.0, 40.0, 40.0, // CC/AA CC/AC CC/AG CC/AT CC/AN CC/A- 
      40.0, 40.0, 8.9, 40.0, 40.0, 40.0, // CC/CA CC/CC CC/CG CC/CT CC/CN CC/C- 
      -0.6, -7.2, -19.9, -4.5, 0, 0, // CC/GA CC/GC CC/GG CC/GT CC/GN CC/G- 
      40.0, 40.0, 13.5, 40.0, 40.0, 40.0, // CC/TA CC/TC CC/TG CC/TT CC/TN CC/T- 
      40.0, 40.0, 14.2, 40.0, 40.0, 40.0, // CC/NA CC/NC CC/NG CC/NT CC/NN CC/N- 
      40.0, 40.0, 14.2, 40.0, 40.0, 40.0, // CC/-A CC/-C CC/-G CC/-T CC/-N CC/-- 
      40.0, 3.7, 40.0, 40.0, 40.0, 40.0, // CG/AA CG/AC CG/AG CG/AT CG/AN CG/A- 
      40.0, -7.2, 40.0, 40.0, 40.0, 40.0, // CG/CA CG/CC CG/CG CG/CT CG/CN CG/C- 
      -13.2, -27.2, -15.3, -11.7, 0, 0, // CG/GA CG/GC CG/GG CG/GT CG/GN CG/G- 
      40.0, -6.1, 40.0, 40.0, 40.0, 40.0, // CG/TA CG/TC CG/TG CG/TT CG/TN CG/T- 
      40.0, 3.7, 40.0, 40.0, 40.0, 40.0, // CG/NA CG/NC CG/NG CG/NT CG/NN CG/N- 
      40.0, 3.7, 40.0, 40.0, 40.0, 40.0, // CG/-A CG/-C CG/-G CG/-T CG/-N CG/-- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // CT/AA CT/AC CT/AG CT/AT CT/AN CT/A- 
      -4.4, 40.0, 40.0, 40.0, 40.0, 40.0, // CT/CA CT/CC CT/CG CT/CT CT/CN CT/C- 
      -21.0, -6.1, -8.0, -15.8, 0, 0, // CT/GA CT/GC CT/GG CT/GT CT/GN CT/G- 
      -6.2, 40.0, 40.0, 40.0, 40.0, 40.0, // CT/TA CT/TC CT/TG CT/TT CT/TN CT/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // CT/NA CT/NC CT/NG CT/NT CT/NN CT/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // CT/-A CT/-C CT/-G CT/-T CT/-N CT/-- 
      40.0, 3.7, 14.2, 8.0, 40.0, 40.0, // CN/AA CN/AC CN/AG CN/AT CN/AN CN/A- 
      -4.4, -7.2, 8.9, 16.4, 40.0, 40.0, // CN/CA CN/CC CN/CG CN/CT CN/CN CN/C- 
      -21.0, -27.2, -19.9, -22.7, 0, 0, // CN/GA CN/GC CN/GG CN/GT CN/GN CN/G- 
      -6.2, -6.1, 13.5, 0.7, 40.0, 40.0, // CN/TA CN/TC CN/TG CN/TT CN/TN CN/T- 
      40.0, 3.7, 14.2, 16.4, 40.0, 40.0, // CN/NA CN/NC CN/NG CN/NT CN/NN CN/N- 
      40.0, 3.7, 14.2, 16.4, 40.0, 40.0, // CN/-A CN/-C CN/-G CN/-T CN/-N CN/-- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // C-/AA C-/AC C-/AG C-/AT C-/AN C-/A- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // C-/CA C-/CC C-/CG C-/CT C-/CN C-/C- 
      -0.6, 3.7, -2.3, -4.5, 3.7, 3.7, // C-/GA C-/GC C-/GG C-/GT C-/GN C-/G- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // C-/TA C-/TC C-/TG C-/TT C-/TN C-/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // C-/NA C-/NC C-/NG C-/NT C-/NN C-/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // C-/-A C-/-C C-/-G C-/-T C-/-N C-/-- 
      40.0, 40.0, 40.0, 0.7, 40.0, 40.0, // GA/AA GA/AC GA/AG GA/AT GA/AN GA/A- 
      -9.8, 14.2, -1.0, -22.2, 14.2, 14.2, // GA/CA GA/CC GA/CG GA/CT GA/CN GA/C- 
      40.0, 40.0, 40.0, 3.6, 40.0, 40.0, // GA/GA GA/GC GA/GG GA/GT GA/GN GA/G- 
      40.0, 40.0, 40.0, -5.3, 40.0, 40.0, // GA/TA GA/TC GA/TG GA/TT GA/TN GA/T- 
      40.0, 40.0, 40.0, 3.6, 40.0, 40.0, // GA/NA GA/NC GA/NG GA/NT GA/NN GA/N- 
      40.0, 40.0, 40.0, 3.6, 40.0, 40.0, // GA/-A GA/-C GA/-G GA/-T GA/-N GA/-- 
      40.0, 40.0, -1.0, 40.0, 40.0, 40.0, // GC/AA GC/AC GC/AG GC/AT GC/AN GC/A- 
      -3.8, 8.9, -24.4, 5.4, 8.9, 8.9, // GC/CA GC/CC GC/CG GC/CT GC/CN GC/C- 
      40.0, 40.0, -15.8, 40.0, 40.0, 40.0, // GC/GA GC/GC GC/GG GC/GT GC/GN GC/G- 
      40.0, 40.0, -12.3, 40.0, 40.0, 40.0, // GC/TA GC/TC GC/TG GC/TT GC/TN GC/T- 
      40.0, 40.0, 0, 40.0, 40.0, 40.0, // GC/NA GC/NC GC/NG GC/NT GC/NN GC/N- 
      40.0, 40.0, 0, 40.0, 40.0, 40.0, // GC/-A GC/-C GC/-G GC/-T GC/-N GC/-- 
      40.0, -2.3, 40.0, 40.0, 40.0, 40.0, // GG/AA GG/AC GG/AG GG/AT GG/AN GG/A- 
      3.2, -19.9, -15.8, 10.4, 10.4, 10.4, // GG/CA GG/CC GG/CG GG/CT GG/CN GG/C- 
      40.0, -15.3, 40.0, 40.0, 40.0, 40.0, // GG/GA GG/GC GG/GG GG/GT GG/GN GG/G- 
      40.0, -8.0, 40.0, 16.3, 40.0, 40.0, // GG/TA GG/TC GG/TG GG/TT GG/TN GG/T- 
      40.0, 0, 40.0, 40.0, 40.0, 40.0, // GG/NA GG/NC GG/NG GG/NT GG/NN GG/N- 
      40.0, 0, 40.0, 40.0, 40.0, 40.0, // GG/-A GG/-C GG/-G GG/-T GG/-N GG/-- 
      -2.3, 40.0, 40.0, 40.0, 40.0, 40.0, // GT/AA GT/AC GT/AG GT/AT GT/AN GT/A- 
      -22.4, 13.5, -12.3, -8.4, 13.5, 13.5, // GT/CA GT/CC GT/CG GT/CT GT/CN GT/C- 
      -9.5, 40.0, 40.0, 40.0, 40.0, 40.0, // GT/GA GT/GC GT/GG GT/GT GT/GN GT/G- 
      -8.3, 40.0, 9.5, 40.0, 40.0, 40.0, // GT/TA GT/TC GT/TG GT/TT GT/TN GT/T- 
      0, 40.0, 40.0, 40.0, 40.0, 40.0, // GT/NA GT/NC GT/NG GT/NT GT/NN GT/N- 
      0, 40.0, 40.0, 40.0, 40.0, 40.0, // GT/-A GT/-C GT/-G GT/-T GT/-N GT/-- 
      -2.3, -2.3, -1.0, 0.7, 40.0, 40.0, // GN/AA GN/AC GN/AG GN/AT GN/AN GN/A- 
      -22.4, -19.9, -24.4, -22.2, 8.9, 8.9, // GN/CA GN/CC GN/CG GN/CT GN/CN GN/C- 
      -9.5, -15.3, -15.8, 3.6, 40.0, 40.0, // GN/GA GN/GC GN/GG GN/GT GN/GN GN/G- 
      -8.3, -8.0, -12.3, -5.3, 40.0, 40.0, // GN/TA GN/TC GN/TG GN/TT GN/TN GN/T- 
      0, 0, 0, 3.6, 40.0, 40.0, // GN/NA GN/NC GN/NG GN/NT GN/NN GN/N- 
      0, 0, 0, 3.6, 40.0, 40.0, // GN/-A GN/-C GN/-G GN/-T GN/-N GN/-- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // G-/AA G-/AC G-/AG G-/AT G-/AN G-/A- 
      3.2, 14.2, -1.0, 10.4, 14.2, 14.2, // G-/CA G-/CC G-/CG G-/CT G-/CN G-/C- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // G-/GA G-/GC G-/GG G-/GT G-/GN G-/G- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // G-/TA G-/TC G-/TG G-/TT G-/TN G-/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // G-/NA G-/NC G-/NG G-/NT G-/NN G-/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // G-/-A G-/-C G-/-G G-/-T G-/-N G-/-- 
      12.9, 8.0, 0.7, -21.3, 12.9, 12.9, // TA/AA TA/AC TA/AG TA/AT TA/AN TA/A- 
      40.0, 40.0, 40.0, 0.7, 40.0, 40.0, // TA/CA TA/CC TA/CG TA/CT TA/CN TA/C- 
      40.0, 40.0, 40.0, -1.7, 40.0, 40.0, // TA/GA TA/GC TA/GG TA/GT TA/GN TA/G- 
      40.0, 40.0, -1.7, -1.5, 40.0, 40.0, // TA/TA TA/TC TA/TG TA/TT TA/TN TA/T- 
      40.0, 40.0, 40.0, 0.7, 40.0, 40.0, // TA/NA TA/NC TA/NG TA/NT TA/NN TA/N- 
      40.0, 40.0, 40.0, 0.7, 40.0, 40.0, // TA/-A TA/-C TA/-G TA/-T TA/-N TA/-- 
      20.2, 16.4, -22.2, 0.7, 20.2, 20.2, // TC/AA TC/AC TC/AG TC/AT TC/AN TC/A- 
      40.0, 40.0, 5.4, 40.0, 40.0, 40.0, // TC/CA TC/CC TC/CG TC/CT TC/CN TC/C- 
      40.0, 40.0, 10.4, 40.0, 40.0, 40.0, // TC/GA TC/GC TC/GG TC/GT TC/GN TC/G- 
      40.0, 40.0, -8.4, 40.0, 40.0, 40.0, // TC/TA TC/TC TC/TG TC/TT TC/TN TC/T- 
      40.0, 40.0, 10.4, 40.0, 40.0, 40.0, // TC/NA TC/NC TC/NG TC/NT TC/NN TC/N- 
      40.0, 40.0, 10.4, 40.0, 40.0, 40.0, // TC/-A TC/-C TC/-G TC/-T TC/-N TC/-- 
      7.4, -22.7, 3.6, -1.7, 7.4, 7.4, // TG/AA TG/AC TG/AG TG/AT TG/AN TG/A- 
      40.0, -4.5, 40.0, 40.0, 40.0, 40.0, // TG/CA TG/CC TG/CG TG/CT TG/CN TG/C- 
      40.0, -11.7, 40.0, -6.2, 40.0, 40.0, // TG/GA TG/GC TG/GG TG/GT TG/GN TG/G- 
      40.0, -15.8, 40.0, 40.0, 40.0, 40.0, // TG/TA TG/TC TG/TG TG/TT TG/TN TG/T- 
      40.0, 0, 40.0, 40.0, 40.0, 40.0, // TG/NA TG/NC TG/NG TG/NT TG/NN TG/N- 
      40.0, 0, 40.0, 40.0, 40.0, 40.0, // TG/-A TG/-C TG/-G TG/-T TG/-N TG/-- 
      -22.2, 0.7, -5.3, -1.5, 0.7, 0.7, // TT/AA TT/AC TT/AG TT/AT TT/AN TT/A- 
      0.2, 40.0, 40.0, 40.0, 40.0, 40.0, // TT/CA TT/CC TT/CG TT/CT TT/CN TT/C- 
      0.9, 40.0, 16.3, 40.0, 40.0, 40.0, // TT/GA TT/GC TT/GG TT/GT TT/GN TT/G- 
      -10.8, 40.0, 40.0, 40.0, 40.0, 40.0, // TT/TA TT/TC TT/TG TT/TT TT/TN TT/T- 
      0.9, 40.0, 40.0, 40.0, 40.0, 40.0, // TT/NA TT/NC TT/NG TT/NT TT/NN TT/N- 
      0.9, 40.0, 40.0, 40.0, 40.0, 40.0, // TT/-A TT/-C TT/-G TT/-T TT/-N TT/-- 
      -22.2, -22.7, -22.2, -21.3, 0.7, 0.7, // TN/AA TN/AC TN/AG TN/AT TN/AN TN/A- 
      0.2, -4.5, 5.4, 0.7, 40.0, 40.0, // TN/CA TN/CC TN/CG TN/CT TN/CN TN/C- 
      0.9, -11.7, 10.4, -6.2, 40.0, 40.0, // TN/GA TN/GC TN/GG TN/GT TN/GN TN/G- 
      -10.8, -15.8, -8.4, -1.5, 40.0, 40.0, // TN/TA TN/TC TN/TG TN/TT TN/TN TN/T- 
      0.9, 0, 10.4, 0.7, 40.0, 40.0, // TN/NA TN/NC TN/NG TN/NT TN/NN TN/N- 
      0.9, 0, 10.4, 0.7, 40.0, 40.0, // TN/-A TN/-C TN/-G TN/-T TN/-N TN/-- 
      20.2, 16.4, 3.6, 0.7, 20.2, 20.2, // T-/AA T-/AC T-/AG T-/AT T-/AN T-/A- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // T-/CA T-/CC T-/CG T-/CT T-/CN T-/C- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // T-/GA T-/GC T-/GG T-/GT T-/GN T-/G- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // T-/TA T-/TC T-/TG T-/TT T-/TN T-/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // T-/NA T-/NC T-/NG T-/NT T-/NN T-/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // T-/-A T-/-C T-/-G T-/-T T-/-N T-/-- 
      12.9, 8.0, 0.7, -21.3, 12.9, 12.9, // NA/AA NA/AC NA/AG NA/AT NA/AN NA/A- 
      -9.8, 14.2, -1.0, -22.2, 14.2, 14.2, // NA/CA NA/CC NA/CG NA/CT NA/CN NA/C- 
      -4.2, 3.7, -2.3, -22.7, 3.7, 3.7, // NA/GA NA/GC NA/GG NA/GT NA/GN NA/G- 
      1.7, 4.6, -2.3, -22.2, 4.6, 4.6, // NA/TA NA/TC NA/TG NA/TT NA/TN NA/T- 
      40.0, 40.0, 40.0, 0.7, 40.0, 40.0, // NA/NA NA/NC NA/NG NA/NT NA/NN NA/N- 
      40.0, 40.0, 40.0, 0.7, 40.0, 40.0, // NA/-A NA/-C NA/-G NA/-T NA/-N NA/-- 
      20.2, 16.4, -22.2, 0.7, 20.2, 20.2, // NC/AA NC/AC NC/AG NC/AT NC/AN NC/A- 
      -3.8, 8.9, -24.4, 5.4, 8.9, 8.9, // NC/CA NC/CC NC/CG NC/CT NC/CN NC/C- 
      -0.6, -7.2, -19.9, -4.5, 0, 0, // NC/GA NC/GC NC/GG NC/GT NC/GN NC/G- 
      14.6, -4.4, -22.4, 0.2, 14.6, 14.6, // NC/TA NC/TC NC/TG NC/TT NC/TN NC/T- 
      40.0, 40.0, 0, 40.0, 40.0, 40.0, // NC/NA NC/NC NC/NG NC/NT NC/NN NC/N- 
      40.0, 40.0, 0, 40.0, 40.0, 40.0, // NC/-A NC/-C NC/-G NC/-T NC/-N NC/-- 
      7.4, -22.7, 3.6, -1.7, 7.4, 7.4, // NG/AA NG/AC NG/AG NG/AT NG/AN NG/A- 
      3.2, -19.9, -15.8, 10.4, 10.4, 10.4, // NG/CA NG/CC NG/CG NG/CT NG/CN NG/C- 
      -13.2, -27.2, -15.3, -11.7, 0, 0, // NG/GA NG/GC NG/GG NG/GT NG/GN NG/G- 
      -2.3, -21.0, -9.5, 0.9, 0.9, 0.9, // NG/TA NG/TC NG/TG NG/TT NG/TN NG/T- 
      40.0, 0, 40.0, 40.0, 40.0, 40.0, // NG/NA NG/NC NG/NG NG/NT NG/NN NG/N- 
      40.0, 0, 40.0, 40.0, 40.0, 40.0, // NG/-A NG/-C NG/-G NG/-T NG/-N NG/-- 
      -22.2, 0.7, -5.3, -1.5, 0.7, 0.7, // NT/AA NT/AC NT/AG NT/AT NT/AN NT/A- 
      -22.4, 13.5, -12.3, -8.4, 13.5, 13.5, // NT/CA NT/CC NT/CG NT/CT NT/CN NT/C- 
      -21.0, -6.1, -8.0, -15.8, 0, 0, // NT/GA NT/GC NT/GG NT/GT NT/GN NT/G- 
      -20.4, -6.2, -8.3, -10.8, 0, 0, // NT/TA NT/TC NT/TG NT/TT NT/TN NT/T- 
      0, 40.0, 40.0, 40.0, 40.0, 40.0, // NT/NA NT/NC NT/NG NT/NT NT/NN NT/N- 
      0, 40.0, 40.0, 40.0, 40.0, 40.0, // NT/-A NT/-C NT/-G NT/-T NT/-N NT/-- 
      -22.2, -22.7, -22.2, -21.3, 0.7, 0.7, // NN/AA NN/AC NN/AG NN/AT NN/AN NN/A- 
      -22.4, -19.9, -24.4, -22.2, 8.9, 8.9, // NN/CA NN/CC NN/CG NN/CT NN/CN NN/C- 
      -21.0, -27.2, -19.9, -22.7, 0, 0, // NN/GA NN/GC NN/GG NN/GT NN/GN NN/G- 
      -20.4, -21.0, -22.4, -22.2, 0, 0, // NN/TA NN/TC NN/TG NN/TT NN/TN NN/T- 
      0, 0, 0, 0.7, 40.0, 40.0, // NN/NA NN/NC NN/NG NN/NT NN/NN NN/N- 
      0, 0, 0, 0.7, 40.0, 40.0, // NN/-A NN/-C NN/-G NN/-T NN/-N NN/-- 
      20.2, 16.4, 3.6, 0.7, 20.2, 20.2, // N-/AA N-/AC N-/AG N-/AT N-/AN N-/A- 
      3.2, 14.2, -1.0, 10.4, 14.2, 14.2, // N-/CA N-/CC N-/CG N-/CT N-/CN N-/C- 
      -0.6, 3.7, -2.3, -4.5, 3.7, 3.7, // N-/GA N-/GC N-/GG N-/GT N-/GN N-/G- 
      14.6, 4.6, -2.3, 0.9, 14.6, 14.6, // N-/TA N-/TC N-/TG N-/TT N-/TN N-/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // N-/NA N-/NC N-/NG N-/NT N-/NN N-/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // N-/-A N-/-C N-/-G N-/-T N-/-N N-/-- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -A/AA -A/AC -A/AG -A/AT -A/AN -A/A- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -A/CA -A/CC -A/CG -A/CT -A/CN -A/C- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -A/GA -A/GC -A/GG -A/GT -A/GN -A/G- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -A/TA -A/TC -A/TG -A/TT -A/TN -A/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -A/NA -A/NC -A/NG -A/NT -A/NN -A/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -A/-A -A/-C -A/-G -A/-T -A/-N -A/-- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -C/AA -C/AC -C/AG -C/AT -C/AN -C/A- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -C/CA -C/CC -C/CG -C/CT -C/CN -C/C- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -C/GA -C/GC -C/GG -C/GT -C/GN -C/G- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -C/TA -C/TC -C/TG -C/TT -C/TN -C/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -C/NA -C/NC -C/NG -C/NT -C/NN -C/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -C/-A -C/-C -C/-G -C/-T -C/-N -C/-- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -G/AA -G/AC -G/AG -G/AT -G/AN -G/A- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -G/CA -G/CC -G/CG -G/CT -G/CN -G/C- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -G/GA -G/GC -G/GG -G/GT -G/GN -G/G- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -G/TA -G/TC -G/TG -G/TT -G/TN -G/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -G/NA -G/NC -G/NG -G/NT -G/NN -G/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -G/-A -G/-C -G/-G -G/-T -G/-N -G/-- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -T/AA -T/AC -T/AG -T/AT -T/AN -T/A- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -T/CA -T/CC -T/CG -T/CT -T/CN -T/C- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -T/GA -T/GC -T/GG -T/GT -T/GN -T/G- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -T/TA -T/TC -T/TG -T/TT -T/TN -T/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -T/NA -T/NC -T/NG -T/NT -T/NN -T/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -T/-A -T/-C -T/-G -T/-T -T/-N -T/-- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -N/AA -N/AC -N/AG -N/AT -N/AN -N/A- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -N/CA -N/CC -N/CG -N/CT -N/CN -N/C- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -N/GA -N/GC -N/GG -N/GT -N/GN -N/G- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -N/TA -N/TC -N/TG -N/TT -N/TN -N/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -N/NA -N/NC -N/NG -N/NT -N/NN -N/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // -N/-A -N/-C -N/-G -N/-T -N/-N -N/-- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // --/AA --/AC --/AG --/AT --/AN --/A- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // --/CA --/CC --/CG --/CT --/CN --/C- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // --/GA --/GC --/GG --/GT --/GN --/G- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // --/TA --/TC --/TG --/TT --/TN --/T- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // --/NA --/NC --/NG --/NT --/NN --/N- 
      40.0, 40.0, 40.0, 40.0, 40.0, 40.0, // --/-A --/-C --/-G --/-T --/-N --/-- 
    ]





    this.seq = seq;
    this.template = template;
    this.T_KELVIN = 273.15
    this.K_mM = 50
    this.ds = 0
    this.dh = 0
    this.Tm = 0
    this.gapCorrect = 0
    this.base = 4000000000
    this.GC_count = 0

    // primer3 default params  
    this.DNA_nM = 50
    this.dmso_conc = 0
    this.dmso_fact = 0.6
    this.formamide_conc = 0.0
    this.divalent = 1.5
    this.dntp = 0.6
    this.monovalent = 50
  }

  static complement(seq) {
    const trantab = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-' };
    return seq.toUpperCase().split('').map(base => trantab[base] || base).join('');
  }

  static isComplement(seq1, seq2) {
    return Primer.complement(seq1) === seq2.toUpperCase();
  }

  static duplex2idx(duplex) {
    const trantab = { 'A': '0', 'C': '1', 'G': '2', 'T': '3', 'N': '4', '-': '5' };
    return parseInt(duplex.toUpperCase().split('').map(base => trantab[base] || '4').join(''), 6);
  }
  symmetry() {
    const seqLen = this.seq.length;
    if (seqLen % 2 === 1) return 0;
    const mp = Math.floor(seqLen / 2);
    for (let i = 0; i < mp; i++) {
      const s = this.seq[i];
      const e = this.seq[seqLen - i - 1];
      const terminalBaseSet = new Set([s, e]);
      if (!(
        (terminalBaseSet.has('A') && terminalBaseSet.has('T')) ||
        (terminalBaseSet.has('C') && terminalBaseSet.has('G'))
      )) return 0;
    }
    return 1;
  }

  divalentToMonovalent() {
    if (this.divalent === 0) this.dntp = 0;
    if (this.divalent < this.dntp) this.divalent = this.dntp;
    return 120 * Math.sqrt(this.divalent - this.dntp);
  }

  calcTm() {
    const sym = this.symmetry();
    if (sym) {
      this.ds -= 1.4;
      this.base /= 4;
    }
    // terminal A T penalty
    for (let [i, j] of [[this.seq[0], this.template[0]], [this.seq[this.seq.length - 1], this.template[this.template.length - 1]]]) {
      if (['A', 'T'].includes(i) || ['A', 'T'].includes(j)) {
        this.ds += 4.1;
        this.dh += 2300;
      } else {
        this.ds -= 2.8;
        this.dh += 100;
      }
    }
    // calc dH and dS by NN
    for (let i = 1; i < this.seq.length; i++) {
      let twoMerPrimer = this.seq[i - 1] + this.seq[i];
      let twoMerTmp = this.template.slice(i - 1, i + 1);
      const idx = Primer.duplex2idx(twoMerPrimer + twoMerTmp);
      if (this.dH_mm[idx] == 8888.0) {
        let GC_num = 0
        for (let index = 0; index < 2; index++) {
          if (twoMerPrimer[index] != Primer.complement(twoMerTmp[index])) {
            if (twoMerTmp[index] != "-") {
              GC_num = (twoMerTmp[index] == "C" || twoMerTmp[index] == "G")
            } else {
              GC_num = (twoMerPrimer[index] == "C" || twoMerPrimer[index] == "G")
            }
            this.gapCorrect -= 2 * (1 - GC_num) + 3 * GC_num
          }
        }
        continue
      }
      this.dh += this.dH_mm[idx];
      this.ds += this.dS_mm[idx];
    }
    // GC correction
    this.GC_count = this.formamide_conc === 0.0 ? 0 : (this.seq.match(/C/g) || []).length + (this.seq.match(/G/g) || []).length;
    // Salt correction
    this.K_mM += this.divalentToMonovalent();
    this.ds += 0.368 * (this.seq.length - 1) * Math.log(this.K_mM / 1000.0);
    this.Tm = this.dh / (this.ds + 1.987 * Math.log(this.DNA_nM / this.base)) - this.T_KELVIN;
    this.Tm -= this.dmso_conc * this.dmso_fact;
    this.Tm += (0.453 * this.GC_count / this.seq.length - 2.88) * this.formamide_conc;
    // Gap correction
    this.Tm += this.gapCorrect
  }

}


// example
// const primer0 = new Primer('CCC-------------------------------ATTGACGTCAAT','GGGATAACCGCAATGATACCCTTGTATGCAGTAATAACTGCAGTTA')
// primer0.calcTm()
// console.log(primer0.Tm)
// const primer = new Primer('TGTAT----CAT-ATGCCAAGTACGCCC----CCT--ATTGACGTCA','ACATGACCCGTATTACGGTCCGCCCGGTAAATGGCAGTAACTGCAGT')
const primer = new Primer('T-TACTGGGCATAATGCCAGGCGGGCCATTTACCGTCATTGACGTCA', 'ACATGACCCGTATTACGGTCCGCCCGGTAAATGGCAGTAACTGCAGT')
primer.calcTm()
console.log(primer.Tm)
// console.log(primer.dh)

