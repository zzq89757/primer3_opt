class Primer {
  constructor(seq, template) {
      // Tables of nearest-neighbor thermodynamics for DNA bases
      
      // Hybridization match table of dH and dS
      this.dH_mm = [
        9999.0,9999.0,9999.0,4700.0,8674.2,  // AA/AA AA/AC AA/AG AA/AT AA/AN 
        9999.0,9999.0,9999.0,7600.0,9399.2,  // AA/CA AA/CC AA/CG AA/CT AA/CN 
        9999.0,9999.0,9999.0,3000.0,8249.2,  // AA/GA AA/GC AA/GG AA/GT AA/GN 
        1200.0,2300.0,-600.0,-7900,-1250.0,  // AA/TA AA/TC AA/TG AA/TT AA/TN 
        7799.2,8074.2,7349.2,1850.0,6268.2,  // AA/NA AA/NC AA/NG AA/NT AA/NN 
        9999.0,9999.0,-2900.0,9999.0,6774.2,  // AC/AA AC/AC AC/AG AC/AT AC/AN 
        9999.0,9999.0,-700.0,9999.0,7324.2,  // AC/CA AC/CC AC/CG AC/CT AC/CN 
        9999.0,9999.0,500.0,9999.0,7624.2,  // AC/GA AC/GC AC/GG AC/GT AC/GN 
        5300.0,0.0,-8400,700.0,-600.0,  // AC/TA AC/TC AC/TG AC/TT AC/TN 
        8824.2,7499.2,-2875.0,7674.2,5280.7,  // AC/NA AC/NC AC/NG AC/NT AC/NN 
        9999.0,-900.0,9999.0,9999.0,7274.2,  // AG/AA AG/AC AG/AG AG/AT AG/AN 
        9999.0,600.0,9999.0,9999.0,7649.2,  // AG/CA AG/CC AG/CG AG/CT AG/CN 
        9999.0,-4000.0,9999.0,9999.0,6499.2,  // AG/GA AG/GC AG/GG AG/GT AG/GN 
        -700.0,-7800,-3100.0,1000.0,-2650.0,  // AG/TA AG/TC AG/TG AG/TT AG/TN 
        7324.2,-3025.0,6724.2,7749.2,4693.2,  // AG/NA AG/NC AG/NG AG/NT AG/NN 
        1200.0,9999.0,9999.0,9999.0,7799.2,  // AT/AA AT/AC AT/AG AT/AT AT/AN 
        5300.0,9999.0,9999.0,9999.0,8824.2,  // AT/CA AT/CC AT/CG AT/CT AT/CN 
        -700.0,9999.0,9999.0,9999.0,7324.2,  // AT/GA AT/GC AT/GG AT/GT AT/GN 
        -7200,-1200.0,-2500.0,-2700.0,-3400.0,  // AT/TA AT/TC AT/TG AT/TT AT/TN 
        -350.0,7199.2,6874.2,6824.2,5136.9,  // AT/NA AT/NC AT/NG AT/NT AT/NN 
        1200.0,-900.0,-2900.0,4700.0,6774.2,  // AN/AA AN/AC AN/AG AN/AT AN/AN 
        5300.0,600.0,-700.0,7600.0,7324.2,  // AN/CA AN/CC AN/CG AN/CT AN/CN 
        -700.0,-4000.0,500.0,3000.0,6499.2,  // AN/GA AN/GC AN/GG AN/GT AN/GN 
        -7200,-7800,-8400,-7900,-3400.0,  // AN/TA AN/TC AN/TG AN/TT AN/TN 
        -350.0,-3025.0,-2875.0,1850.0,4693.2,  // AN/NA AN/NC AN/NG AN/NT AN/NN 
        9999.0,9999.0,9999.0,3400.0,8349.2,  // CA/AA CA/AC CA/AG CA/AT CA/AN 
        9999.0,9999.0,9999.0,6100.0,9024.2,  // CA/CA CA/CC CA/CG CA/CT CA/CN 
        -900.0,1900.0,-700.0,-8500,-2050.0,  // CA/GA CA/GC CA/GG CA/GT CA/GN 
        9999.0,9999.0,9999.0,1000.0,7749.2,  // CA/TA CA/TC CA/TG CA/TT CA/TN 
        7274.2,7974.2,7324.2,500.0,5768.2,  // CA/NA CA/NC CA/NG CA/NT CA/NN 
        9999.0,9999.0,5200.0,9999.0,8799.2,  // CC/AA CC/AC CC/AG CC/AT CC/AN 
        9999.0,9999.0,3600.0,9999.0,8399.2,  // CC/CA CC/CC CC/CG CC/CT CC/CN 
        600.0,-1500.0,-8000,-800.0,-2425.0,  // CC/GA CC/GC CC/GG CC/GT CC/GN 
        9999.0,9999.0,5200.0,9999.0,8799.2,  // CC/TA CC/TC CC/TG CC/TT CC/TN 
        7649.2,7124.2,1500.0,7299.2,5893.2,  // CC/NA CC/NC CC/NG CC/NT CC/NN 
        9999.0,1900.0,9999.0,9999.0,7974.2,  // CG/AA CG/AC CG/AG CG/AT CG/AN 
        9999.0,-1500.0,9999.0,9999.0,7124.2,  // CG/CA CG/CC CG/CG CG/CT CG/CN 
        -4000.0,-10600,-4900.0,-4100.0,-5900.0,  // CG/GA CG/GC CG/GG CG/GT CG/GN 
        9999.0,-1500.0,9999.0,9999.0,7124.2,  // CG/TA CG/TC CG/TG CG/TT CG/TN 
        6499.2,-2925.0,6274.2,6474.2,4080.7,  // CG/NA CG/NC CG/NG CG/NT CG/NN 
        9999.0,9999.0,9999.0,9999.0,9999.0,  // CT/AA CT/AC CT/AG CT/AT CT/AN 
        0.0,9999.0,9999.0,9999.0,7499.2,  // CT/CA CT/CC CT/CG CT/CT CT/CN 
        -7800,-1500.0,-2800.0,-5000.0,-4275.0,  // CT/GA CT/GC CT/GG CT/GT CT/GN 
        -1200.0,9999.0,9999.0,9999.0,7199.2,  // CT/TA CT/TC CT/TG CT/TT CT/TN 
        249.8,7124.2,6799.2,6249.2,5105.6,  // CT/NA CT/NC CT/NG CT/NT CT/NN 
        9999.0,1900.0,5200.0,3400.0,7974.2,  // CN/AA CN/AC CN/AG CN/AT CN/AN 
        0.0,-1500.0,3600.0,6100.0,7124.2,  // CN/CA CN/CC CN/CG CN/CT CN/CN 
        -7800,-10600,-8000,-8500,-5900.0,  // CN/GA CN/GC CN/GG CN/GT CN/GN 
        -1200.0,-1500.0,5200.0,1000.0,7124.2,  // CN/TA CN/TC CN/TG CN/TT CN/TN 
        249.8,-2925.0,1500.0,500.0,4080.7,  // CN/NA CN/NC CN/NG CN/NT CN/NN 
        9999.0,9999.0,9999.0,700.0,7674.2,  // GA/AA GA/AC GA/AG GA/AT GA/AN 
        -2900.0,5200.0,-600.0,-8200,-1625.0,  // GA/CA GA/CC GA/CG GA/CT GA/CN 
        9999.0,9999.0,9999.0,1600.0,7899.2,  // GA/GA GA/GC GA/GG GA/GT GA/GN 
        9999.0,9999.0,9999.0,-1300.0,7174.2,  // GA/TA GA/TC GA/TG GA/TT GA/TN 
        6774.2,8799.2,7349.2,-1800.0,5280.7,  // GA/NA GA/NC GA/NG GA/NT GA/NN 
        9999.0,9999.0,-600.0,9999.0,7349.2,  // GC/AA GC/AC GC/AG GC/AT GC/AN 
        -700.0,3600.0,-9800,2300.0,-1150.0,  // GC/CA GC/CC GC/CG GC/CT GC/CN 
        9999.0,9999.0,-6000.0,9999.0,5999.2,  // GC/GA GC/GC GC/GG GC/GT GC/GN 
        9999.0,9999.0,-4400.0,9999.0,6399.2,  // GC/TA GC/TC GC/TG GC/TT GC/TN 
        7324.2,8399.2,-5200.0,8074.2,4649.4,  // GC/NA GC/NC GC/NG GC/NT GC/NN 
        9999.0,-700.0,9999.0,9999.0,7324.2,  // GG/AA GG/AC GG/AG GG/AT GG/AN 
        500.0,-8000,-6000.0,3300.0,-2550.0,  // GG/CA GG/CC GG/CG GG/CT GG/CN 
        9999.0,-4900.0,9999.0,9999.0,6274.2,  // GG/GA GG/GC GG/GG GG/GT GG/GN 
        9999.0,-2800.0,9999.0,5800.0,5749.5,  // GG/TA GG/TC GG/TG GG/TT GG/TN 
        7624.2,-4100.0,5999.2,7274.5,4199.5,  // GG/NA GG/NC GG/NG GG/NT GG/NN 
        -600.0,9999.0,9999.0,9999.0,7349.2,  // GT/AA GT/AC GT/AG GT/AT GT/AN 
        -8400,5200.0,-4400.0,-2200.0,-2450.0,  // GT/CA GT/CC GT/CG GT/CT GT/CN 
        -3100.0,9999.0,9999.0,9999.0,6724.2,  // GT/GA GT/GC GT/GG GT/GT GT/GN 
        -2500.0,9999.0,4100.0,9999.0,5399.5,  // GT/TA GT/TC GT/TG GT/TT GT/TN 
        -3650.0,8799.2,4924.5,6949.2,4255.8,  // GT/NA GT/NC GT/NG GT/NT GT/NN 
        -600.0,-700.0,-600.0,700.0,7324.2,  // GN/AA GN/AC GN/AG GN/AT GN/AN 
        -8400,-8000,-9800,-8200,-2550.0,  // GN/CA GN/CC GN/CG GN/CT GN/CN 
        -3100.0,-4900.0,-6000.0,1600.0,5999.2,  // GN/GA GN/GC GN/GG GN/GT GN/GN 
        -2500.0,-2800.0,-4400.0,-1300.0,5399.5,  // GN/TA GN/TC GN/TG GN/TT GN/TN 
        -3650.0,-4100.0,-5200.0,-1800.0,4199.5,  // GN/NA GN/NC GN/NG GN/NT GN/NN 
        4700.0,3400.0,700.0,-7200,400.0,  // TA/AA TA/AC TA/AG TA/AT TA/AN 
        9999.0,9999.0,9999.0,1200.0,7799.2,  // TA/CA TA/CC TA/CG TA/CT TA/CN 
        9999.0,9999.0,9999.0,-100.0,7474.2,  // TA/GA TA/GC TA/GG TA/GT TA/GN 
        9999.0,9999.0,-100.0,200.0,5024.5,  // TA/TA TA/TC TA/TG TA/TT TA/TN 
        8674.2,8349.2,5149.5,-1475.0,5174.5,  // TA/NA TA/NC TA/NG TA/NT TA/NN 
        7600.0,6100.0,-8200,1200.0,1675.0,  // TC/AA TC/AC TC/AG TC/AT TC/AN 
        9999.0,9999.0,2300.0,9999.0,8074.2,  // TC/CA TC/CC TC/CG TC/CT TC/CN 
        9999.0,9999.0,3300.0,9999.0,8324.2,  // TC/GA TC/GC TC/GG TC/GT TC/GN 
        9999.0,9999.0,-2200.0,9999.0,6949.2,  // TC/TA TC/TC TC/TG TC/TT TC/TN 
        9399.2,9024.2,-1200.0,7799.2,6255.7,  // TC/NA TC/NC TC/NG TC/NT TC/NN 
        3000.0,-8500,1600.0,-100.0,-1000.0,  // TG/AA TG/AC TG/AG TG/AT TG/AN 
        9999.0,-800.0,9999.0,9999.0,7299.2,  // TG/CA TG/CC TG/CG TG/CT TG/CN 
        9999.0,-4100.0,9999.0,-1400.0,3624.5,  // TG/GA TG/GC TG/GG TG/GT TG/GN 
        9999.0,-5000.0,9999.0,9999.0,6249.2,  // TG/TA TG/TC TG/TG TG/TT TG/TN 
        8249.2,-4600.0,7899.2,4624.5,4043.2,  // TG/NA TG/NC TG/NG TG/NT TG/NN 
        -7900,1000.0,-1300.0,200.0,-2000.0,  // TT/AA TT/AC TT/AG TT/AT TT/AN 
        700.0,9999.0,9999.0,9999.0,7674.2,  // TT/CA TT/CC TT/CG TT/CT TT/CN 
        1000.0,9999.0,5800.0,9999.0,6699.5,  // TT/GA TT/GC TT/GG TT/GT TT/GN 
        -2700.0,9999.0,9999.0,9999.0,6824.2,  // TT/TA TT/TC TT/TG TT/TT TT/TN 
        -2225.0,7749.2,6124.5,7549.2,4799.5,  // TT/NA TT/NC TT/NG TT/NT TT/NN 
        -7900,-8500,-8200,-7200,-2000.0,  // TN/AA TN/AC TN/AG TN/AT TN/AN 
        700.0,-800.0,2300.0,1200.0,7299.2,  // TN/CA TN/CC TN/CG TN/CT TN/CN 
        1000.0,-4100.0,3300.0,-1400.0,3624.5,  // TN/GA TN/GC TN/GG TN/GT TN/GN 
        -2700.0,-5000.0,-2200.0,200.0,5024.5,  // TN/TA TN/TC TN/TG TN/TT TN/TN 
        -2225.0,-4600.0,-1200.0,-1475.0,4043.2,  // TN/NA TN/NC TN/NG TN/NT TN/NN 
        4700.0,3400.0,700.0,-7200,400.0,  // NA/AA NA/AC NA/AG NA/AT NA/AN 
        -2900.0,5200.0,-600.0,-8200,-1625.0,  // NA/CA NA/CC NA/CG NA/CT NA/CN 
        -900.0,1900.0,-700.0,-8500,-2050.0,  // NA/GA NA/GC NA/GG NA/GT NA/GN 
        1200.0,2300.0,-600.0,-7900,-1250.0,  // NA/TA NA/TC NA/TG NA/TT NA/TN 
        6774.2,7974.2,5149.5,-1800.0,5174.5,  // NA/NA NA/NC NA/NG NA/NT NA/NN 
        7600.0,6100.0,-8200,1200.0,1675.0,  // NC/AA NC/AC NC/AG NC/AT NC/AN 
        -700.0,3600.0,-9800,2300.0,-1150.0,  // NC/CA NC/CC NC/CG NC/CT NC/CN 
        600.0,-1500.0,-8000,-800.0,-2425.0,  // NC/GA NC/GC NC/GG NC/GT NC/GN 
        5300.0,0.0,-8400,700.0,-600.0,  // NC/TA NC/TC NC/TG NC/TT NC/TN 
        7324.2,7124.2,-5200.0,7299.2,4649.4,  // NC/NA NC/NC NC/NG NC/NT NC/NN 
        3000.0,-8500,1600.0,-100.0,-1000.0,  // NG/AA NG/AC NG/AG NG/AT NG/AN 
        500.0,-8000,-6000.0,3300.0,-2550.0,  // NG/CA NG/CC NG/CG NG/CT NG/CN 
        -4000.0,-10600,-4900.0,-4100.0,-5900.0,  // NG/GA NG/GC NG/GG NG/GT NG/GN 
        -700.0,-7800,-3100.0,1000.0,-2650.0,  // NG/TA NG/TC NG/TG NG/TT NG/TN 
        6499.2,-4600.0,5999.2,4624.5,4043.2,  // NG/NA NG/NC NG/NG NG/NT NG/NN 
        -7900,1000.0,-1300.0,200.0,-2000.0,  // NT/AA NT/AC NT/AG NT/AT NT/AN 
        -8400,5200.0,-4400.0,-2200.0,-2450.0,  // NT/CA NT/CC NT/CG NT/CT NT/CN 
        -7800,-1500.0,-2800.0,-5000.0,-4275.0,  // NT/GA NT/GC NT/GG NT/GT NT/GN 
        -7200,-1200.0,-2500.0,-2700.0,-3400.0,  // NT/TA NT/TC NT/TG NT/TT NT/TN 
        -3650.0,7124.2,4924.5,6249.2,4255.8,  // NT/NA NT/NC NT/NG NT/NT NT/NN 
        -7900,-8500,-8200,-7200,-2000.0,  // NN/AA NN/AC NN/AG NN/AT NN/AN 
        -8400,-8000,-9800,-8200,-2550.0,  // NN/CA NN/CC NN/CG NN/CT NN/CN 
        -7800,-10600,-8000,-8500,-5900.0,  // NN/GA NN/GC NN/GG NN/GT NN/GN 
        -7200,-7800,-8400,-7900,-3400.0,  // NN/TA NN/TC NN/TG NN/TT NN/TN 
        -3650.0,-4600.0,-5200.0,-1800.0,4043.2,  // NN/NA NN/NC NN/NG NN/NT NN/NN 
        ]

     
        this.dS_mm = [
        25.0,25.0,25.0,12.9,22.0,  // AA/AA AA/AC AA/AG AA/AT AA/AN 
        25.0,25.0,25.0,20.2,23.8,  // AA/CA AA/CC AA/CG AA/CT AA/CN 
        25.0,25.0,25.0,7.4,20.6,  // AA/GA AA/GC AA/GG AA/GT AA/GN 
        1.7,4.6,-2.3,-22.2,-4.6,  // AA/TA AA/TC AA/TG AA/TT AA/TN 
        19.2,19.9,18.2,4.6,15.5,  // AA/NA AA/NC AA/NG AA/NT AA/NN 
        25.0,25.0,-9.8,25.0,16.3,  // AC/AA AC/AC AC/AG AC/AT AC/AN 
        25.0,25.0,-3.8,25.0,17.8,  // AC/CA AC/CC AC/CG AC/CT AC/CN 
        25.0,25.0,3.2,25.0,19.6,  // AC/GA AC/GC AC/GG AC/GT AC/GN 
        14.6,-4.4,-22.4,0.2,-3.0,  // AC/TA AC/TC AC/TG AC/TT AC/TN 
        22.4,17.6,-8.2,18.8,12.7,  // AC/NA AC/NC AC/NG AC/NT AC/NN 
        25.0,-4.2,25.0,25.0,17.7,  // AG/AA AG/AC AG/AG AG/AT AG/AN 
        25.0,-0.6,25.0,25.0,18.6,  // AG/CA AG/CC AG/CG AG/CT AG/CN 
        25.0,-13.2,25.0,25.0,15.4,  // AG/GA AG/GC AG/GG AG/GT AG/GN 
        -2.3,-21.0,-9.5,0.9,-8.0,  // AG/TA AG/TC AG/TG AG/TT AG/TN 
        18.2,-9.8,16.4,19.0,10.9,  // AG/NA AG/NC AG/NG AG/NT AG/NN 
        1.7,25.0,25.0,25.0,19.2,  // AT/AA AT/AC AT/AG AT/AT AT/AN 
        14.6,25.0,25.0,25.0,22.4,  // AT/CA AT/CC AT/CG AT/CT AT/CN 
        -2.3,25.0,25.0,25.0,18.2,  // AT/GA AT/GC AT/GG AT/GT AT/GN 
        -20.4,-6.2,-8.3,-10.8,-11.4,  // AT/TA AT/TC AT/TG AT/TT AT/TN 
        -1.6,17.2,16.7,16.0,12.1,  // AT/NA AT/NC AT/NG AT/NT AT/NN 
        1.7,-4.2,-9.8,12.9,16.3,  // AN/AA AN/AC AN/AG AN/AT AN/AN 
        14.6,-0.6,-3.8,20.2,17.8,  // AN/CA AN/CC AN/CG AN/CT AN/CN 
        -2.3,-13.2,3.2,7.4,15.4,  // AN/GA AN/GC AN/GG AN/GT AN/GN 
        -20.4,-21.0,-22.4,-22.2,-11.4,  // AN/TA AN/TC AN/TG AN/TT AN/TN 
        -1.6,-9.8,-8.2,4.6,10.9,  // AN/NA AN/NC AN/NG AN/NT AN/NN 
        25.0,25.0,25.0,8.0,20.8,  // CA/AA CA/AC CA/AG CA/AT CA/AN 
        25.0,25.0,25.0,16.4,22.8,  // CA/CA CA/CC CA/CG CA/CT CA/CN 
        -4.2,3.7,-2.3,-22.7,-6.4,  // CA/GA CA/GC CA/GG CA/GT CA/GN 
        25.0,25.0,25.0,0.7,18.9,  // CA/TA CA/TC CA/TG CA/TT CA/TN 
        17.7,19.7,18.2,0.6,14.0,  // CA/NA CA/NC CA/NG CA/NT CA/NN 
        25.0,25.0,14.2,25.0,22.3,  // CC/AA CC/AC CC/AG CC/AT CC/AN 
        25.0,25.0,8.9,25.0,21.0,  // CC/CA CC/CC CC/CG CC/CT CC/CN 
        -0.6,-7.2,-19.9,-4.5,-8.0,  // CC/GA CC/GC CC/GG CC/GT CC/GN 
        25.0,25.0,13.5,25.0,22.1,  // CC/TA CC/TC CC/TG CC/TT CC/TN 
        18.6,17.0,4.2,17.6,14.3,  // CC/NA CC/NC CC/NG CC/NT CC/NN 
        25.0,3.7,25.0,25.0,19.7,  // CG/AA CG/AC CG/AG CG/AT CG/AN 
        25.0,-7.2,25.0,25.0,17.0,  // CG/CA CG/CC CG/CG CG/CT CG/CN 
        -13.2,-27.2,-15.3,-11.7,-16.8,  // CG/GA CG/GC CG/GG CG/GT CG/GN 
        25.0,-6.1,25.0,25.0,17.2,  // CG/TA CG/TC CG/TG CG/TT CG/TN 
        15.4,-9.2,14.9,15.8,9.2,  // CG/NA CG/NC CG/NG CG/NT CG/NN 
        25.0,25.0,25.0,25.0,25.0,  // CT/AA CT/AC CT/AG CT/AT CT/AN 
        -4.4,25.0,25.0,25.0,17.6,  // CT/CA CT/CC CT/CG CT/CT CT/CN 
        -21.0,-6.1,-8.0,-15.8,-12.7,  // CT/GA CT/GC CT/GG CT/GT CT/GN 
        -6.2,25.0,25.0,25.0,17.2,  // CT/TA CT/TC CT/TG CT/TT CT/TN 
        -1.6,17.2,16.8,14.8,11.8,  // CT/NA CT/NC CT/NG CT/NT CT/NN 
        25.0,3.7,14.2,8.0,19.7,  // CN/AA CN/AC CN/AG CN/AT CN/AN 
        -4.4,-7.2,8.9,16.4,17.0,  // CN/CA CN/CC CN/CG CN/CT CN/CN 
        -21.0,-27.2,-19.9,-22.7,-16.8,  // CN/GA CN/GC CN/GG CN/GT CN/GN 
        -6.2,-6.1,13.5,0.7,17.2,  // CN/TA CN/TC CN/TG CN/TT CN/TN 
        -1.6,-9.2,4.2,0.6,9.2,  // CN/NA CN/NC CN/NG CN/NT CN/NN 
        25.0,25.0,25.0,0.7,18.9,  // GA/AA GA/AC GA/AG GA/AT GA/AN 
        -9.8,14.2,-1.0,-22.2,-4.7,  // GA/CA GA/CC GA/CG GA/CT GA/CN 
        25.0,25.0,25.0,3.6,19.6,  // GA/GA GA/GC GA/GG GA/GT GA/GN 
        25.0,25.0,25.0,-5.3,17.4,  // GA/TA GA/TC GA/TG GA/TT GA/TN 
        16.3,22.3,18.5,-5.8,12.8,  // GA/NA GA/NC GA/NG GA/NT GA/NN 
        25.0,25.0,-1.0,25.0,18.5,  // GC/AA GC/AC GC/AG GC/AT GC/AN 
        -3.8,8.9,-24.4,5.4,-3.5,  // GC/CA GC/CC GC/CG GC/CT GC/CN 
        25.0,25.0,-15.8,25.0,14.8,  // GC/GA GC/GC GC/GG GC/GT GC/GN 
        25.0,25.0,-12.3,25.0,15.7,  // GC/TA GC/TC GC/TG GC/TT GC/TN 
        17.8,21.0,-13.4,20.1,11.4,  // GC/NA GC/NC GC/NG GC/NT GC/NN 
        25.0,-2.3,25.0,25.0,18.2,  // GG/AA GG/AC GG/AG GG/AT GG/AN 
        3.2,-19.9,-15.8,10.4,-5.5,  // GG/CA GG/CC GG/CG GG/CT GG/CN 
        25.0,-15.3,25.0,25.0,14.9,  // GG/GA GG/GC GG/GG GG/GT GG/GN 
        25.0,-8.0,25.0,16.3,14.6,  // GG/TA GG/TC GG/TG GG/TT GG/TN 
        19.6,-11.4,14.8,19.2,10.5,  // GG/NA GG/NC GG/NG GG/NT GG/NN 
        -2.3,25.0,25.0,25.0,18.2,  // GT/AA GT/AC GT/AG GT/AT GT/AN 
        -22.4,13.5,-12.3,-8.4,-7.4,  // GT/CA GT/CC GT/CG GT/CT GT/CN 
        -9.5,25.0,25.0,25.0,16.4,  // GT/GA GT/GC GT/GG GT/GT GT/GN 
        -8.3,25.0,9.5,25.0,12.8,  // GT/TA GT/TC GT/TG GT/TT GT/TN 
        -10.6,22.1,11.8,16.6,10.0,  // GT/NA GT/NC GT/NG GT/NT GT/NN 
        -2.3,-2.3,-1.0,0.7,18.2,  // GN/AA GN/AC GN/AG GN/AT GN/AN 
        -22.4,-19.9,-24.4,-22.2,-7.4,  // GN/CA GN/CC GN/CG GN/CT GN/CN 
        -9.5,-15.3,-15.8,3.6,14.8,  // GN/GA GN/GC GN/GG GN/GT GN/GN 
        -8.3,-8.0,-12.3,-5.3,12.8,  // GN/TA GN/TC GN/TG GN/TT GN/TN 
        -10.6,-11.4,-13.4,-5.8,10.0,  // GN/NA GN/NC GN/NG GN/NT GN/NN 
        12.9,8.0,0.7,-21.3,0.1,  // TA/AA TA/AC TA/AG TA/AT TA/AN 
        25.0,25.0,25.0,0.7,18.9,  // TA/CA TA/CC TA/CG TA/CT TA/CN 
        25.0,25.0,25.0,-1.7,18.3,  // TA/GA TA/GC TA/GG TA/GT TA/GN 
        25.0,25.0,-1.7,-1.5,11.7,  // TA/TA TA/TC TA/TG TA/TT TA/TN 
        22.0,20.8,12.2,-6.0,12.3,  // TA/NA TA/NC TA/NG TA/NT TA/NN 
        20.2,16.4,-22.2,0.7,3.8,  // TC/AA TC/AC TC/AG TC/AT TC/AN 
        25.0,25.0,5.4,25.0,20.1,  // TC/CA TC/CC TC/CG TC/CT TC/CN 
        25.0,25.0,10.4,25.0,21.4,  // TC/GA TC/GC TC/GG TC/GT TC/GN 
        25.0,25.0,-8.4,25.0,16.6,  // TC/TA TC/TC TC/TG TC/TT TC/TN 
        23.8,22.8,-3.7,18.9,15.5,  // TC/NA TC/NC TC/NG TC/NT TC/NN 
        7.4,-22.7,3.6,-1.7,-3.4,  // TG/AA TG/AC TG/AG TG/AT TG/AN 
        25.0,-4.5,25.0,25.0,17.6,  // TG/CA TG/CC TG/CG TG/CT TG/CN 
        25.0,-11.7,25.0,-6.2,8.0,  // TG/GA TG/GC TG/GG TG/GT TG/GN 
        25.0,-15.8,25.0,25.0,14.8,  // TG/TA TG/TC TG/TG TG/TT TG/TN 
        20.6,-13.7,19.6,10.5,9.3,  // TG/NA TG/NC TG/NG TG/NT TG/NN 
        -22.2,0.7,-5.3,-1.5,-7.1,  // TT/AA TT/AC TT/AG TT/AT TT/AN 
        0.2,25.0,25.0,25.0,18.8,  // TT/CA TT/CC TT/CG TT/CT TT/CN 
        0.9,25.0,16.3,25.0,16.8,  // TT/GA TT/GC TT/GG TT/GT TT/GN 
        -10.8,25.0,25.0,25.0,16.0,  // TT/TA TT/TC TT/TG TT/TT TT/TN 
        -8.0,18.9,15.2,18.4,11.1,  // TT/NA TT/NC TT/NG TT/NT TT/NN 
        -22.2,-22.7,-22.2,-21.3,-7.1,  // TN/AA TN/AC TN/AG TN/AT TN/AN 
        0.2,-4.5,5.4,0.7,17.6,  // TN/CA TN/CC TN/CG TN/CT TN/CN 
        0.9,-11.7,10.4,-6.2,8.0,  // TN/GA TN/GC TN/GG TN/GT TN/GN 
        -10.8,-15.8,-8.4,-1.5,11.7,  // TN/TA TN/TC TN/TG TN/TT TN/TN 
        -8.0,-13.7,-3.7,-6.0,9.3,  // TN/NA TN/NC TN/NG TN/NT TN/NN 
        12.9,8.0,0.7,-21.3,0.1,  // NA/AA NA/AC NA/AG NA/AT NA/AN 
        -9.8,14.2,-1.0,-22.2,-4.7,  // NA/CA NA/CC NA/CG NA/CT NA/CN 
        -4.2,3.7,-2.3,-22.7,-6.4,  // NA/GA NA/GC NA/GG NA/GT NA/GN 
        1.7,4.6,-2.3,-22.2,-4.6,  // NA/TA NA/TC NA/TG NA/TT NA/TN 
        16.3,19.7,12.2,-6.0,12.3,  // NA/NA NA/NC NA/NG NA/NT NA/NN 
        20.2,16.4,-22.2,0.7,3.8,  // NC/AA NC/AC NC/AG NC/AT NC/AN 
        -3.8,8.9,-24.4,5.4,-3.5,  // NC/CA NC/CC NC/CG NC/CT NC/CN 
        -0.6,-7.2,-19.9,-4.5,-8.0,  // NC/GA NC/GC NC/GG NC/GT NC/GN 
        14.6,-4.4,-22.4,0.2,-3.0,  // NC/TA NC/TC NC/TG NC/TT NC/TN 
        17.8,17.0,-13.4,17.6,11.4,  // NC/NA NC/NC NC/NG NC/NT NC/NN 
        7.4,-22.7,3.6,-1.7,-3.4,  // NG/AA NG/AC NG/AG NG/AT NG/AN 
        3.2,-19.9,-15.8,10.4,-5.5,  // NG/CA NG/CC NG/CG NG/CT NG/CN 
        -13.2,-27.2,-15.3,-11.7,-16.8,  // NG/GA NG/GC NG/GG NG/GT NG/GN 
        -2.3,-21.0,-9.5,0.9,-8.0,  // NG/TA NG/TC NG/TG NG/TT NG/TN 
        15.4,-13.7,14.8,10.5,9.2,  // NG/NA NG/NC NG/NG NG/NT NG/NN 
        -22.2,0.7,-5.3,-1.5,-7.1,  // NT/AA NT/AC NT/AG NT/AT NT/AN 
        -22.4,13.5,-12.3,-8.4,-7.4,  // NT/CA NT/CC NT/CG NT/CT NT/CN 
        -21.0,-6.1,-8.0,-15.8,-12.7,  // NT/GA NT/GC NT/GG NT/GT NT/GN 
        -20.4,-6.2,-8.3,-10.8,-11.4,  // NT/TA NT/TC NT/TG NT/TT NT/TN 
        -10.6,17.2,11.8,14.8,10.0,  // NT/NA NT/NC NT/NG NT/NT NT/NN 
        -22.2,-22.7,-22.2,-21.3,-7.1,  // NN/AA NN/AC NN/AG NN/AT NN/AN 
        -22.4,-19.9,-24.4,-22.2,-7.4,  // NN/CA NN/CC NN/CG NN/CT NN/CN 
        -21.0,-27.2,-19.9,-22.7,-16.8,  // NN/GA NN/GC NN/GG NN/GT NN/GN 
        -20.4,-21.0,-22.4,-22.2,-11.4,  // NN/TA NN/TC NN/TG NN/TT NN/TN 
        -10.6,-13.7,-13.4,-6.0,9.2,  // NN/NA NN/NC NN/NG NN/NT NN/NN 
        ]
        this.seq = seq;
        this.template = template;
        this.T_KELVIN = 273.15
        this.K_mM = 50
        this.ds = 0
        this.dh = 0
        this.Tm = 0
        this.base = 4000000000
        this.GC_count = 0
        
        // primer3 default params  
        this.DNA_nM = 50
        this.dmso_conc = 0
        this.dmso_fact = 0.6
        this.formamide_conc = 0.8
        this.divalent = 1.5
        this.dntp = 0.6
        this.monovalent = 50
        
        // cal Tm
        // this.calc_tm()
      }

  static complement(seq) {
    const trantab = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N' };
    return seq.toUpperCase().split('').map(base => trantab[base] || base).join('');
    }

  static isComplement(seq1, seq2) {
      return Primer.complement(seq1) === seq2.toUpperCase();
    }

  static duplex2idx(duplex) {
      const trantab = { 'A': '0', 'C': '1', 'G': '2', 'T': '3', 'N': '4' };
      return parseInt(duplex.toUpperCase().split('').map(base => trantab[base] || base).join(''), 5);
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

      for (let [i, j] of [[this.seq[0], this.template[0]], [this.seq[this.seq.length - 1], this.template[this.template.length - 1]]]) {
          if (['A', 'T'].includes(i) || ['A', 'T'].includes(j)) {
              this.ds += 4.1;
              this.dh += 2300;
          } else {
              this.ds -= 2.8;
              this.dh += 100;
          }
      }

      for (let i = 1; i < this.seq.length; i++) {
          const twoMerPrimer = this.seq[i - 1] + this.seq[i];
          const twoMerTmp = this.template.slice(i - 1, i + 1);
          const idx = Primer.duplex2idx(twoMerPrimer + twoMerTmp);
          this.dh += this.dH_mm[idx];
          this.ds += this.dS_mm[idx];
      }

      this.GC_count = this.formamide_conc === 0.0 ? 0 : (this.seq.match(/C/g) || []).length + (this.seq.match(/G/g) || []).length;
      this.K_mM += this.divalentToMonovalent();
      this.ds += 0.368 * (this.seq.length - 1) * Math.log(this.K_mM / 1000.0);

      this.Tm = this.dh / (this.ds + 1.987 * Math.log(this.DNA_nM / this.base)) - this.T_KELVIN;
      this.Tm -= this.dmso_conc * this.dmso_fact;
      this.Tm += (0.453 * this.GC_count / this.seq.length - 2.88) * this.formamide_conc;
  }
  
    }

    const primer = new Primer("CGATGTGCTAGTTAGTTCNN","GCTACACGATCAATCAAGTA")
    primer.calcTm()
    console.log(primer.Tm)
  