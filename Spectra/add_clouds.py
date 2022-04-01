import pandas as pd
import numpy as np
import os
from tqdm import tqdm
pd.options.mode.chained_assignment = None  # default='warn'


def add_clouds_to_gcm_output(path, runname, planet_name, grav, MTLX, CLOUDS, MOLEF, aerosol_layers, INITIAL_NTAU):
    column_names = ['lat' , 'lon', 'level' , 'altitude(m)',
                    'pressure(bars)', 'temp(k)',
                    'EW vel(m/s)','NS vel','vert vel']

    df = pd.read_csv(path + planet_name + '.txt', delim_whitespace=True, names=column_names)

    df['tau1'] = 0.
    df['g01']  = 0.
    df['pi01'] = 0.

    df['tau2'] = 0.
    df['g02']  = 0.
    df['pi02'] = 0.

    df['tau3'] = 0.
    df['g03']  = 0.
    df['pi03'] = 0.

    df['tau4'] = 0.
    df['g04']  = 0.
    df['pi04'] = 0.

    df['tau5'] = 0.
    df['g05']  = 0.
    df['pi05'] = 0.
    
    df['tau6'] = 0.
    df['g06']  = 0.
    df['pi06'] = 0.
    
    df['tau7'] = 0.
    df['g07']  = 0.
    df['pi07'] = 0.
    
    df['tau8'] = 0.
    df['g08']  = 0.
    df['pi08'] = 0.
    
    df['tau9'] = 0.
    df['g09']  = 0.
    df['pi09'] = 0.
    
    df['tau10'] = 0.
    df['g010']  = 0.
    df['pi010'] = 0.
    
    df['tau11'] = 0.
    df['g011']  = 0.
    df['pi011'] = 0.
    
    df['tau12'] = 0.
    df['g012']  = 0.
    df['pi012'] = 0.
    
    df['tau13'] = 0.
    df['g013']  = 0.
    df['pi013'] = 0.


    input_pressure_array_cgs               = np.asarray([100.239698379016,132.141657859,174.196630916642,229.635882539686,302.719049572289,399.061426988061,526.065415224947,693.48927854176,914.196914554296,1205.14624297868,1588.69215575258,2094.30413981254,2760.83054489446,3639.48346981948,4797.77360895418,6324.69739007194,8337.57495379246,10991.064366056,14489.0446644606,19100.2807641699,25179.0738256863,33192.4837413202,43756.2153614193,57681.9257539164,76039.5873180106,100239.698379016,132141.657858999,174196.630916641,229635.882539685,302719.049572288,399061.426988059,526065.415224945,693489.278541758,914196.914554293,1205146.24297867,1588692.15575257,2094304.13981253,2760830.54489445,3639483.46981947,4797773.60895417,6324697.39007192,8337574.95379243,10991064.366056,14489044.6644605,19100280.7641698,25179073.8256862,33192483.7413201,43756215.3614192,57681925.7539162,76039587.3180104])
    input_particle_size_array_in_meters    = np.asarray([1.00000000E-07, 1.15139540E-07, 1.32571137E-07,1.52641797E-07, 1.75751062E-07,2.02358965E-07, 2.32995181E-07, 2.68269580E-07, 3.08884360E-07, 3.55648031E-07, 4.09491506E-07,4.71486636E-07, 5.42867544E-07, 6.25055193E-07, 7.19685673E-07, 8.28642773E-07, 9.54095476E-07,1.09854114E-06, 1.26485522E-06, 1.45634848E-06, 1.67683294E-06, 1.93069773E-06, 2.22299648E-06,2.55954792E-06, 2.94705170E-06, 3.39322177E-06, 3.90693994E-06, 4.49843267E-06, 5.17947468E-06,5.96362332E-06, 6.86648845E-06, 7.90604321E-06, 9.10298178E-06, 1.04811313E-05, 1.20679264E-05,1.38949549E-05, 1.59985872E-05, 1.84206997E-05, 2.12095089E-05, 2.44205309E-05, 2.81176870E-05,3.23745754E-05, 3.72759372E-05, 4.29193426E-05, 4.94171336E-05, 5.68986603E-05, 6.55128557E-05,7.54312006E-05, 8.68511374E-05, 1.00000000E-04])
    input_temperature_array                = np.asarray([500.00000000, 551.02040816, 602.04081633, 653.06122449, 704.08163265,755.10204082, 806.12244898, 857.14285714, 908.16326531, 959.18367347, 1010.20408163, 1061.2244898,1112.24489796, 1163.26530612, 1214.28571429, 1265.30612245, 1316.32653061, 1367.34693878, 1418.36734694,1469.3877551, 1520.40816327, 1571.42857143, 1622.44897959, 1673.46938776, 1724.48979592, 1775.51020408,1826.53061224, 1877.55102041, 1928.57142857, 1979.59183673, 2030.6122449, 2081.63265306, 2132.65306122,2183.67346939, 2234.69387755, 2285.71428571, 2336.73469388, 2387.75510204, 2438.7755102, 2489.79591837,2540.81632653, 2591.83673469, 2642.85714286, 2693.87755102, 2744.89795918, 2795.91836735, 2846.93877551,2897.95918367, 2948.97959184, 3000.00000000])
    particle_size_vs_layer_array_in_meters = np.asarray([0.1000E-6, 0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6, 0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1023E-6,0.1060E-6, 0.1108E-6,0.1170E-6,0.1250E-6,0.1360E-6,0.1500E-6,0.1950E-6,0.2285E-6,0.2723E-6,0.3300E-6,0.4060E-6, 0.5060E-6,0.6387E-6,0.8130E-6,1.0430E-6,1.3458E-6,1.7450E-6,2.2710E-6,2.9660E-6,3.8800E-6,5.0870E-6, 6.6767E-6,8.7720E-6,11.536E-6,15.1780E-6,19.9800E-6,26.3100E-6,34.6500E-6,45.6500E-6,60.1540E-6, 79.2700E-6])


    DENSITY = [1.98e3,4.09e3,1.86e3,4.0e3,5.22e3,2.65e3,3.27e3,5.76e3,8.9e3,7.9e3,3.34e3,3.98e3,3.95e3]
    FMOLW   = [31.59,41.30,33.07,36.87,64.40,25.46,59.61,28.37,24.87,23.66,72.99,50.83,43.20]
    CORFACT = [1,0,1.0,1.0,1.0,1.0,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000]

    TconKCl     = [617.032,621.573,626.038,630.552,635.053, 639.555,644.050,648.556,653.049,657.552,662.043,666.609,671.436,676.532,681.879,687.233,692.462,697.665,702.916,708.306,713.767,719.366,725.024,730.775,736.460,742.266,748.065,753.932,759.779,765.571, 771.346,777.201,783.301,789.715,796.379,803.117,809.863,816.737,823.798,831.052,838.426,845.980,853.873,862.074,870.494,878.913,887.351,895.768,904.198,912.867,917.385]
    TconZnS     = [708.296,712.010,716.704,719.507,722.309,727.329,730.718,736.323,739.126,744.731,747.534,753.139,755.942,761.547,764.350,769.955,772.758,778.363,783.969,788.079,792.377,797.982,803.587,806.390, 811.995,815.183,823.206,828.812,834.364,840.022, 845.628,851.233,856.839,862.444,868.049,873.655,879.260,884.865,890.471,896.076,901.682,907.287,915.695,921.897,929.708,936.547,943.722,949.327,955.732,963.063, 966.143]
    TconNa2S    = [776.182,781.957,787.715,793.430,799.135,804.761,810.359, 815.868, 821.297, 826.731, 832.157,837.761,843.767,850.179,856.911,863.688,870.458,877.165,884.006, 890.930,897.948, 905.112,912.216,919.374,926.464,933.659,940.930,948.253,955.902,963.944,972.372,980.941,989.400,998.113,1007.19,1016.45,1025.51,1035.47,1044.52,1054.07,1063.67,1073.49,1083.32,1093.15,1103.16,1113.41,1123.84,1134.53,1145.92,1158.12,1164.38]
    TCONMnS     = [1113.43, 1120.07,1126.64,1133.21,1139.9,1146.46,1153.15,  1159.65, 1166.25, 1173.01, 1180.09, 1187.73,1194.69, 1202.12, 1209.26, 1216.50, 1223.92, 1231.43,1239.19, 1247.17, 1255.39,1263.78, 1272.32,1281.40,1289.69, 1297.76, 1306.00, 1315.24, 1324.18, 1333.27,1342.44, 1351.39, 1360.54, 1369.99, 1379.72, 1389.42,1399.22, 1409.04,  1418.99,1428.77, 1438.60, 1449.11,1459.19, 1469.78, 1481.06, 1492.70, 1504.21, 1515.49,1527.84,1540.17,1545.90]
    TconCr      = [1213.17,1219.05,1224.98,1231.07,1237.15,1243.21,1249.26,1255.35,1261.51,1267.83,1274.22,1280.63,1287.04,1293.56,1300.89,1309.34,1318.81,1329.04,1339.04,1349.79,1361.13,1373.13,1385.33,1397.48,1409.69,1421.78,1434.01,1446.16,1458.55,1471.19,1483.84,1496.49,1508.99,1522.14,1536.54,1552.17,1568.58,1585.09,1601.61,1618.14,1634.62,1651.02,1667.41,1683.80,1700.20,1716.57,1732.89,1749.26,1765.73,1783.15,1792.19]
    TCONSiO2    = [1334.63,1342.58,1350.30,1358.48,1366.64,1374.85,1383.15,1391.59,1400.10,1408.68, 1417.25,1425.87,1434.53, 1443.14, 1451.71,1460.28,1468.90,1477.44,1486.12,1494.77,1503.91,1513.72 ,1524.07, 1534.57,1544.98,1555.25, 1565.35,1575.40, 1585.44,1595.55,1606.04,1617.06, 1628.55,1640.25,1651.98, 1664.07,1676.82,1689.97,1703.53, 1717.16, 1731.36, 1746.01,1761.10,1776.94, 1793.70,1811.36,1829.60, 1848.37,1867.54, 1887.00, 1896.46]
    TCONMg2SiO4 = [1370.00,1378.31,1386.62,1395.03,1403.56,1412.29,1421.17, 1430.18, 1439.17, 1448.16, 1457.24,1466.52,1475.94,1485.57, 1495.23, 1505.09, 1515.04, 1525.21,1535.33,1545.80, 1556.37, 1567.12,1578.02, 1589.13,1600.29, 1611.67,1623.20, 1634.84,1646.61,1658.58,1670.74, 1683.03, 1695.59, 1708.44,1721.29,1734.41,1747.86, 1761.64, 1775.79, 1789.95, 1804.36,1819.11,1834.34, 1850.41, 1867.40, 1885.24,1903.85,1923.11,1943.09, 1963.67, 1974.17]
    TconVO      = [1363.07,1371.53, 1379.98, 1389.36,1399.69,1409.65,1419.40, 1429.94, 1439.11,1450.21, 1458.82,1470.07, 1481.33,1492.58,1503.61,1513.72,1526.34,1536.53,1547.03,1557.31,1569.87,1580.06,1591.07, 1603.14,1616.05,1627.63, 1639.06, 1652.94, 1664.91, 1677.84, 1690.74, 1703.56,1716.58,1729.50, 1742.92,1756.97, 1771.03, 1783.97, 1796.90,1810.39, 1824.44,1838.46, 1852.55, 1867.09, 1880.66, 1897.51,1914.24, 1929.95, 1945.67, 1964.18,1972.03]
    TconNi      = [1315.67, 1323.99,1333.24,1343.43,1353.31,1363.35,1373.31,1383.34,1393.43,1403.31, 1412.86,1421.18,1432.30,1443.43,1454.55,1465.97,1478.76,1491.77,1504.48, 1515.78,1529.72,1540.84,1554.77, 1568.23,1580.96,1593.81,1607.69,1622.14,1635.56,1650.51,1666.23,1681.67,1697.25,1713.65,1728.42,1746.81,1766.68,1786.23,1800.17,1817.54,1838.73,1857.11,1878.31, 1896.69,1917.47,1936.26,1957.45,1973.04,1994.22, 2018.22, 2028.81]
    TCONFe      = [1362.14,1373.43, 1384.60,1395.83, 1407.00,1418.26,1430.00,1442.48, 1455.61, 1469.10, 1482.58,1496.08,1509.58,1523.08,1536.57,1550.07,1563.57,1577.13,1590.57,1604.74,1619.69, 1635.41, 1651.59, 1667.75,1684.03,1700.80, 1718.31, 1736.60, 1755.27, 1773.98,1792.93, 1812.32,1832.10, 1852.28,1872.54, 1892.90,1913.24, 1934.27, 1956.41,1978.37, 2008.05, 2030.80,2051.13, 2081.89, 2103.84, 2132.13,2157.98, 2190.91,2221.92, 2247.77, 2263.48]
    TconCa2SiO4 = [1508.24,1518.14,1527.48,1536.15,1544.81,1556.36,1565.02, 1573.68, 1585.23, 1593.89, 1604.95,1614.10,1625.65,1634.31,1645.86, 1656.04, 1666.07, 1677.62,1689.17, 1700.23,1712.27,1722.50, 1735.37,1746.92,1758.47,1772.69,1784.45,1796.00, 1810.24, 1821.99,1835.26,1847.97,1862.41,1876.85,1891.14,1903.65,1919.06,1934.48, 1949.03,1963.46,1980.70,1995.22,2011.50, 2026.98, 2044.31,2060.60,2078.90,2097.18,2118.35,2139.54, 2150.11]
    TconCaTiO3  = [1600.35,1609.68,1616.74,1626.47,1633.19,1643.94, 1652.50,1662.08,1672.50,1681.30,1692.50,1703.94,1712.95,1723.94,1735.02,1744.61,1755.37,1766.81,1778.25,1788.78,1798.49,1810.72,1821.12,1832.96,1845.32,1855.43,1867.42,1879.89,1892.37,1904.05,1917.31,1929.79,1941.23,1952.67,1966.98,1979.67,1992.73,2007.04,2021.35,2035.31,2047.77,2063.12,2078.59,2096.68,2112.02,2130.14,2144.45,2161.64,2182.03,2204.63, 2213.67]
    TCONAl2O3   = [1685.09, 1693.08, 1700.92,1708.84, 1716.78,1724.79,1732.91, 1741.42,1750.37,1759.75, 1769.31,1779.06,1789.37,1799.94, 1810.65,1820.94,1830.99,1841.08,1851.02, 1861.02,1871.26,1881.74,1892.36,1903.03,1913.69,1924.39,1935.05,1945.95,1957.45,1969.45,1980.72,1989.12,1995.54, 2002.05,2010.76,2021.69,2033.01,2044.36, 2055.67,2066.99,2078.33,2089.65,2100.99,2112.62,2124.88,2137.43,2150.32,2163.28, 2176.89, 2191.32,2198.76]

    qe0  = np.loadtxt('SCATTERING_DATA/KCl_ir_2310_qext.txt')
    g00  = np.loadtxt('SCATTERING_DATA/KCl_ir_2310_gg.txt')
    pi00 = np.loadtxt('SCATTERING_DATA/KCl_ir_2310_pi0.txt')

    qe1  = np.loadtxt('SCATTERING_DATA/ZnS_ir_2310_qext.txt')
    g01  = np.loadtxt('SCATTERING_DATA/ZnS_ir_2310_gg.txt')
    pi01 = np.loadtxt('SCATTERING_DATA/ZnS_ir_2310_pi0.txt')

    qe2  = np.loadtxt('SCATTERING_DATA/Na2S_ir_2310_qext.txt')
    g02  = np.loadtxt('SCATTERING_DATA/Na2S_ir_2310_gg.txt')
    pi02 = np.loadtxt('SCATTERING_DATA/Na2S_ir_2310_pi0.txt')

    qe3  = np.loadtxt('SCATTERING_DATA/MnS_ir_2310_qext.txt')
    g03  = np.loadtxt('SCATTERING_DATA/MnS_ir_2310_gg.txt')
    pi03 = np.loadtxt('SCATTERING_DATA/MnS_ir_2310_pi0.txt')

    qe4  = np.loadtxt('SCATTERING_DATA/Cr_ir_2310_qext.txt')
    g04  = np.loadtxt('SCATTERING_DATA/Cr_ir_2310_gg.txt')
    pi04 = np.loadtxt('SCATTERING_DATA/Cr_ir_2310_pi0.txt')

    qe5  = np.loadtxt('SCATTERING_DATA/SiO2_ir_2310_qext.txt')
    g05  = np.loadtxt('SCATTERING_DATA/SiO2_ir_2310_gg.txt')
    pi05 = np.loadtxt('SCATTERING_DATA/SiO2_ir_2310_pi0.txt')

    qe6  = np.loadtxt('SCATTERING_DATA/Mg2SiO4_ir_2310_qext.txt')
    g06  = np.loadtxt('SCATTERING_DATA/Mg2SiO4_ir_2310_gg.txt')
    pi06 = np.loadtxt('SCATTERING_DATA/Mg2SiO4_ir_2310_pi0.txt')

    qe7  = np.loadtxt('SCATTERING_DATA/VO_ir_2310_qext.txt')
    g07  = np.loadtxt('SCATTERING_DATA/VO_ir_2310_gg.txt')
    pi07 = np.loadtxt('SCATTERING_DATA/VO_ir_2310_pi0.txt')

    qe8  = np.loadtxt('SCATTERING_DATA/Ni_ir_2310_qext.txt')
    g08  = np.loadtxt('SCATTERING_DATA/Ni_ir_2310_gg.txt')
    pi08 = np.loadtxt('SCATTERING_DATA/Ni_ir_2310_pi0.txt')

    qe9  = np.loadtxt('SCATTERING_DATA/Fe_ir_2310_qext.txt')
    g09  = np.loadtxt('SCATTERING_DATA/Fe_ir_2310_gg.txt')
    pi09 = np.loadtxt('SCATTERING_DATA/Fe_ir_2310_pi0.txt')

    qe10  = np.loadtxt('SCATTERING_DATA/CaSiO4_ir_2310_qext.txt')
    g010  = np.loadtxt('SCATTERING_DATA/CaSiO4_ir_2310_gg.txt')
    pi010 = np.loadtxt('SCATTERING_DATA/CaSiO4_ir_2310_pi0.txt')

    qe11  = np.loadtxt('SCATTERING_DATA/CaTiO3_ir_2310_qext.txt')
    g011  = np.loadtxt('SCATTERING_DATA/CaTiO3_ir_2310_gg.txt')
    pi011 = np.loadtxt('SCATTERING_DATA/CaTiO3_ir_2310_pi0.txt')

    qe12  = np.loadtxt('SCATTERING_DATA/Al2O3_ir_2310_qext.txt')
    g012  = np.loadtxt('SCATTERING_DATA/Al2O3_ir_2310_gg.txt')
    pi012 = np.loadtxt('SCATTERING_DATA/Al2O3_ir_2310_pi0.txt')

    QE_OPPR  = np.array([qe0, qe1,qe2,qe3,qe4,qe5,qe6,qe7,qe8,qe9,qe10,qe11,qe12])
    PI0_OPPR = np.array([pi00, pi01,pi02,pi03,pi04,pi05,pi06,pi07,pi08,pi09,pi010,pi011,pi012])
    G0_OPPR  = np.array([g00, g01,g02,g03,g04,g05,g06,g07,g08,g09,g010,g011,g012])

    Tconds = np.array([TconKCl, TconZnS,TconNa2S,TCONMnS,TconCr,TCONSiO2,TCONMg2SiO4, TconVO,TconNi,TCONFe,TconCa2SiO4,TconCaTiO3,TCONAl2O3])

    G = grav

    print ("Adding Clouds, this just fills in 0s unless MOLEF is specified")
    if CLOUDS == 1:
        max_cloud_level1 = 0
        max_cloud_level2 = 0
        max_cloud_level3 = 0
        max_cloud_level4 = 0
        max_cloud_level5 = 0
        max_cloud_level6 = 0
        max_cloud_level7 = 0
        max_cloud_level8 = 0
        max_cloud_level9 = 0
        max_cloud_level10 = 0
        max_cloud_level11 = 0
        max_cloud_level12 = 0
        max_cloud_level13 = 0
        for z in tqdm(range(len(df))):
            i = len(df) - z - 1
            layer_index   = np.abs(input_pressure_array_cgs - df['pressure(bars)'][i]*1e6).argmin()
            temp_loc      = np.abs(input_temperature_array  - df['temp(k)'][i]).argmin()
            particle_size = particle_size_vs_layer_array_in_meters[layer_index]
            size_loc      = np.abs(input_particle_size_array_in_meters  - particle_size).argmin()

            # This is in SI
            dpg           = (df['pressure(bars)'][i]*1e5/G)
            CLOUD_INDEX = 0
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau1'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau1'][i] > 0):
                df['g01'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi01'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau1'][i] != 0):
                max_cloud_level1   = max(layer_index - aerosol_layers, 0, max_cloud_level1)
            if (layer_index <= max_cloud_level1):
                df['tau1'][i] = 0
                df['g01'][i]  = 0
                df['pi01'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level1 = 0


            CLOUD_INDEX = 1
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau2'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau2'][i] > 0):
                df['g02'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi02'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau2'][i] != 0):
                max_cloud_level2   = max(layer_index - aerosol_layers, 0, max_cloud_level2)
            if (layer_index <= max_cloud_level2):
                df['tau2'][i] = 0
                df['g02'][i]  = 0
                df['pi02'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level2 = 0

            CLOUD_INDEX = 2
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau3'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau3'][i] > 0):
                df['g03'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi03'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau3'][i] != 0):
                max_cloud_level3   = max(layer_index - aerosol_layers, 0, max_cloud_level3)
            if (layer_index <= max_cloud_level3):
                df['tau3'][i] = 0
                df['g03'][i]  = 0
                df['pi03'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level3 = 0

            CLOUD_INDEX = 3
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau4'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau4'][i] > 0):
                df['g04'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi04'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau4'][i] != 0):
                max_cloud_level4   = max(layer_index - aerosol_layers, 0, max_cloud_level4)
            if (layer_index <= max_cloud_level4):
                df['tau4'][i] = 0
                df['g04'][i]  = 0
                df['pi04'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level4 = 0

            CLOUD_INDEX = 4
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau5'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau5'][i] > 0):
                df['g05'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi05'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau5'][i] != 0):
                max_cloud_level5   = max(layer_index - aerosol_layers, 0, max_cloud_level5)
            if (layer_index <= max_cloud_level5):
                df['tau5'][i] = 0
                df['g05'][i]  = 0
                df['pi05'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level5 = 0


            CLOUD_INDEX = 5
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau6'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau6'][i] > 0):
                df['g06'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi06'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau6'][i] != 0):
                max_cloud_level6   = max(layer_index - aerosol_layers, 0, max_cloud_level6)
            if (layer_index <= max_cloud_level6):
                df['tau6'][i] = 0
                df['g06'][i]  = 0
                df['pi06'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level6 = 0

            CLOUD_INDEX = 6
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau7'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau7'][i] > 0):
                df['g07'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi07'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau7'][i] > 1e-10):
                max_cloud_level7 = max(layer_index - aerosol_layers, 0.0, max_cloud_level7)
            if (layer_index <= max_cloud_level7):
                df['tau7'][i] = 0
                df['g07'][i]  = 0
                df['pi07'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level7 = 0


            CLOUD_INDEX = 7
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau8'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau8'][i] > 0):
                df['g08'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi08'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau8'][i] != 0):
                max_cloud_level8 = max(layer_index - aerosol_layers, 0, max_cloud_level8)
            if (layer_index <= max_cloud_level8):
                df['tau8'][i] = 0
                df['g08'][i]  = 0
                df['pi08'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level8 = 0

            CLOUD_INDEX = 8
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau9'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau9'][i] > 0):
                df['g09'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi09'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau9'][i] != 0):
                max_cloud_level9   = max(layer_index - aerosol_layers, 0, max_cloud_level9)
            if (layer_index <= max_cloud_level9):
                df['tau9'][i] = 0
                df['g09'][i]  = 0
                df['pi09'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level9 = 0

            CLOUD_INDEX = 9
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau10'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau10'][i] > 0):
                df['g010'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi010'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau10'][i] != 0):
                max_cloud_level10   = max(layer_index - aerosol_layers, 0, max_cloud_level10)
            if (layer_index <= max_cloud_level10):
                df['tau10'][i] = 0
                df['g010'][i]  = 0
                df['pi010'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level10 = 0

            CLOUD_INDEX = 10
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau11'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau11'][i] > 0):
                df['g011'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi011'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau11'][i] != 0):
                max_cloud_level11   = max(layer_index - aerosol_layers, 0, max_cloud_level11)
            if (layer_index <= max_cloud_level11):
                df['tau11'][i] = 0
                df['g011'][i]  = 0
                df['pi011'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level11 = 0

            CLOUD_INDEX = 11
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau12'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau12'][i] > 0):
                df['g012'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi012'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau12'][i] != 0):
                max_cloud_level12   = max(layer_index - aerosol_layers, 0, max_cloud_level12)
            if (layer_index <= max_cloud_level12):
                df['tau12'][i] = 0
                df['g012'][i]  = 0
                df['pi012'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level12 = 0

            CLOUD_INDEX = 12
            cond_fact = (min(max((Tconds[CLOUD_INDEX][layer_index]-df['temp(k)'][i]) / 10.0, 0.0), 1.0))
            df['tau13'][i] = dpg * MOLEF[CLOUD_INDEX]*3./4./particle_size/DENSITY[CLOUD_INDEX]*FMOLW[CLOUD_INDEX]*cond_fact*MTLX*CORFACT[layer_index]*QE_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau13'][i] > 0):
                df['g013'][i]  = G0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
                df['pi013'][i] = PI0_OPPR[CLOUD_INDEX][size_loc][temp_loc]
            if (df['tau13'][i] != 0):
                max_cloud_level13   = max(layer_index - aerosol_layers, 0, max_cloud_level13)

            if (layer_index <= max_cloud_level13):
                df['tau13'][i] = 0
                df['g013'][i]  = 0
                df['pi013'][i] = 0
            if (layer_index == INITIAL_NTAU-1):
                max_cloud_level13 = 0
    else:
        pass

    planet_file_with_clouds = path + planet_name + '_with_clouds.txt'
    np.savetxt(planet_file_with_clouds, df.values,
               fmt=' '.join(['%5.2f']*2 + ['%3d']*1 + ['%9.2E']*6 + ['%9.2E']*39 + ['\t']))
    return None
