import os, gdal, sys
import numpy as np
import pandas as pd
import pygwr
import shutil

for year in range(2019,2020):
    for day in range(97,365,8):
        day = "{0:03}".format(day)
        if not os.path.exists("I:\\RemoteSensingData\\AMSR\\4_Downscale\\1_Downscaled\\%s\\D_x_SOILM3_%s%s_dwnscld.tif"%(year, year, day)):
##            print "I:\\RemoteSensingData\\AMSR\\4_Downscale\\1_Downscaled\\%s\\D_x_SOILM3_%s%s_dwnscld.tif"%(year, year, day)
            soilfilenm = "D_x_SOILM3_%s%s"%(year, day)
            soilpath ="I:\\RemoteSensingData\\AMSR\\3_Interpolated\\2_ClippedInterpolatedAt05\\%s\\"%(year)

            ndvifilenm = "MOD09_NDVI.A%s%s"%(year, day)
            ndvipath = "I:\\RemoteSensingData\\ModisReflectance\\MOD09A1\\6_NdviReplaced\\Envelope\\%s\\"%(year)

            lstfilenm = "MOD11A2_%s%s_B5"%(year, day)
            lstpath = "I:\\RemoteSensingData\\ModisLST\\MOD11A2\\5_LstReplaceNoData\\%s\\"%(year)

            if os.path.exists('C:\\Final\\temp_soil'):
                shutil.rmtree('C:\\Final\\temp_soil')

            if os.path.isfile('%s%s.tif' %(soilpath,soilfilenm)):        
                temppath = "C:\\Final\\temp_soil\\"
                if not os.path.exists(temppath):
                    os.makedirs(temppath)

                shppath = "C:\\Final\\temp_soil\\shp\\"
                if not os.path.exists(shppath):
                    os.mkdir(shppath)

                rmsmplpath = "C:\\Final\\temp_soil\\rmsmpl\\"
                if not os.path.exists(rmsmplpath):
                    os.mkdir(rmsmplpath)
                    
                rmsmplpath1 = "C:\\Final\\temp_soil\\rmsmpl1\\"
                if not os.path.exists(rmsmplpath1):
                    os.mkdir(rmsmplpath1)
                    
                csvgwrpath = "C:\\Final\\temp_soil\\csvforgwr\\"
                if not os.path.exists(csvgwrpath):
                    os.mkdir(csvgwrpath)

                gwrpath = "C:\\Final\\temp_soil\\gwr_result\\"
                if not os.path.exists(gwrpath):
                    os.mkdir(gwrpath)

                residual = "C:\\Final\\temp_soil\\residualimg\\"
                if not os.path.exists(residual):
                    os.mkdir(residual)

                param1 = "C:\\Final\\temp_soil\\parameter1\\"
                if not os.path.exists(param1):
                    os.mkdir(param1)

                param2 = "C:\\Final\\temp_soil\\parameter2\\"
                if not os.path.exists(param2):
                    os.mkdir(param2)

                param3 = "C:\\Final\\temp_soil\\parameter3\\"
                if not os.path.exists(param3):
                    os.mkdir(param3)

                output = "I:\\RemoteSensingData\\AMSR\\4_Downscale\\1_Downscaled\\%s\\"%(year)
                if not os.path.exists(output):
                    os.makedirs(output)

                clipped = "C:\\Final\\temp_soil\\clipped\\"
                if not os.path.exists(clipped):
                    os.mkdir(clipped)

                soishp ="C:/Python27/Scripts/dmsScript/SoiShapeFile/MaharashtraShp/MaharashtraSoi_GCS.shp"
                soishp1 ="C:/Python27/Scripts/dmsScript/MaharashtraEnvelopeExtendedAt05/extendedEnvelope.shp"
                
                if not os.path.exists(soishp):
                    os.mkdir(soishp)
                
                print "Converting soil tiff to xyz"
                os.system('gdal_translate -q -of xyz -co ADD_HEADER_LINE=YES -a_nodata -28672 %s%s.tif %s%s.xyz'%(soilpath,  soilfilenm, temppath, soilfilenm ))

                print "Renaming xyz to csv"
                if not os.path.isfile('%s%s.csv' %(temppath,soilfilenm)):
                    os.rename('%s%s.xyz'%(temppath,soilfilenm), '%s%s.csv'%(temppath,soilfilenm))

                print "Converting soil.csv to point soil.shp"
                query = "SELECT X, Y, Z FROM %s WHERE Z != '-28672'"%(soilfilenm)
                os.system('ogr2ogr -f "ESRI Shapefile" -sql "%s" %s%s.shp %s%s.csv -oo X_POSSIBLE_NAMES=X* -oo Y_POSSIBLE_NAMES=Y* -a_srs EPSG:4326 -nlt POINT -skipfailures'%(query,temppath, soilfilenm,temppath, soilfilenm ))

                print "Resampling ndvi to 0.25degree"
                os.system('gdal_translate -a_ullr 72 22.5 81.5 15 %s%s.tif %s%s_R.tif'%(ndvipath, ndvifilenm, rmsmplpath1, ndvifilenm))
                os.system('gdal_translate -tr 0.25 0.25 -a_nodata -28672 %s%s_R.tif %s%s_R.tif'%(rmsmplpath1, ndvifilenm, rmsmplpath, ndvifilenm))
                


                print "Resampling lst to 0.25degree"
                os.system('gdal_translate -a_ullr 72 22.5 81.5 15 %s%s.tif %s%s_R.tif'%(lstpath, lstfilenm, rmsmplpath1, lstfilenm))
                os.system('gdal_translate -tr 0.25 0.25 %s%s_R.tif %s%s_R.tif'%(rmsmplpath1, lstfilenm, rmsmplpath, lstfilenm))

                print "Extracting values ndvi to soil.shp"
                os.system('python extract_values.py %s%s.shp %s%s_R.tif %s%s_R.tif '%(temppath,soilfilenm, rmsmplpath, ndvifilenm, rmsmplpath, lstfilenm))

                print "Convert shp to csv"
                query2 = "SELECT X,Y,Z,MOD09_NDVI,MOD11A2_20 FROM %s WHERE MOD09_NDVI != '-28672' AND MOD11A2_20 != '0' "%(soilfilenm)
                os.system('ogr2ogr -f "CSV" -lco GEOMETRY=AS_XY -overwrite -sql "%s" %s%s.csv %s%s.shp'%(query2,csvgwrpath,soilfilenm,temppath, soilfilenm))
    ##            os.system('ogr2ogr -f "CSV" -lco GEOMETRY=AS_XY -overwrite %s%s.csv %s%s.shp'%(csvgwrpath,soilfilenm,temppath, soilfilenm)) 

    ##            try:
    ##                import pygwr
    ##            except:
    ##                sys.path.append(os.path.abspath('../'))
    ##                import pygwr

                print "Starting..."

                data = pd.read_csv('%s%s.csv'%(csvgwrpath,soilfilenm))

                ndvi = data['MOD09_NDVI'].values
                lst = data['MOD11A2_20'].values
                X = data['X'].values
                Y = data['Y'].values

                ind = [ndvi,lst]

                xy = [X,Y]

                dep = data['Z'].values     # grid_code is the dependent variable

                y = np.array(dep)
                y = np.vstack(dep)

                print "size of y is"
                print len(y)

                x=[]  # independent variables
                for i in range(len(ind[0])):
                    x.append([ind[0][i],ind[1][i]])
                x = np.array(x)

                print "size of ind is"
                print len(x)

                g=[] # locations
                for i in range(len(X)):
                    g.append([X[i],Y[i]])
                g = np.array(g)

                #Create our GWR model
                model = pygwr.GWR(targets=y, samples=x, locations=g)

                #Make the global model first
                print "Estimating global model..."
                globalfit = model.global_regression()
                print "Result for global model:"
                print globalfit.summary()

                #Make the bandwidth selection using simple interval search
                #We use AICc as selection criterion
                print "Estimating optimal bandwidth..."
                bwmin, bwmax, bwstep = 1, 1, 1
                opt_bw, opt_aicc = None, np.inf     # initial values (AICc = infinity)
                for bw in range(bwmin, bwmax+bwstep, bwstep):
                    aicc = model.aicc(bw)   # calculate AICc (and AIC, BIC, deviance and K)
                    print "   Bandwidth: %i -- AICc: %f" % (bw, aicc['aicc'])
                    if aicc['aicc'] < opt_aicc: opt_bw, opt_aicc = bw, aicc['aicc']
                print "   Optimal bandwidth is: %i" % opt_bw

                #Estimate the GWR model at all data points
                print "Estimating GWR model at all data points..."
                gwr_result = model.estimate_at_target_locations(bandwidth=opt_bw)

                #Write the result into a result file
                gwr_result.export_csv('%s%s.csv'%(gwrpath,soilfilenm))

                print "CREate shape file from gwr result csv"
                os.system('ogr2ogr -f "ESRI Shapefile" %s%s.shp %s%s.csv -oo X_POSSIBLE_NAMES=EST_PT_1* -oo Y_POSSIBLE_NAMES=EST_PT_2* -a_srs EPSG:4326 -nlt POINT -skipfailures'%(shppath, soilfilenm,gwrpath,soilfilenm))

                print "add field to shape file"
                os.system('ogrinfo %s%s.shp -sql "ALTER TABLE %s ADD COLUMN residual FLOAT(5,4)"'%(shppath, soilfilenm, soilfilenm))

                print "calculate residual in new column"
                os.system('ogrinfo %s%s.shp -dialect SQLite -sql "UPDATE %s SET residual = TARGET - PREDICT"'%(shppath, soilfilenm,soilfilenm))

                print "idw the residual at 500m"
                os.system('gdal_grid -l %s -zfield residual -a invdist:power=2:smoothing=0.0:radius1=0.0:radius2=0.0:angle=0.0:max_points=12.0:nodata=0 -outsize 2174, 1850 -txe 72.216014 81.29835548648649 -tye 14.926851 22.652061 %s%s.shp %sresidual_%s%s.tif'%( soilfilenm, shppath, soilfilenm, residual, year, day))

                 #create param1 out of shp

                os.system('gdal_grid -zfield PARAM_1 -a nearest -a_srs EPSG:4326 -txe 72 81.5 -tye 15 22.5 -outsize 38 30 %s%s.shp %sparam1_%s%s.tif'%(shppath, soilfilenm,param1, year,day))
                
                os.system('gdal_grid -zfield PARAM_2 -a nearest -a_srs EPSG:4326 -txe 72 81.5 -tye 15 22.5 -outsize 38 30 %s%s.shp %sparam2_%s%s.tif'%(shppath, soilfilenm,param2, year,day))

                os.system('gdal_grid -zfield PARAM_3 -a nearest -a_srs EPSG:4326 -txe 72 81.5 -tye 15 22.5 -outsize 38 30 %s%s.shp %sparam3_%s%s.tif'%(shppath, soilfilenm,param3, year,day))
                
                os.system('gdal_translate -outsize 2174, 1850 -a_nodata -28672 %sparam1_%s%s.tif %sparam1%s%s.tif'%(param1,year, day, rmsmplpath, year, day)) #resample param1 to 500m

                os.system('gdal_translate -outsize 2174, 1850 -a_nodata -28672 %sparam2_%s%s.tif %sparam2%s%s.tif'%(param2,year, day, rmsmplpath, year, day)) #resample param2 to 500m

                os.system('gdal_translate -outsize 2174, 1850 -a_nodata -28672 %sparam3_%s%s.tif %sparam3%s%s.tif'%(param3,year, day, rmsmplpath, year, day)) #resample param3 to 500m
                            
                os.system('gdalwarp -dstnodata -28672 -t_srs "EPSG:4326" -cutline %s -tr 0.0041763689, 0.0041766316 -te 72.645855 15.606164 80.898360 22.034000 -overwrite %sparam1%s%s.tif %sparam1_%s%s.tif'%(soishp,rmsmplpath, year,day,  clipped, year, day)) 
            
                os.system('gdalwarp -dstnodata -28672 -t_srs "EPSG:4326" -cutline %s -tr 0.0041763689, 0.0041766316 -te 72.645855 15.606164 80.898360 22.034000 -overwrite %sparam2%s%s.tif %sparam2_%s%s.tif'%(soishp,rmsmplpath, year,day,  clipped, year, day))
                
                os.system('gdalwarp -dstnodata -28672 -t_srs "EPSG:4326" -cutline %s -tr 0.0041763689, 0.0041766316 -te 72.645855 15.606164 80.898360 22.034000 -overwrite %sparam3%s%s.tif %sparam3_%s%s.tif'%(soishp,rmsmplpath, year,day,  clipped, year, day))
                
                os.system('gdalwarp -dstnodata -28672 -t_srs "EPSG:4326" -cutline %s -te 72.645855 15.606164 80.898360 22.034000 -overwrite %sresidual_%s%s.tif %sresidual_%s%s.tif'%(soishp, residual, year, day,clipped, year, day))
                
                os.system('gdalwarp -dstnodata -28672 -t_srs "EPSG:4326" -cutline %s -te 72.645855 15.606164 80.898360 22.034000 -overwrite %s%s.tif %s%s.tif'%(soishp, ndvipath, ndvifilenm, clipped, ndvifilenm))
                
                os.system('gdalwarp -dstnodata -28672 -t_srs "EPSG:4326" -cutline %s -te 72.645855 15.606164 80.898360 22.034000 -overwrite %s%s.tif %s%s.tif'%(soishp, lstpath, lstfilenm, clipped, lstfilenm))
        
                print "use raster calculation to give output for final downscaled image"
                os.system('python gdal_calc.py -A %sparam1_%s%s.tif -B %sparam2_%s%s.tif -C %sparam3_%s%s.tif -D %s%s.tif -E %s%s.tif -F %sresidual_%s%s.tif --outfile=%s%s_dwnscld.tif --calc="((C + (A * D) + (B * E)) + F)"'%(clipped, year,day,clipped, year,day,clipped,  year,day,clipped, ndvifilenm, clipped, lstfilenm, clipped,year,day, temppath, soilfilenm))

                os.system('gdalwarp -srcnodata 1.79769313486e+308 -dstnodata -28672 -t_srs "EPSG:4326" -cutline %s -tr 0.0041760103, 0.0041760103 -te 72.6458547146811 15.6023908288528 80.8983604630547 22.0339995598726 -overwrite %s%s_dwnscld.tif %s%s_dwnscld.tif' \
                          %(soishp, temppath, soilfilenm, output, soilfilenm ))

print "done"

