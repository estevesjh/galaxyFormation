# Build the bias calibs:

obsDate=20211104
bias0=$(printf "%05d" 1)
bias1=$(printf "%05d" 50)

dark0=$(printf "%05d" 1)
dark1=$(printf "%05d" 50)

defect0=$(printf "%05d" 1)
defect1=$(printf "%05d" 50)
defect2=$(printf "%05d" 50)
defect3=$(printf "%05d" 50)

pipetask run -j 8 -d "detector IN (0) AND exposure IN (${obsDate}${bias0}..${obsDate}${bias1}) \
	 AND instrument='LATISS' " \
	 -b /repo/main/butler.yaml -i LATISS/raw/all,LATISS/calib \
	 -o u/jesteves/latiss/bias_${obsDate} \
	 -p $CP_PIPE_DIR/pipelines/Latiss/cpBias.yaml --register-dataset-types

# it didn't run. I got bug
# Trying changing the output 
#export PYTHONPATH=$PYTHONPATH:/usr/lib64/python3.6/site-packages

# Ingest the bias calibs:

butler certify-calibrations /repo/main u/jesteves/latiss/bias_${obsDate} \
       u/jesteves/calib/latiss/bias_${obsDate} --begin-date 1980-01-01 \
       --end-date 2050-01-01 bias

# Build the dark calibs:

pipetask run -j 8 -d "detector IN (0) AND exposure IN (${obsDate}${dark0}..${obsDate}${dark1}) \
	 AND instrument='LATISS' " \
	 -b /repo/main/butler.yaml -i LATISS/raw/all,LATISS/calib \
	 -o u/jesteves/latiss/dark_${obsDate} \
	 -p $CP_PIPE_DIR/pipelines/Latiss/cpDark.yaml \
	 -c isr:doDefect=False --register-dataset-types

# Ingest the dark calibs:

butler certify-calibrations /repo/main u/jesteves/latiss/dark_${obsDate} \
       u/jesteves/calib/latiss/dark_${obsDate} --begin-date 1980-01-01 \
       --end-date 2050-01-01 dark

# Start the chained collection:

butler collection-chain /repo/main u/jesteves/calib/latiss/calib.${obsDate} \
       u/jesteves/calib/latiss/bias_${obsDate} \
       u/jesteves/calib/latiss/dark_${obsDate} 

# Build the defect calibs:

pipetask run -j 8 -d "detector IN (0) AND exposure IN \
	 (${obsDate}${defect0}..${obsDate}${defect1}, ${obsDate}${defect2}..${obsDate}${defect3}) \
	 AND instrument='LATISS' " \
	 -b /repo/main/butler.yaml \
	 -i LATISS/raw/all,LATISS/calib,u/jesteves/calib/latiss/calib.${obsDate} \
	 -o u/jesteves/latiss/defects_${obsDate} \
	 -p $CP_PIPE_DIR/pipelines/Latiss/findDefects.yaml \
	 --register-dataset-types

# Ingest the defect calibs:

butler certify-calibrations /repo/main u/jesteves/latiss/defects_${obsDate} \
       u/jesteves/calib/latiss/defects_${obsDate} --begin-date 1980-01-01 \
       --end-date 2050-01-01 defects

# Add the defects to the chained collection:

butler collection-chain /repo/main --mode=extend \
       u/jesteves/calib/latiss/calib.${obsDate} \
       u/jesteves/calib/latiss/defects_${obsDate}

# Build the flat calibs:

pipetask run -j 8 -d "detector IN (0) AND exposure IN (${obsDate}${flat0}..${obsDate}${flat1}) \
	 AND instrument='LATISS' " \
	 -b /repo/main/butler.yaml \
	 -i LATISS/raw/all,LATISS/calib,u/jesteves/calib/latiss/calib.${obsDate} \
	 -o u/jesteves/latiss/flat_${obsDate} \
	 -p $CP_PIPE_DIR/pipelines/Latiss/cpFlat.yaml --register-dataset-types

# Ingest the flat calibs:

butler certify-calibrations /repo/main u/jesteves/latiss/flat_${obsDate} \
       u/jesteves/calib/latiss/flat_${obsDate} --begin-date 1980-01-01 \
       --end-date 2050-01-01 flat

# Add the flat to the chained collection:

butler collection-chain /repo/main --mode=extend \
       u/jesteves/calib/latiss/calib.${obsDate} \
       u/jesteves/calib/latiss/flat_${obsDate} 

# Run the PTC: 

pipetask run -j 32 -d "detector IN (0) AND instrument='LATISS' AND \
	 exposure IN (${obsDate}00096..${obsDate}00135) AND exposure.observation_type='flat'" \
	 -b /repo/main \
	 -i LATISS/raw/all,LATISS/calib,LATISS/calib,u/jesteves/calib/latiss/calib.${obsDate} \
	 -o u/jesteves/latiss/ptc_${obsDate} \
	 -p /project/jesteves/BOT_LSSTCam/pipelines/measurePhotonTransferCurve.yaml \
	 --register-dataset-types

Plot the PTC:

plotPhotonTransferCurve.py \
	/repo/main/u/jesteves/latiss/ptc_${obsDate}/20220113T011324Z/ptc/ptc_LATISS_RXX_S00_u_cslage_latiss_ptc_${obsDate}_20220113T011324Z.fits \
	--detNum=0 --outDir=/repo/main/u/jesteves/latiss/ptc_${obsDate}/plots
