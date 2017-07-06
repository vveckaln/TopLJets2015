#!/bin/bash

WHAT=$1; 
ERA=$2;
if [ "$#" -ne 2 ]; then 
    echo "steerTOPWidthAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/PLOT/WWW>";
    echo "        SHAPES          - create shapes files from plotter files";
    echo "        MERGE_SHAPES    - merge output of shapes generation";
    echo "        MERGE_DATACARDS - merge all datacards for one LFS/WID into one WID datacard";
    echo "        WORKSPACE       - create a workspace with the TopHypoTest model";
    echo "        SCAN            - scans the likelihood";
    echo "        CLs             - calculates CLs using higgs_combine";
    echo "        TOYS            - plots toys and saves locally";
    echo "        QUANTILES       - plots the quantiles chart for all distributions, by lepton final state";
    echo "        WWW             - moves the analysis plots to an afs-web-based area";
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N

queue=2nw
outdir=/afs/cern.ch/work/e/ecoleman/TopWidth/
cardsdir=${outdir}/datacards
wwwdir=~/www/TopWidth_${ERA}/
CMSSW_7_4_7dir=~/CMSSW_7_4_7/src/
CMSSW_7_6_3dir=~/CMSSW_8_0_26_patch1/src/

unblind=false
nPseudo=1000

wid=(0p2w 0p4w 0p6w 0p8w 
    1p0w 1p2w 1p4w 1p6w 
    1p8w 2p0w 2p2w 2p4w 
    2p6w 2p8w 3p0w 3p5w 4p0w)
lbCat=(highpt lowpt)
lfs=(EE EM MM)
cat=(1b 2b)
dists=(incmlb)

nuisances=(jes 
    jesrate jer pu btag les ltag 
    trigEE trigEM trigMM 
    selEE selEM selMM 
    toppt MEqcdscale PDF Herwig amcnloFxFx 
    Mtop ttPartonShower tWttinterf tWQCDScale 
    DYnorm_thEE DYnorm_thEM DYnorm_thMM) 

RED='\e[31m'
NC='\e[0m'

function join { local IFS=','; echo "$*"; }

# Helpers: getting comma-separated lists of variables
distStr="$(join ${dists[@]})"
lfsStr="$(join ${lfs[@]})"
widStr="$(join ${wid[@]})"

case $WHAT in
############################### SHAPES ####################################
    SHAPES ) # to get the shapes file / datacards for the full analysis
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
    
        rm ${cardsdir}/*.root

        for index in `seq 0 17 204` ; do
            min=$index
            max=$(($index+16))

            nohup python test/TopWidthAnalysis/createShapesFromPlotter.py \
                    -s tbart,Singletop \
                    --dists ${distStr} \
                    -o ${outdir}/datacards/ \
                    -n shapes$index \
                    -i ${outdir}/analysis/plots/plotter.root \
                    --systInput ${outdir}/analysis/plots/syst_plotter.root \
                    --min $min --max $max > shapes${index}.txt & 
        done
    ;;
############################### MERGE_SHAPES ##############################
    MERGE_SHAPES ) # hadd the shapes files since they are split
        hadd ${outdir}/datacards/shapes.root ${outdir}/datacards/shapes*.root 
    ;;
############################### MERGE_DATACARDS ###########################
    MERGE_DATACARDS ) # get all the datacards you could possibly desire for the analysis
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

        # for a given width, merge all
        allcmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
        for tlbCat in ${lbCat[*]} ; do

            # for a given width and lbcat, merge all
            lbccmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
            for tlfs in ${lfs[*]} ; do

                # for a given width, lfs, and lbcat, merge all
                lfscmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
                for tCat in ${cat[*]} ; do
                    cardname="${tlbCat}${tlfs}${tCat}=${cardsdir}"
                    cardname="${cardname}/datacard__${twid}_${tlbCat}${tlfs}${tCat}_${dist}.dat"
                    allcmd="${allcmd} ${cardname} "
                    lbccmd="${lbccmd} ${cardname} "
                    lfscmd="${lfscmd} ${cardname} "
                done

                lfscmd="${lfscmd} > ${cardsdir}/datacard__${twid}_${tlbCat}${tlfs}_${dist}.dat"
            done

            lbccmd="${lbccmd} > ${cardsdir}/datacard__${twid}_${tlbCat}_${dist}.dat"
        done

        allcmd="${allcmd} > ${cardsdir}/datacard__${twid}_${dist}.dat"
            
        eval $allcmd
        eval $lbccmd
        eval $lfscmd

        done
        done
    ;;
############################### WORKSPACE #################################
    WORKSPACE ) # generate combine workspaces
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        mkdir ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do
            
            # All datacards
            echo "Creating workspace for ${twid}${dist}" 
            text2workspace.py ${cardsdir}/datacard__${twid}_${dist}.dat -P \
                HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest \
                -m 172.5 --PO verbose --PO altSignal=${twid} --PO muFloating \
                -o ${outdir}/${twid}_${dist}.root 
        done
        done
    ;;
############################### SCAN ######################################
    SCAN ) # perform likelihood scans on Asimov datasets 
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

            # All datacards
            cmd="combine ${outdir}/${twid}_${dist}.root -M MultiDimFit" 
            cmd="${cmd} -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200"
            cmd="${cmd} --expectSignal=1 --setPhysicsModelParameters x=0,r=1"
            cmd="${cmd} -n x0_scan_Asimov_${twid}_${dist}"
            if [[ ${unblind} == false ]] ; then 
                echo "Analysis is blinded"
                cmd="${cmd} -t -1"
            fi
            
            bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}" 
        done
        done
    ;;
############################### CLs #######################################
    CLs ) # get CLs statistics from combine
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

            # pre-fit 
            echo "Making CLs for ${twid} ${dist}"
            cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
            cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T 2000 -i 2 --fork 8"
            cmd="${cmd} --clsAcc 0 --fullBToys  --saveWorkspace --generateExt=1 --generateNuis=0"
            cmd="${cmd} --expectedFromGrid 0.5 -n x_pre-fit_exp__${twid}_${dist}"
            #cmd="${cmd} &> ${outdir}/x_pre-fit_exp__${twid}_${dist}.log"

            bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${outdir}/" "${cmd}" 

            # post-fit expected 
            echo "Making CLs for ${twid} ${dist}"
            cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
            cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 8"
            cmd="${cmd} --clsAcc 0 --fullBToys  --saveWorkspace --saveToys --frequentist"
            cmd="${cmd} --expectedFromGrid 0.5 -n x_post-fit_exp__${twid}_${dist}"
            #cmd="${cmd} &> ${outdir}/x_post-fit_exp__${twid}_${dist}.log"

            bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${outdir}/" "${cmd}" 

            # post-fit observed
            if [[ ${unblind} == true ]] ; then 
                echo "Making CLs for ${twid} ${dist}"
                cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
                cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 8"
                cmd="${cmd} --clsAcc 0 --fullBToys --saveWorkspace --saveToys --frequentist"
                cmd="${cmd} -n x_post-fit_obs__${twid}_${dist}"
                #cmd="${cmd} &> ${outdir}/x_post-fit_obs__${twid}_${dist}.log"

                bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${outdir}/" "${cmd}" 
            fi
        done
        done
    ;;
############################### TOYS ######################################
    TOYS ) # get toys distributions from the pseudoexperiments
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

            ## pre-fit expected 
            #rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
            #rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_pre-fit_exp__${twid}_${dist}.qvals.root"
            #rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",${nPseudo},\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"pre\\\"\)"

            #cmd=""
            #cmd="${cmd}root -l -q -b"
            #cmd="${cmd} ${outdir}"
            #cmd="${cmd}/higgsCombinex_pre-fit_exp__${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root"
            #cmd="${cmd} ${rootcmds}"

            ##bsub -q ${queue} \
            #sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
            #    "${outdir}/" "${cmd}"

            ## post-fit expected 
            rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
            rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_post-fit_exp__${twid}_${dist}.qvals.root"
            rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",${nPseudo},\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"post\\\"\)"

            cmd=""
            cmd="${cmd}root -l -q -b"
            cmd="${cmd} ${outdir}"
            cmd="${cmd}/higgsCombinex_post-fit_exp__${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root"
            cmd="${cmd} ${rootcmds}"

            #bsub -q ${queue} \
            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}"

            ## post-fit observed 
            if [[ ${unblind} == true ]] ; then 
                rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
                rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_post-fit_obs__${twid}_${dist}.qvals.root"
                rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",${nPseudo},\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"obs\\\"\)"

                cmd=""
                cmd="${cmd}root -l -q -b"
                cmd="${cmd} ${outdir}"
                cmd="${cmd}/higgsCombinex_post-fit_obs__${twid}_${dist}.HybridNew.mH172.5.8192.root"
                cmd="${cmd} ${rootcmds}"

                bsub -q ${queue} \
                    sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                        "${outdir}/" "${cmd}"
            fi
        done
        done
    ;;
############################### WWW #######################################
    WWW )
        mkdir -p ${wwwdir}/ana_${ERA}
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana_${ERA} 
        cp ${outdir}/*.{png,pdf} ${wwwdir}/ana_${ERA} 
        cp test/index.php ${wwwdir}/ana_${ERA}
    ;;
############################### QUANTILES #################################
    QUANTILES ) # plot quantiles distributions of all toys, get CLsPlot
            
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh` 

        cd ${outdir}
        rm statsPlots.root
        for dist in ${dists[*]} ; do

            # Quantiles plot with post-fit information 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --dist ${dist}  \
                --prep pre

            # Quantiles plot with post-fit information 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --dist ${dist}  \
                --prep post #\
                #--unblind

            # Quantiles plot with post-fit information 
            #python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
            #    -i ${outdir}/ -o ${outdir}/ \
            #    --wid ${widStr} \
            #    --dist ${dist}  \
            #    --prep obs #\
            #    #--unblind
            
            # Get CLs plots for pre-fit expectations
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --prep pre \
                --dist ${dist}
            
            # Get CLs plots for post-fit expectations 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --dist ${dist} \
                --prep post \
                --addPre #\
                #--unblind
            
            # Get CLs plots for post-fit expectations 
            #python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
            #    -i ${outdir}/ -o ${outdir}/ \
            #    --wid ${widStr} \
            #    --dist ${dist} \
            #    --prep obs #\
            #    #--unblind

            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getCLsFromFit.py \
                -i ${outdir}/ \
                --dist ${dist} \
                --prep pre

            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getCLsFromFit.py \
                -i ${outdir}/ \
                --dist ${dist} \
                --prep post #\
                #--unblind
        done
        
       # # Get Separation plots for all dists 
       # python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
       #     -i ${outdir}/ -o ${outdir}/ \
       #     --dist ${distStr} \
       #     --doAll \
       #     --unblind
    ;;
############################### PLOT_NUIS #################################
    PLOT_NUIS )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh` 
        cd - 

        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do
            nwid=${twid/p/.}
            nwid=${nwid/w/}

            #python test/TopWidthAnalysis/getNuisances.py \
            #    -i ${outdir}/higgsCombinex_pre-fit_exp__${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root \
            #    -o ${outdir}/ \
            #    -n preNuisances_${twid}_${dist} \
            #    --extraText "#Gamma_{Alt.} = ${nwid} #times #Gamma_{SM}"

            python test/TopWidthAnalysis/getNuisances.py \
                -i ${outdir}/higgsCombinex_post-fit_exp__${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root \
                -o ${outdir}/ \
                -n postNuisances_${twid}_${dist} \
                --extraText "#Gamma_{Alt.} = ${nwid} #times #Gamma_{SM}"
        done
        done
    ;;
############################### PLOT_SYSTS #################################
    PLOT_SYSTS )
        catList="highptEE1b,highptEE2b,lowptEE1b,lowptEE2b,"
        catList="${catList}highptEM1b,highptEM2b,lowptEM1b,lowptEM2b,"
        catList="${catList}highptMM1b,highptMM2b,lowptMM1b,lowptMM2b"

        python test/TopWidthAnalysis/getShapeUncPlots.py \
            -i ${outdir}/withSel/datacards/shapes.root \
            -o ${outdir}/withSel/datacards/plots \
            --uncs btag,ltag,pu \
            --altProc tbart2p0w \
            --cats ${catList} 

        python test/TopWidthAnalysis/getShapeUncPlots.py \
            -i ${outdir}/withSel/datacards/shapes.root \
            -o ${outdir}/withSel/datacards/plots \
            --uncs jes,les,jer \
            --altProc tbart2p0w \
            --cats ${catList} 

        python test/TopWidthAnalysis/getShapeUncPlots.py \
            -i ${outdir}/withSel/datacards/shapes.root \
            -o ${outdir}/withSel/datacards/plots \
            --uncs Mtop \
            --altProc tbart2p0w \
            --cats ${catList} 

        python test/TopWidthAnalysis/getShapeUncPlots.py \
            -i ${outdir}/withSel/datacards/shapes.root \
            -o ${outdir}/withSel/datacards/plots \
            --uncs hdamp,toppt \
            --altProc tbart2p0w \
            --cats ${catList} 

        python test/TopWidthAnalysis/getShapeUncPlots.py \
            -i ${outdir}/withSel/datacards/shapes.root \
            -o ${outdir}/withSel/datacards/plots \
            --uncs ISR,FSR,UE \
            --altProc tbart2p0w \
            --cats ${catList} 
    ;;
############################### NUISANCES #################################
    NUISANCES ) # run once for each nuisance, doing a full analysis to understand effects of systematics 
        for dist in ${dists[*]} ; do
            systList=""
            i=0
        for syst in ${nuisances[*]} ; do
            if [[ "${syst}" != "${nuisances[0]}" ]]; then
               systList="${systList}," 
            fi
            systList="${systList}${syst}"

            echo "Frozen systematics are now: ${systList[@]}"

            mkdir ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}
            cd ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}
            
            # start writing job 
            for twid in ${wid[*]}; do

                # run likelihood scan
                cmd="combine ${outdir}/${twid}_${dist}.root -M MultiDimFit" 
                cmd="${cmd} -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200"
                cmd="${cmd} --expectSignal=1 --setPhysicsModelParameters x=0,r=1"
                cmd="${cmd} -n x0_scan_Asimov_${twid}_${dist}"
                cmd="${cmd} --freezeNuisances ${systList}"
                if [[ ${unblind} == false ]] ; then
                  echo "Analysis is blinded"
                  cmd="${cmd} -t -1"   
                fi
                
                # launch job
                bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}" 

                # run CLs pre-fit
                cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
                cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
                cmd="${cmd} --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0"
                cmd="${cmd} --expectedFromGrid 0.5 -n x_pre-fit_exp_${twid}_${dist}"
                cmd="${cmd} --freezeNuisances ${systList}"
                if [[ ${unblind} == false ]] ; then
                  echo "Analysis is blinded"
                  cmd="${cmd} -t -1"   
                fi
                
                # launch job
                bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}" 

                # run CLs post-fit exp
                cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
                cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
                cmd="${cmd} --clsAcc 0 --fullBToys  --frequentist"
                cmd="${cmd} --expectedFromGrid 0.5 -n x_post-fit_exp_${twid}_${dist}"
                cmd="${cmd} --freezeNuisances ${systList}"
                if [[ ${unblind} == false ]] ; then
                  echo "Analysis is blinded"
                  cmd="${cmd} -t -1"   
                fi
                
                # launch job
                bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}" 

                # run CLs post-fit obs
                cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
                cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
                cmd="${cmd} --clsAcc 0 --fullBToys  --frequentist"
                cmd="${cmd} -n x_post-fit_obs_${twid}_${dist}"
                cmd="${cmd} --freezeNuisances ${systList}"
                if [[ ${unblind} == false ]] ; then
                  echo "Analysis is blinded"
                  cmd="${cmd} -t -1"   
                fi
                
                # launch job
                bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}" 

            done
            
            let "i += 1"
        done    
        done
    ;;
############################### NUISANCES_TOYS ############################
    NUISANCES_TOYS )
        for dist in ${dists[*]} ; do
            i=0
        for syst in ${nuisances[*]} ; do
            cd ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}
            
            # start writing job 
            for twid in ${wid[*]}; do

                ## pre-fit expected 
                rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
                rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_pre-fit_exp_${twid}_${dist}.qvals.root"
                rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"pre\\\"\)"

                cmd=""
                cmd="${cmd}root -l -q -b"
                cmd="${cmd} ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}"
                cmd="${cmd}/higgsCombinex_pre-fit_exp_${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root"
                cmd="${cmd} ${rootcmds}"

                #bsub -q ${queue} \
                sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}" 
            

                ## post-fit expected 
                rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
                rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_post-fit_exp_${twid}_${dist}.qvals.root"
                rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"post\\\"\)"

                cmd=""
                cmd="${cmd}root -l -q -b"
                cmd="${cmd} ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}"
                cmd="${cmd}/higgsCombinex_post-fit_exp_${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root"
                cmd="${cmd} ${rootcmds}"

                #bsub -q ${queue} \
                sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}"
            

                ## post-fit observed 
                rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
                rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_post-fit_obs_${twid}_${dist}.qvals.root"
                rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"obs\\\"\)"

                cmd=""
                cmd="${cmd}root -l -q -b"
                cmd="${cmd} ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}"
                cmd="${cmd}/higgsCombinex_post-fit_obs_${twid}_${dist}.HybridNew.mH172.5.8192.root"
                cmd="${cmd} ${rootcmds}"

                #bsub -q ${queue} \
                sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}"
            
            done

            let "i += 1" 
        done
        done
    ;;
############################### NUISANCES_FITS ############################
    NUISANCE_FITS )
        cd ${CMSSW_7_4_7dir}/
        eval `scramv1 runtime -sh`

        for dist in ${dists[*]} ; do
            i=0
        for syst in ${nuisances[*]} ; do
            dirName=${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}

            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
                -i ${dirName}/ \
                -o ${dirName}/ \
                --dist ${dist} \
                --wid ${widStr} \
                --unblind

            let "i += 1" 

        done
        done

        for dist in ${dists[*]} ; do
            i=0
        for syst in ${nuisances[*]} ; do
            dirName=${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}


            echo " "
            echo "$i ||| $dist $syst : "

            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getCLsFromFit.py \
                -i $dirName \
                --dist ${dist} \
                --unblind
            
            let "i += 1" 

        done
        done
    ;;
esac
