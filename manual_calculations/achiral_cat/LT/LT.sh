#!/bin/bash
#SBATCH -N 1
#SBATCH -t 240:00:00
#SBATCH --ntasks-per-node=16
#SBATCH --partition=tc

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"
echo "== Scratch dir. : ${TMPDIR}"


module load ams/2021.102
module list


$AMSBIN/ams <<eor>$SLURM_SUBMIT_DIR/$SLURM_JOB_NAME.out


Task PESScan
PESScan
    ScanCoordinate
        nPoints 15
        Distance 11 2 3 1.54
    End
End
TransitionStateSearch
    ReactionCoordinate
        Distance 11 2 1.0
    End
End
System
    Atoms
        H 0.09039377530115 -1.80988464988162 -0.51555853096461 
        C -0.90781844755927 -0.00771283202526 -0.04943451243937 
        C 0.18266547595065 -0.75219646514015 -0.26512257667393 
        C 1.55928917734665 -0.22013974416922 -0.18291419018735 
        O 1.86635737553854 0.93365884237376 0.07180470331268 
        H -0.82431619256681 1.044912935058 0.20242593978238 
        H 2.35514434118519 -0.97166537611335 -0.378325752864 
        H -1.90537397931056 -0.43378138173538 -0.11863168636987 
        I 0.75455898220339 3.9760638549207 0.75978564902957 
        I -0.28217983426102 6.48163757873695 1.31970909811973 
        C -1.094387020976904 -0.0899331687638257 2.943629466769813 
        H -0.1799355420761556 0.4098111783162233 3.227238493949433 
        H -1.934020332111253 0.563855388716794 3.127956417976025 
        H -1.206614272626499 -0.9928108204201239 3.525537526627025 
    End
    BondOrders
         1 3 1.0
         2 3 2.0
         2 6 1.0
         2 8 1.0
         2 11 1.0
         3 4 1.0
         4 5 2.0
         4 7 1.0
         9 10 1.0
         12 11 1.0
         13 11 1.0
         14 11 1.0
    End
End

Engine ADF
    Basis
        Type TZ2P
        Core None
    End
    SpinPolarization 1.0
    XC
        GGA OLYP
    End
    Unrestricted Yes
    NumericalQuality VeryGood
EndEngine
eor

cp ams.results/ams.rkf $SLURM_SUBMIT_DIR/$SLURM_JOB_NAME.ams.rkf
cp ams.results/adf.rkf $SLURM_SUBMIT_DIR/$SLURM_JOB_NAME.adf.rkf
cp ams.results/ams.log $SLURM_SUBMIT_DIR/$SLURM_JOB_NAME.log

rm -r ams.results
