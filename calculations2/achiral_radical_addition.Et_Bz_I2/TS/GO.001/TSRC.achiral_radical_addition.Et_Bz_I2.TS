#!/bin/bash
#SBATCH -N 1
#SBATCH -t 240:00:00
#SBATCH --ntasks-per-node=64
#SBATCH --partition=tc

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"
echo "== Scratch dir. : ${TMPDIR}"


module load ams
module list



echo "========= PREOPTIMIZING ========="
# This is the pre-optimization step
# First use low level of theory (DZP, numerical quality normal)
# Extract the structure using amsreport and save it to a file
# Read in the structure skipping the first row and inputting it in the next run
$AMSBIN/ams <<eor>$SLURM_SUBMIT_DIR/${SLURM_JOB_NAME}_preopt.out

Task TransitionStateSearch 

TransitionStateSearch
  ReactionCoordinate
    Distance 2 10 -1

  End
End

System
  Atoms
         H 	 1.58532161	 2.14308939	-0.10189678
    C 	 3.45079260	 1.12686407	-0.05828611
    C 	 2.09550181	 1.18120576	-0.06109689
    C 	 1.22560214	 0.00000000	 0.00000000
    O 	 0.00000000	 0.00000000	 0.00000000
    H 	 3.78579749	 0.21247831	-0.55308344
    H 	 3.97989967	 2.19315844	 2.34801160
    H 	 4.87732711	 0.57303906	 2.12414324
    H 	 3.03822396	 0.58473702	 2.43521204
    C 	 3.94535932	 1.12219372	 2.19026597
    C 	 1.07263683	-2.46617540	-0.46781104
    C 	 1.97916614	-1.34142182	 0.06579367
    H 	 0.80200244	-2.25514150	-1.49183336
    H 	 0.17904993	-2.52526016	 0.13585399
    H 	 1.60111027	-3.40691278	-0.42167002
    H 	 2.87275303	-1.28233707	-0.53787136
    H 	 2.24980053	-1.55245572	 1.08981599
    C 	 4.23925088	 2.16360119	-0.87993670
    C 	 5.67721306	 1.73717355	-0.97574419
    C 	 6.00881484	 0.42943200	-1.60426084
    C 	 7.28867359	 0.02970764	-1.69965419
    C 	 8.37649407	 0.89114923	-1.17843100
    C 	 8.08401183	 2.07267373	-0.60999190
    C 	 6.67235557	 2.51756094	-0.50276945
    H 	 3.80307952	 2.22423818	-1.89911270
    H 	 4.13415659	 3.15606459	-0.39091823
    H 	 5.21620450	-0.20581793	-1.98783866
    H 	 7.52857262	-0.92423361	-2.15829341
    H 	 9.40874294	 0.56485625	-1.25732544
    H 	 8.87928151	 2.70457092	-0.22787697
    H 	 6.45749935	 3.47611092	-0.03976405
    I 	-0.86240238	 2.87337121	 0.00000000
    I 	-1.65839978	 5.52549284	 0.00000000

  End
End

Engine ADF
  
  BASIS
     type TZ2P
     core none
  END
  
  SYMMETRY NOSYM
  
  SCF
     ITERATIONS 99
     Converge 1.0e-6 
  END
  
  XC
    GGA BLYP
    DISPERSION GRIMME3 BJDAMP
  END
  
  NumericalQuality NORMAL
  
  SpinPolarization 1
  Unrestricted Yes

ENDEngine

eor

cp ams.results/ams.rkf $SLURM_SUBMIT_DIR/${SLURM_JOB_NAME}_preopt.ams.rkf
cp ams.results/adf.rkf $SLURM_SUBMIT_DIR/${SLURM_JOB_NAME}_preopt.adf.rkf
cp ams.results/ams.log $SLURM_SUBMIT_DIR/${SLURM_JOB_NAME}_preopt.log
sleep 5
rm -r ams.results

amsreport $SLURM_SUBMIT_DIR/${SLURM_JOB_NAME}_preopt.adf.rkf geometry-a* > $SLURM_SUBMIT_DIR/preopt.xyz
PREOPT_ATOMS=$(awk -v c=1 'i++ < c {next};1' $SLURM_SUBMIT_DIR/preopt.xyz)


echo "========= OPTIMIZING ========="
# HIGH QUALITY OPTIMIZATION
$AMSBIN/ams <<eor>$SLURM_SUBMIT_DIR/$SLURM_JOB_NAME.out

Task TransitionStateSearch 

TransitionStateSearch
  ReactionCoordinate
        Distance 2 10 -1

  End
End

System
  Atoms
     ${PREOPT_ATOMS}
  End
End

Properties
  NormalModes Yes
End
NormalModes
  ReScanFreqRange -10000000.0 20.0 
End
PESPointCharacter
    NegativeFrequenciesTolerance -30
End

Engine ADF
  
  BASIS
     type TZ2P
     core none
  END
  
  SYMMETRY NOSYM
  
  SCF
     ITERATIONS 99
     Converge 1.0e-6 
  END
  
  XC
    GGA BLYP
    DISPERSION GRIMME3 BJDAMP
  END
  
  NumericalQuality Good
  
  SpinPolarization 1
  Unrestricted Yes

ENDEngine

eor

cp ams.results/ams.rkf $SLURM_SUBMIT_DIR/${SLURM_JOB_NAME}.ams.rkf
cp ams.results/adf.rkf $SLURM_SUBMIT_DIR/${SLURM_JOB_NAME}.adf.rkf
cp ams.results/ams.log $SLURM_SUBMIT_DIR/${SLURM_JOB_NAME}.log
sleep 5
rm -r ams.results

amsreport $SLURM_SUBMIT_DIR/${SLURM_JOB_NAME}.adf.rkf geometry-a* > $SLURM_SUBMIT_DIR/output.xyz