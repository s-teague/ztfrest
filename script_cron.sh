#!/bin/bash

source /scr2/ztfrest/anaconda3/bin/conda activate ztfrest2
export PATH=/scr2/ztfrest/anaconda3/envs/ztfrest2/bin/:$PATH

cd /scr2/ztfrest/ZTF/sarah_test
python filter_kn.py --doWriteDb --doForcePhot --doAuxFp --doCLU

source /scr2/ztfrest/anaconda3/bin/conda activate ztfrest
export PATH=/scr2/ztfrest/anaconda3/envs/ztfrest/bin/:$PATH

python trigger_forcephotztf.py --doWriteDb

source /scr2/ztfrest/anaconda3/bin/conda activate ztfrest2
export PATH=/scr2/ztfrest/anaconda3/envs/ztfrest2/bin/:$PATH

python final_select_fpztf.py --doWriteDb

echo "Finished."
echo "Time:"
echo `date`
