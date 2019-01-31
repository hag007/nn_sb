import shutil
import os
import constants

def clear_cache(dataset_name):
    cache_dir=os.path.join(constants.DATASETS_DIR,dataset_name,"cache")
    try:
        shutil.rmtree(cache_dir, ignore_errors=True)
        os.makedirs(cache_dir)
    except:
        print "sth wrong with cache cleaning. continue anyway..."
        pass
