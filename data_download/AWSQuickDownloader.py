import os
import json
from datetime import datetime
from concurrent import futures
from concurrent.futures import ThreadPoolExecutor
from functools import partial
#from multiprocessing import Pool, RLock, freeze_support
from random import random
from threading import RLock as TRLock
import uuid

from multiprocess import Pool, RLock, freeze_support
import s3fs
import boto3
import botocore
import botocore.config as botoconfig
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map, thread_map

def aws_download_multithread_worker(save_dir: str, 
                         bucket: str,
                         boto3: boto3.resource, 
                         s3_file: str,
                         file_name: str, 
                         progress_pos: int):

    file_size = int(boto3.Object(bucket, s3_file).content_length)
    pretty_file_name = os.path.basename(s3_file)
    file_path = os.path.join(save_dir, file_name)

    try:
        with tqdm(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=pretty_file_name, total=file_size, leave=None, position=progress_pos) as progress:
            boto3.Bucket(bucket).download_file(s3_file, file_path, Callback=progress.update)
            progress.close()    
            progress.display(f"{pretty_file_name}: Finished downloading.", pos=progress_pos)
    except Exception as e:
        print(e)

def aws_download_multithread(save_dir, bucket, keys, file_names):
    '''
    Thin wrapper for multithreaded downloading.
    '''
    boto3_session = boto3.Session()
    boto3_client = boto3_session.resource("s3", config=botoconfig.Config(signature_version=botocore.UNSIGNED))

    tqdm.set_lock(TRLock())
    try:
        with ThreadPoolExecutor(initializer=tqdm.set_lock, initargs=(tqdm.get_lock(),)) as executor:
            executor.map(partial(aws_download_multithread_worker, save_dir, bucket, boto3_client), keys, file_names, range(1, len(keys)+1, 1))
    except Exception as e:
        print(e)

def list_aws_keys(bucket: str, 
                   key_pattern: str,
                   glob_match=False):
    '''
    Function to find valid s3 files based on a bucket name and a key pattern.
    Allows for unix-style glob matching.

    Parameters
    ----------
    bucket  : str
        Bucket name to use.
    key_pattern  : str
        Key pattern string to use. If glob_match is False, must be a
        path to a folder containing files (not subfolders), otherwise
        nothing will be downloaded.
        If glob_match is True, uses standard terminology.
    glob_match  : bool, optional [default: False]
        Turns on glob-style matching for key names.

    Returns
    -------
    out  : A list of valid file keys for the given bucket.
    '''
    s3fs_client = s3fs.S3FileSystem(anon=True)

    if not glob_match:
        return [key.replace(f"{bucket}/", "") for key in s3fs_client.ls(f"{bucket}/{key_pattern}")]
    else:
        return [key.replace(f"{bucket}/", "") for key in s3fs_client.glob(f"{bucket}/{key_pattern}")]

if __name__ == "__main__":        

    bucket = "noaa-goes17"
    keys = []
    file_names = []
    save_dir = '/Users/rpurciel/Documents/Utah Dust Storm/Satellite/G17/'
    for day in range(205, 207, 1):

        for hour in range(0, 25, 1):
            hr_str = str(hour).zfill(2)

            protokey = f"ABI-L2-MCMIPC/2021/{day}/{hr_str}/*"
            files = list_aws_keys(bucket, protokey, glob_match=True)

            for file in files:
                keys += [file]
                file_names += [file[file.rfind('/')+1:]]

    aws_download_multithread(save_dir, bucket, keys, file_names)

