import os
import hashlib
import json

def verify_hash(file):
    file_path = os.path.abspath(file)
    dir_path = os.path.dirname(file_path)
    file_name = os.path.basename(file_path)
    with open(file_path, 'rb') as f:
        hash = hashlib.md5(f.read()).hexdigest()
    saved_hash = None
    # Check wether dir_path/hash.json exists
    if(os.path.exists(dir_path + '/hash.json')):
        saved_hash = json.load(open(dir_path + '/hash.json'))
    else:
        saved_hash = {file_name: hash}
        json.dump(saved_hash, open(dir_path + '/hash.json', 'w'))
        return False

    if(file_name in saved_hash.keys()):
        if(hash == saved_hash[file_name]):
            return True
        else:
            saved_hash[file_name] = hash
            json.dump(saved_hash, open(dir_path + '/hash.json', 'w'))
            return False
    else:
        saved_hash[file_name] = hash
        json.dump(saved_hash, open(dir_path + '/hash.json', 'w'))
        return False

