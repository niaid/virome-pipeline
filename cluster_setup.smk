import yaml

############ cluster config ###############
def default_clust(cc):
    """adds __default__ values from cluster config to other
       rules lacking those keys"""
    dkeys = set(cc["__default__"].keys())
    rules = [k for k in cc if k != "__default__"]
    for r in rules:
        toadd = dkeys.difference(cc[r].keys())
        if toadd:
            cc[r].update(dict((k, cc["__default__"][k]) for k in toadd))
    return(cc)

with open(config["clusterfile"]) as f:
    clust_conf = default_clust(yaml.safe_load(f))


## temp directory
TMP = clust_conf["__default__"]["tmpdir"]
shell.prefix("export TMPDIR={TMP}; set -xe;")
