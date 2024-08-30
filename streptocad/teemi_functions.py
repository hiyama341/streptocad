
from Bio import SeqIO  # For reading Genbank and FASTA files
import json  # For handling JSON data
import requests  # For making HTTP requests to the NEB API

def read_genbank_files(path):
    """Reads single Genbank files.
    Parameters
    ----------
    path: str
        path to the genbank file you want to read.

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord
    """

    ncbi_hits = []
    for seq_record in SeqIO.parse(path, format="gb"):
        ncbi_hits.append(seq_record)

    return ncbi_hits



def read_fasta_files(path):
    """Reads FASTA files.
    Parameters
    ----------
    path: str
        path to the fasta file you want to read.

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord
    """

    ncbi_hits = []
    for seq_record in SeqIO.parse(path, format="fasta"):
        ncbi_hits.append(seq_record)

    return ncbi_hits



def primer_tm_neb(primer, conc=0.5, prodcode="q5-0"):
    """Calculates a single primers melting temp from NEB.

    Parameters
    ----------
    primer1 : str
    conc : float
    prodcode : str
        find product codes on nebswebsite: https://tmapi.neb.com/docs/productcodes

    Returns
    -------
    tm : int
        primer melting temperature

    """

    url = "https://tmapi.neb.com/tm/batch"
    seqpairs = [[primer]]

    input = {"seqpairs": seqpairs, "conc": conc, "prodcode": prodcode}
    headers = {"content-type": "application/json"}
    res = requests.post(url, data=json.dumps(input), headers=headers)

    r = json.loads(res.content)

    if r["success"]:
        for row in r["data"]:
            return row["tm1"]
    else:
        print("request failed")
        print(r["error"][0])



def primer_ta_neb(primer1, primer2, conc=0.5, prodcode="q5-0"):
    """Calculates primer pair melting temp TA,  from NEB.

    Parameters
    ----------
    primer1 : str
        first primer to be used for finding the optimal ta
    primer2 : str
        second primer to be used for finding the optimal ta
    conc : float
    prodcode : str
        find product codes on nebswebsite: https://tmapi.neb.com/docs/productcodes

    Returns
    -------
    ta : int
        primer pair annealing temp

    """

    url = "https://tmapi.neb.com/tm/batch"
    seqpairs = [[primer1, primer2]]

    input = {"seqpairs": seqpairs, "conc": conc, "prodcode": prodcode}
    headers = {"content-type": "application/json"}
    res = requests.post(url, data=json.dumps(input), headers=headers)

    r = json.loads(res.content)

    if r["success"]:
        for row in r["data"]:
            return row["ta"]

    else:
        print("request failed")
        print(r["error"][0])


