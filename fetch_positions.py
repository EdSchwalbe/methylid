"""Fetch genomic coordinates for marker rsIDs from the Ensembl REST API.

Needed because real VCFs often carry no rsIDs: the 1000 Genomes phase 3 release
files, for instance, have "." in the ID column for every variant. With a
position table, read_vcf_dosage() can join on chrom:pos instead.

Usage: fetch_positions.py FEATURE_TABLE OUT.tsv [grch37|grch38]
Writes: rsid, chrom, pos, alleles  (one row per rsID that resolves)
"""
import json
import sys
import time
import urllib.error
import urllib.request

HOSTS = {"grch37": "https://grch37.rest.ensembl.org",
         "grch38": "https://rest.ensembl.org"}
BATCH = 200
PRIMARY = {str(c) for c in range(1, 23)} | {"X", "Y", "MT"}


def post(host, ids, attempt=0):
    req = urllib.request.Request(
        host + "/variation/homo_sapiens",
        data=json.dumps({"ids": ids}).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.load(r)
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
        if attempt >= 4:
            print(f"  giving up on a batch of {len(ids)}: {e}", file=sys.stderr)
            return {}
        time.sleep(2 ** attempt)
        return post(host, ids, attempt + 1)


def main():
    feat_path, out_path = sys.argv[1], sys.argv[2]
    build = (sys.argv[3] if len(sys.argv) > 3 else "grch37").lower()
    host = HOSTS[build]

    with open(feat_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = header.index("rsid")
        rsids = sorted({l.split("\t")[col].strip() for l in fh if l.strip()})
    rsids = [r for r in rsids if r.startswith("rs")]
    print(f"{len(rsids)} rsIDs to resolve against {build}", file=sys.stderr)

    rows, missing = [], 0
    for i in range(0, len(rsids), BATCH):
        chunk = rsids[i:i + BATCH]
        data = post(host, chunk)
        for rs in chunk:
            v = data.get(rs)
            maps = [m for m in (v or {}).get("mappings", [])
                    if str(m.get("seq_region_name")) in PRIMARY]
            if not maps:
                missing += 1
                continue
            m = maps[0]
            rows.append((rs, str(m["seq_region_name"]), int(m["start"]),
                         m.get("allele_string", "")))
        print(f"  {min(i + BATCH, len(rsids))}/{len(rsids)}", end="\r", file=sys.stderr)
        time.sleep(0.1)

    with open(out_path, "w") as out:
        out.write("rsid\tchrom\tpos\talleles\n")
        for r in sorted(rows, key=lambda x: (x[1], x[2])):
            out.write("\t".join(map(str, r)) + "\n")
    print(f"\nwrote {out_path}: {len(rows)} resolved, {missing} unresolved",
          file=sys.stderr)


if __name__ == "__main__":
    main()
