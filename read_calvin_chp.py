"""Minimal reader for Affymetrix Calvin (Command Console) genotype CHP files.

Extracts probeset ID, call code and confidence from Birdseed CHPs (as deposited
in GEO, e.g. GSE37384). Call codes: 6=AA, 7=BB, 8=AB, 11=NoCall.

Usage: read_calvin_chp.py OUT.tsv.gz PROBESET_LIST CHP.gz [CHP.gz ...]
Writes long format: sample, probeset, call, confidence (only probesets listed),
plus OUT.summary.tsv.gz with the affymetrix-chipsummary-* header values
(computed_gender, het_rate, chrX het rate, ...).
"""
import gzip
import struct
import sys

COL_FMT = {0: "b", 1: "B", 2: "h", 3: "H", 4: "i", 5: "I", 6: "f"}


class Buf:
    def __init__(self, data):
        self.d, self.p = data, 0

    def take(self, fmt):
        v = struct.unpack_from(">" + fmt, self.d, self.p)
        self.p += struct.calcsize(">" + fmt)
        return v[0]

    def raw(self):
        n = self.take("i")
        s = self.d[self.p:self.p + n]
        self.p += n
        return s

    def wstr(self):
        n = self.take("i")
        s = self.d[self.p:self.p + 2 * n].decode("utf-16-be")
        self.p += 2 * n
        return s


def param_value(raw, ptype):
    if ptype == "text/x-calvin-float":
        return struct.unpack(">f", raw[:4])[0]
    if ptype in ("text/x-calvin-integer-32", "text/x-calvin-unsigned-integer-32"):
        return struct.unpack(">i", raw[:4])[0]
    if ptype == "text/ascii":
        return raw.split(b"\0")[0].decode("ascii")
    if ptype == "text/plain":
        return raw.decode("utf-16-be").split("\0")[0]
    return None


def read_params(b, summary):
    for _ in range(b.take("i")):
        name = b.wstr(); raw = b.raw(); ptype = b.wstr()
        if name.startswith("affymetrix-chipsummary-"):
            summary[name[len("affymetrix-chipsummary-"):]] = param_value(raw, ptype)


def skip_generic_header(b, summary):
    b.raw()  # data type id
    b.raw()  # file id
    b.wstr()  # creation time
    b.wstr()  # locale
    read_params(b, summary)
    for _ in range(b.take("i")):
        skip_generic_header(b, summary)


def read_genotype_set(data):
    b = Buf(data)
    b.take("B"); b.take("B")
    n_groups = b.take("i")
    b.take("I")
    summary = {}
    skip_generic_header(b, summary)
    for _ in range(n_groups):
        next_group = b.take("I"); first_set = b.take("I")
        n_sets = b.take("i"); b.wstr()
        b.p = first_set
        for _ in range(n_sets):
            first_elem = b.take("I"); next_set = b.take("I")
            name = b.wstr()
            read_params(b, summary)
            cols = []
            for _ in range(b.take("I")):
                cols.append((b.wstr(), b.take("B"), b.take("i")))
            n_rows = b.take("I")
            if name == "Genotype":
                return parse_rows(data, first_elem, cols, n_rows), summary
            b.p = next_set
        b.p = next_group
    raise ValueError("no Genotype data set found")


def parse_rows(data, start, cols, n_rows):
    width = sum(c[2] for c in cols)
    out = {}
    for r in range(n_rows):
        off = start + r * width
        rec = {}
        for cname, ctype, csize in cols:
            if ctype == 7:
                n = struct.unpack_from(">i", data, off)[0]
                rec[cname] = data[off + 4:off + 4 + n].decode("ascii")
            elif ctype == 8:
                n = struct.unpack_from(">i", data, off)[0]
                rec[cname] = data[off + 4:off + 4 + 2 * n].decode("utf-16-be")
            else:
                rec[cname] = struct.unpack_from(">" + COL_FMT[ctype], data, off)[0]
            off += csize
        out[rec["ProbeSetName"]] = (rec["Call"], rec["Confidence"])
    return out


def main():
    out_path, list_path, *chps = sys.argv[1:]
    wanted = set(l.strip() for l in open(list_path) if l.strip())
    summaries = {}
    with gzip.open(out_path, "wt") as out:
        out.write("sample\tprobeset\tcall\tconfidence\n")
        for f in chps:
            opener = gzip.open if f.endswith(".gz") else open
            with opener(f, "rb") as fh:
                calls, summary = read_genotype_set(fh.read())
            sample = f.split("/")[-1].split(".")[0]
            summaries[sample] = summary
            for ps in wanted & calls.keys():
                c, conf = calls[ps]
                out.write(f"{sample}\t{ps}\t{c}\t{conf:.5f}\n")
    keys = sorted({k for v in summaries.values() for k in v})
    with gzip.open(out_path.replace(".tsv.gz", "") + ".summary.tsv.gz", "wt") as out:
        out.write("sample\t" + "\t".join(keys) + "\n")
        for smp, v in summaries.items():
            out.write(smp + "\t" + "\t".join(str(v.get(k, "")) for k in keys) + "\n")


if __name__ == "__main__":
    main()

