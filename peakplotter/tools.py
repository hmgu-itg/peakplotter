import shlex
import subprocess as sp
from pathlib import Path
from typing import List


class Plink:
    def __init__(self, memory = 30000):
        self.memory = memory
    
    @property
    def _cmd(self):
        return f'plink --memory {self.memory}'
    
    def __call__(self, command):
        return sp.run(shlex.split(f'{self._cmd} {command}'), stdout = sp.PIPE, stderr = sp.PIPE)
    
    def dry_call(self, command):
        return f'{self._cmd} {command}'
    
    def extract_genotypes(self, bfile, chrom, start, end, out):
        return sp.run(
            shlex.split(f'{self._cmd} --bfile {bfile} --chr {chrom} --from-bp {start} --to-bp {end} --out {out} --make-bed'), stdout = sp.PIPE, stderr = sp.PIPE
            )
    
    def merge(self, file: str, bfiles: List[str], chrom, start, end, out: str):
        with open(file, 'w') as f:
            for bfile in bfiles:
                f.write(f'{bfile}\n')
        process = sp.run(shlex.split(f'{self._cmd} --merge-list {file} --chr {chrom} --from-bp {start} --to-bp {end} --out {out} --make-bed'), stdout = sp.PIPE, stderr = sp.PIPE)
        # os.remove(file)
        return process
    
    def exclude(self, bfile, exclude, out):
        return sp.run(shlex.split(f'{self._cmd} --allow-no-sex --bfile {bfile} --exclude {exclude} --make-bed --out {out}'), stdout = sp.PIPE, stderr = sp.PIPE)
    
    def ld(self, bfile, ld_snp, ext_flank_kb, out):
        return sp.run(shlex.split(f'{self._cmd} --allow-no-sex --bfile {bfile} --r2 --ld-snp {ld_snp}  --ld-window-kb {ext_flank_kb} --ld-window 999999 --ld-window-r2 0 --out {out}'), stdout = sp.PIPE, stderr = sp.PIPE)


def plink_exclude_across_bfiles(plink: Plink, bfiles: List[Path], missnp_file: Path) -> List[Path]:
    outlist = list()
    for bfile in bfiles:
        excluded_bfile = Path(f'{bfile}.excluded')
        ps = plink.exclude(bfile, missnp_file, excluded_bfile)
        print(ps.stdout.decode())
        print(ps.stderr.decode())
        if ps.returncode==0:
            outlist.append(excluded_bfile)
        else:
            print(f'[WARNING] All variants in {bfile} were excluded')
    return outlist
