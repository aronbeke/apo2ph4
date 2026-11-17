#!/usr/bin/env python3
"""
apo2ph4_ftmap.py -- drive the Apo2ph4 pharmacophore generation from FTMap/ATLAS output.

The stock Apo2ph4 workflow obtains posed fragments by docking a fragment library
with AutoDock Vina (apo2ph4_prepare_and_dock.sh).  FTMap/ATLAS already provides
posed probe molecules clustered into hotspots, so the docking step is replaced
while the rest of the workflow (interaction pharmacophore generation + grid
guided feature selection) is reused unchanged.

Per hotspot ("site") this script:

  1. collects the probe molecules of the selected consensus cluster(s),
  2. prepares the receptor and runs AutoGrid4 to obtain the grid maps that
     apo2ph4_generate_ph4.py requires (no docking is performed),
  3. writes feature_count.txt from the probe molecules (scoring normalisation),
  4. builds one pseudo-complex per probe and generates interaction
     pharmacophores (scripts/gen_ia_ph4s.py or scripts/ftmap_gen_ia_ph4s.py),
  5. runs apo2ph4_generate_ph4.py to produce the final apo pharmacophore,
  6. optionally converts the result to further CDPKit formats and writes a
     machine readable summary for downstream processing.

Nothing in the upstream repository is modified.

Usage examples
--------------
    # strongest hotspot only (default)
    python3 apo2ph4_ftmap.py atlas_output/1h24

    # the CS1/CS2/CS4 ensemble, merged into one site
    python3 apo2ph4_ftmap.py atlas_output/1h24 --clusters 1,2,4 --group-within 10

    # every consensus cluster, one model each, also written as .cdf
    python3 apo2ph4_ftmap.py atlas_output/1h24 --clusters all --formats pml,cdf
"""

import argparse
import json
import math
import os
import os.path as path
import re
import shutil
import subprocess
import sys

SCRIPT_DIR = path.dirname(path.abspath(__file__))
SCRIPTS_DIR = path.join(SCRIPT_DIR, 'scripts')
DATA_DIR = path.join(SCRIPT_DIR, 'data')

# grid map file names expected by apo2ph4_generate_ph4.py
REQUIRED_MAPS = ('receptor.C.map', 'receptor.e.map', 'receptor.HD.map',
                 'receptor.OA.map', 'receptor.A.map')

CONSENSUS_RE = re.compile(r'consensus\.(\d+)\.(\d+)\.sdf$')

LIG_RESNAME = 'LIG'
LIG_RESSEQ = 999
LIG_CHAIN = 'Z'


# --------------------------------------------------------------------------
# minimal SDF (V2000) handling
# --------------------------------------------------------------------------

class SDRecord:
    """One molecule of an SD-file, kept both parsed and verbatim."""

    def __init__(self, block):
        self.block = block
        lines = block.splitlines()

        if len(lines) < 4:
            raise ValueError('truncated SD record')

        self.title = lines[0].strip()

        counts = lines[3]

        if 'V3000' in counts:
            raise ValueError('V3000 SD-files are not supported')

        n_atoms = int(counts[0:3])
        n_bonds = int(counts[3:6])

        self.atoms = []   # (element, x, y, z)
        self.bonds = []   # (atom_idx_1, atom_idx_2, order)   1-based indices

        for line in lines[4:4 + n_atoms]:
            self.atoms.append((line[31:34].strip(),
                               float(line[0:10]), float(line[10:20]), float(line[20:30])))

        for line in lines[4 + n_atoms:4 + n_atoms + n_bonds]:
            self.bonds.append((int(line[0:3]), int(line[3:6]), int(line[6:9])))

        self.props = {}

        for i, line in enumerate(lines):
            if line.startswith('> <') and line.rstrip().endswith('>'):
                key = line.strip()[3:-1]
                value = lines[i + 1].strip() if i + 1 < len(lines) else ''
                self.props[key] = value

    def centroid(self):
        n = len(self.atoms)
        return [sum(a[i + 1] for a in self.atoms) / n for i in range(3)]


def read_sdf(file_path):
    records = []

    with open(file_path, 'r') as f:
        block = []

        for line in f:
            if line.startswith('$$$$'):
                if block:
                    records.append(SDRecord(''.join(block)))

                block = []
            else:
                block.append(line)

        if block and ''.join(block).strip():
            records.append(SDRecord(''.join(block)))

    return records


def write_sdf(records, file_path):
    with open(file_path, 'w') as f:
        for rec in records:
            block = rec.block

            if not block.endswith('\n'):
                block += '\n'

            f.write(block)
            f.write('$$$$\n')


# --------------------------------------------------------------------------
# PDB writing
# --------------------------------------------------------------------------

def _atom_name(element, serial):
    name = '%s%d' % (element.upper(), serial)
    name = name[:4]

    # one-letter elements are indented by one column in the PDB atom name field
    return (' %-3s' % name) if len(element) == 1 else ('%-4s' % name)


def pdb_hetatm_lines(record, resname=LIG_RESNAME, resseq=LIG_RESSEQ,
                     chain=LIG_CHAIN, first_serial=1):
    """Render an SD record as HETATM + CONECT lines (bond orders preserved)."""
    lines = []

    for i, (element, x, y, z) in enumerate(record.atoms):
        serial = first_serial + i
        lines.append('HETATM%5d %s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s'
                     % (serial, _atom_name(element, i + 1), resname, chain, resseq,
                        x, y, z, 1.0, 0.0, element.upper()))

    # CONECT records; a bond of order n is encoded by listing the partner n times,
    # which CDPKit's PDB reader turns back into the correct bond order
    # (Biomol.ControlParameter.PDB_DEDUCE_BOND_ORDERS_FROM_CONECT_RECORDS).
    partners = {}

    for a1, a2, order in record.bonds:
        order = max(1, min(3, order if order in (1, 2, 3) else 1))
        partners.setdefault(a1, []).extend([a2] * order)
        partners.setdefault(a2, []).extend([a1] * order)

    for atom_idx in sorted(partners):
        bonded = partners[atom_idx]

        for chunk_start in range(0, len(bonded), 4):
            chunk = bonded[chunk_start:chunk_start + 4]
            line = 'CONECT%5d' % (first_serial + atom_idx - 1)

            for p in chunk:
                line += '%5d' % (first_serial + p - 1)

            lines.append(line)

    return lines


def read_receptor_atom_lines(receptor_pdb):
    """ATOM records of the receptor, HETATM/CONECT/anisotropic records dropped."""
    keep = []

    with open(receptor_pdb, 'r') as f:
        for line in f:
            if line.startswith('ATOM  '):
                keep.append(line.rstrip('\n'))
            elif line.startswith('TER'):
                keep.append('TER')

    if not keep:
        sys.exit('Error: no ATOM records found in receptor file \'%s\'' % receptor_pdb)

    return keep


def renumber_pdb_lines(lines, first_serial):
    out = []
    serial = first_serial

    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            out.append('%s%5d%s' % (line[:6], serial, line[11:]))
            serial += 1
        else:
            out.append(line)

    return out, serial


def write_complex_pdb(record, receptor_lines, out_path):
    """Pseudo co-crystal complex: probe atoms FIRST, then the receptor.

    scripts/gen_ia_ph4s.py identifies the ligand as the first atom of the file
    plus every atom sharing its residue code and sequence number, so the probe
    must come first and must use a residue code/number absent from the protein.
    """
    lig_lines = pdb_hetatm_lines(record)
    n_lig_atoms = len(record.atoms)
    rec_lines, _ = renumber_pdb_lines(receptor_lines, n_lig_atoms + 1)

    with open(out_path, 'w') as f:
        f.write('REMARK   FTMap probe \'%s\' + receptor (apo2ph4_ftmap.py)\n' % record.title)

        for line in lig_lines:
            f.write(line + '\n')

        for line in rec_lines:
            f.write(line + '\n')

        f.write('END\n')


def write_dummy_ligand_pdb(center, out_path):
    """Benzene from data/benzene.sdf translated to `center`; defines the grid box."""
    benzene_sdf = path.join(DATA_DIR, 'benzene.sdf')

    if not path.isfile(benzene_sdf):
        sys.exit('Error: \'%s\' not found' % benzene_sdf)

    rec = read_sdf(benzene_sdf)[0]
    cen = rec.centroid()
    shift = [center[i] - cen[i] for i in range(3)]
    rec.atoms = [(e, x + shift[0], y + shift[1], z + shift[2]) for (e, x, y, z) in rec.atoms]

    with open(out_path, 'w') as f:
        f.write('REMARK   grid box definition dummy ligand (apo2ph4_ftmap.py)\n')

        for line in pdb_hetatm_lines(rec, resname='BNZ', resseq=1, chain='Z'):
            f.write(line + '\n')

        f.write('END\n')


# --------------------------------------------------------------------------
# site selection
# --------------------------------------------------------------------------

class Cluster:
    def __init__(self, file_path, rank, size, records):
        self.file_path = file_path
        self.rank = rank
        self.size = size
        self.records = records

    @property
    def centroid(self):
        pts = [a[1:] for rec in self.records for a in rec.atoms]
        n = len(pts)

        return [sum(p[i] for p in pts) / n for i in range(3)]

    @property
    def best_energy(self):
        energies = []

        for rec in self.records:
            try:
                energies.append(float(rec.props.get('ATLAS_ENERGY')))
            except (TypeError, ValueError):
                pass

        return min(energies) if energies else None


def discover_clusters(atlas_dir):
    clusters = []

    for name in sorted(os.listdir(atlas_dir)):
        m = CONSENSUS_RE.search(name)

        if not m:
            continue

        file_path = path.join(atlas_dir, name)

        try:
            records = read_sdf(file_path)
        except Exception as e:
            sys.exit('Error: parsing \'%s\' failed: %s' % (file_path, e))

        if not records:
            continue

        clusters.append(Cluster(file_path, int(m.group(1)), int(m.group(2)), records))

    if not clusters:
        sys.exit('Error: no *_consensus.<rank>.<size>.sdf files found in \'%s\'' % atlas_dir)

    clusters.sort(key=lambda c: c.rank)

    return clusters


def select_clusters(clusters, spec):
    by_rank = {c.rank: c for c in clusters}

    if spec == 'all':
        return list(clusters)

    if spec == 'top':
        return [clusters[0]]

    m = re.fullmatch(r'top:(\d+)', spec)

    if m:
        return clusters[:int(m.group(1))]

    selected = []

    for token in spec.split(','):
        token = token.strip()

        if not token:
            continue

        try:
            rank = int(token)
        except ValueError:
            sys.exit('Error: invalid --clusters value \'%s\'' % token)

        if rank not in by_rank:
            sys.exit('Error: consensus cluster %d not present (available: %s)'
                     % (rank, ', '.join(str(c.rank) for c in clusters)))

        selected.append(by_rank[rank])

    if not selected:
        sys.exit('Error: --clusters selected nothing')

    return selected


def group_clusters(selected, cutoff):
    """Single-linkage grouping of clusters whose centroids are within `cutoff`."""
    if cutoff <= 0 or len(selected) == 1:
        return [[c] for c in selected]

    groups = [[c] for c in selected]

    merged = True

    while merged:
        merged = False

        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                close = False

                for a in groups[i]:
                    for b in groups[j]:
                        d = math.dist(a.centroid, b.centroid)

                        if d <= cutoff:
                            close = True
                            break

                    if close:
                        break

                if close:
                    groups[i] = groups[i] + groups[j]
                    del groups[j]
                    merged = True
                    break

            if merged:
                break

    for g in groups:
        g.sort(key=lambda c: c.rank)

    groups.sort(key=lambda g: g[0].rank)

    return groups


# --------------------------------------------------------------------------
# external command handling
# --------------------------------------------------------------------------

class Runner:
    def __init__(self, dry_run=False, verbose=False):
        self.dry_run = dry_run
        self.verbose = verbose
        self.env = self._build_env()

    @staticmethod
    def _build_env():
        env = os.environ.copy()
        conda_prefix = env.get('CONDA_PREFIX')

        if conda_prefix:
            # MGLTools utilities (prepare_receptor4.py etc.), cf. apo2ph4_prepare_and_dock.sh
            utils = path.join(conda_prefix, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24')

            if path.isdir(utils):
                env['PATH'] = env.get('PATH', '') + os.pathsep + utils

            try:
                version = subprocess.check_output(['obabel', '-V'], env=env,
                                                  stderr=subprocess.DEVNULL).decode().split()[2]
                env['BABEL_LIBDIR'] = path.join(conda_prefix, 'lib', 'openbabel', version)
                env['BABEL_DATADIR'] = path.join(conda_prefix, 'share', 'openbabel', version)
            except Exception:
                pass

        return env

    def run(self, cmd, cwd=None, capture=False, log_path=None):
        printable = ' '.join(cmd)

        if self.dry_run:
            print('  [dry-run] (cd %s && %s)' % (cwd or os.getcwd(), printable))
            return ''

        if self.verbose:
            print('  $ %s' % printable)

        try:
            proc = subprocess.run(cmd, cwd=cwd, env=self.env, text=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT if capture else None)
        except FileNotFoundError:
            sys.exit('Error: executable \'%s\' not found -- is the apo2ph4 conda '
                     'environment active?' % cmd[0])

        output = proc.stdout or ''

        if log_path and output:
            with open(log_path, 'w') as f:
                f.write(output)

        if proc.returncode != 0:
            if output and not self.verbose:
                sys.stderr.write(output)

            sys.exit('Error: command failed (exit %d): %s' % (proc.returncode, printable))

        return output


# --------------------------------------------------------------------------
# per-site workflow
# --------------------------------------------------------------------------

def prepare_grids(site_dir, receptor_pdb, center, args, runner):
    """prepare_receptor4 / prepare_ligand4 / prepare_gpf4 / autogrid4 -- no docking."""
    existing = [path.join(site_dir, m) for m in REQUIRED_MAPS]

    if args.skip_grids and all(path.isfile(p) for p in existing):
        print('  - grid maps already present, skipped (--skip-grids)')
        return

    shutil.copyfile(receptor_pdb, path.join(site_dir, 'receptor.pdb'))
    write_dummy_ligand_pdb(center, path.join(site_dir, 'gridbox_ligand.pdb'))

    runner.run(['prepare_receptor4.py', '-r', 'receptor.pdb', '-o', 'receptor.pdbqt'],
               cwd=site_dir, capture=True)
    runner.run(['prepare_ligand4.py', '-l', 'gridbox_ligand.pdb', '-o', 'gridbox_ligand.pdbqt'],
               cwd=site_dir, capture=True)

    gpf_cmd = ['prepare_gpf4.py',
               '-l', 'gridbox_ligand.pdbqt',
               '-r', 'receptor.pdbqt',
               '-o', 'receptor.gpf',
               '-p', 'npts=%d,%d,%d' % (args.npts, args.npts, args.npts),
               '-p', 'ligand_types=C, OA, HD, A',
               '-p', 'spacing=%.3f' % args.spacing]

    if args.grid_center:
        gpf_cmd += ['-p', 'gridcenter=%.3f,%.3f,%.3f' % tuple(args.grid_center)]
    else:
        gpf_cmd += ['-y']

    runner.run(gpf_cmd, cwd=site_dir, capture=True)
    runner.run(['autogrid4', '-p', 'receptor.gpf', '-l', 'receptor.glg'],
               cwd=site_dir, capture=True)

    if runner.dry_run:
        return

    missing = [m for m in REQUIRED_MAPS if not path.isfile(path.join(site_dir, m))]

    if missing:
        sys.exit('Error: AutoGrid did not produce %s in \'%s\' -- see receptor.glg'
                 % (', '.join(missing), site_dir))


def build_complexes(site_dir, records, receptor_pdb, runner):
    complex_dir = path.join(site_dir, 'complexes')
    os.makedirs(complex_dir, exist_ok=True)

    if runner.dry_run:
        print('  [dry-run] would write %d pseudo-complexes to %s' % (len(records), complex_dir))
        return path.join(site_dir, 'complex_list.list')

    receptor_lines = read_receptor_atom_lines(receptor_pdb)
    paths = []

    for i, rec in enumerate(records):
        out_path = path.join(complex_dir, 'probe_complex%03d.pdb' % i)
        write_complex_pdb(rec, receptor_lines, out_path)
        paths.append(path.abspath(out_path))

    list_path = path.join(site_dir, 'complex_list.list')

    with open(list_path, 'w') as f:
        for p in paths:
            f.write(p + '\n')

    return list_path


def convert_pharmacophore(pml_path, formats):
    """Write the final model in further CDPKit formats (e.g. .cdf) next to the .pml."""
    extra = [f.strip().lstrip('.').lower() for f in formats.split(',') if f.strip()]
    extra = [f for f in extra if f != 'pml']

    if not extra:
        return []

    try:
        import CDPL.Pharm as Pharm
        import CDPL.Base as Base
    except ImportError:
        print('  ! CDPKit not importable, skipping format conversion')
        return []

    written = []
    stem = path.splitext(pml_path)[0]

    for fmt in extra:
        out_path = '%s.%s' % (stem, fmt)

        try:
            reader = Pharm.PMLPharmacophoreReader(Base.FileIOStream(pml_path, 'r'))
            pharm = Pharm.BasicPharmacophore()

            if not reader.read(pharm):
                print('  ! reading \'%s\' failed, skipping .%s' % (pml_path, fmt))
                continue

            writer = Pharm.FeatureContainerWriter(out_path)

            if not writer.write(pharm):
                print('  ! writing \'%s\' failed' % out_path)
                continue

            writer.close()
            written.append(out_path)
        except Exception as e:
            print('  ! conversion to .%s failed: %s' % (fmt, e))

    return written


def parse_ph4_report(output):
    """Extract the per-feature table printed by apo2ph4_generate_ph4.py."""
    features = []
    in_table = False

    for line in output.splitlines():
        if line.startswith('Type'):
            in_table = True
            continue

        if not in_table:
            continue

        fields = line.split()

        if len(fields) < 7:
            continue

        try:
            features.append({'type': fields[0],
                             'x': float(fields[1]), 'y': float(fields[2]), 'z': float(fields[3]),
                             'score': float(fields[4]),
                             'buriedness_psp': int(fields[5]),
                             'grid_energy': float(fields[6])})
        except ValueError:
            continue

    return features


def process_site(site_index, group, out_dir, receptor_pdb, args, runner):
    ranks = [c.rank for c in group]
    site_name = 'site_%02d_cs%s' % (site_index, '-'.join('%03d' % r for r in ranks))
    site_dir = path.join(out_dir, site_name)
    os.makedirs(site_dir, exist_ok=True)

    records = [rec for c in group for rec in c.records]
    pts = [a[1:] for rec in records for a in rec.atoms]
    center = [sum(p[i] for p in pts) / len(pts) for i in range(3)]
    extent = [max(p[i] for p in pts) - min(p[i] for p in pts) for i in range(3)]

    print('\n=== %s ===' % site_name)
    print('  consensus clusters : %s' % ', '.join('CS%d (%d probes)' % (c.rank, c.size) for c in group))
    print('  probe molecules    : %d' % len(records))
    print('  centroid           : %.2f %.2f %.2f' % tuple(center))
    print('  probe extent       : %.1f x %.1f x %.1f A' % tuple(extent))

    box = args.npts * args.spacing

    if max(extent) > box:
        print('  ! probe extent exceeds the %.1f A grid box -- increase --npts or --spacing' % box)

    probes_sdf = path.join(site_dir, 'probes.sdf')

    if not runner.dry_run:
        write_sdf(records, probes_sdf)

    print('  - preparing receptor and AutoGrid maps (no docking)')
    prepare_grids(site_dir, receptor_pdb, center, args, runner)

    print('  - counting probe pharmacophore features (feature_count.txt)')
    runner.run([sys.executable, path.join(SCRIPTS_DIR, 'count_features.py'), 'probes.sdf'],
               cwd=site_dir, capture=True)

    ph4_dir = path.join(site_dir, 'pharmacophores')
    os.makedirs(ph4_dir, exist_ok=True)
    frag_ph4s = path.join('pharmacophores', 'pharmacophores.pml')

    if args.ph4_mode == 'pdb':
        print('  - building %d pseudo-complexes' % len(records))
        build_complexes(site_dir, records, receptor_pdb, runner)
        print('  - generating interaction pharmacophores (gen_ia_ph4s.py)')
        runner.run([sys.executable, path.join(SCRIPTS_DIR, 'gen_ia_ph4s.py'),
                    'complex_list.list', frag_ph4s],
                   cwd=site_dir, capture=True)
    else:
        print('  - generating probe pharmacophores (%s mode)' % args.ph4_mode)
        cmd = [sys.executable, path.join(SCRIPTS_DIR, 'ftmap_gen_ia_ph4s.py'),
               '-r', 'receptor.pdb', '-p', 'probes.sdf', '-o', frag_ph4s,
               '-m', args.ph4_mode]

        if args.hb_max_length is not None:
            cmd += ['--hb-max-length', str(args.hb_max_length)]

        out = runner.run(cmd, cwd=site_dir, capture=True)

        for line in out.splitlines():
            if line.startswith(('  ', '- Done', '- Relaxed')):
                print('   ' + line.strip())

    out_pml = path.join('pharmacophores', '%s.pml' % site_name)
    model_name = args.name if args.name else '%s %s' % (args.label or 'apo2ph4-ftmap', site_name)

    gen_cmd = [sys.executable, path.join(SCRIPT_DIR, 'apo2ph4_generate_ph4.py'),
               '-i', frag_ph4s, '-o', out_pml, '-g', '.',
               '--name', model_name,
               '-n', str(args.num_features),
               '-t', str(args.distance_threshold),
               '-H', str(args.max_hydrophobic),
               '-D', str(args.max_HBD),
               '-P', str(args.max_PI),
               '-N', str(args.max_NI),
               '--H_energy', str(args.H_energy),
               '--PI_energy', str(args.PI_energy),
               '--NI_energy', str(args.NI_energy),
               '-v']

    if args.ph4_args:
        gen_cmd += args.ph4_args.split()

    print('  - generating the apo pharmacophore')
    output = runner.run(gen_cmd, cwd=site_dir, capture=True,
                        log_path=path.join(site_dir, 'ph4_report.txt'))

    if runner.dry_run:
        return None

    print(output.rstrip())

    abs_pml = path.join(site_dir, out_pml)
    converted = convert_pharmacophore(abs_pml, args.formats)

    summary = {
        'site': site_name,
        'consensus_clusters': [{'rank': c.rank,
                                'size': c.size,
                                'centroid': [round(v, 3) for v in c.centroid],
                                'best_atlas_energy': c.best_energy,
                                'source_file': path.abspath(c.file_path)} for c in group],
        'n_probes': len(records),
        'probe_titles': [rec.title for rec in records],
        'site_centroid': [round(v, 3) for v in center],
        'probe_extent': [round(v, 3) for v in extent],
        'grid_box_edge_A': round(box, 2),
        'ph4_mode': args.ph4_mode,
        'pharmacophore_pml': path.abspath(abs_pml),
        'pharmacophore_other_formats': [path.abspath(p) for p in converted],
        'fragment_pharmacophores_pml': path.abspath(path.join(site_dir, frag_ph4s)),
        'features': parse_ph4_report(output),
    }

    with open(path.join(site_dir, 'site_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    if not args.keep_complexes:
        shutil.rmtree(path.join(site_dir, 'complexes'), ignore_errors=True)

    return summary


# --------------------------------------------------------------------------

def find_receptor(atlas_dir, explicit):
    if explicit:
        if not path.isfile(explicit):
            sys.exit('Error: receptor file \'%s\' not found' % explicit)

        return explicit

    candidates = [n for n in sorted(os.listdir(atlas_dir)) if n.endswith('_atlas_protein.pdb')]

    if candidates:
        return path.join(atlas_dir, candidates[0])

    sys.exit('Error: no *_atlas_protein.pdb in \'%s\'; pass --receptor explicitly' % atlas_dir)


def parse_args():
    p = argparse.ArgumentParser(
        description='Generate Apo2ph4 pharmacophores from FTMap/ATLAS hotspots.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('atlas_dir', help='directory holding the FTMap/ATLAS output')

    p.add_argument('-o', '--out-dir', default='ftmap_ph4',
                   help='output directory (one sub-directory per site)')
    p.add_argument('-r', '--receptor', default=None,
                   help='receptor PDB (default: the *_atlas_protein.pdb of ATLAS_DIR)')

    sel = p.add_argument_group('hotspot selection')
    sel.add_argument('-c', '--clusters', default='top',
                     help="'top', 'top:N', 'all' or a comma separated list of "
                          'consensus cluster ranks (e.g. 0,1,2)')
    sel.add_argument('--group-within', type=float, default=8.0,
                     help='merge selected clusters whose centroids are within this '
                          'distance (A) into a single site; 0 disables merging')

    grid = p.add_argument_group('grid generation')
    grid.add_argument('--npts', type=int, default=100, help='AutoGrid points per dimension')
    grid.add_argument('--spacing', type=float, default=0.300, help='AutoGrid spacing in A')
    grid.add_argument('--grid-center', type=float, nargs=3, default=None,
                      metavar=('X', 'Y', 'Z'),
                      help='explicit grid centre (default: the probe centroid of the site)')
    grid.add_argument('--skip-grids', action='store_true',
                      help='reuse existing grid maps if present')

    ph4 = p.add_argument_group('pharmacophore generation')
    ph4.add_argument('--ph4-mode', choices=('ligand', 'interaction', 'pdb'), default='ligand',
                     help="'ligand': all probe features + environment exclusion volumes; "
                          "'interaction': only probe features with a valid interaction "
                          "partner; 'pdb': the same, via pseudo-complex PDBs and the "
                          'upstream gen_ia_ph4s.py')
    ph4.add_argument('--hb-max-length', type=float, default=None, metavar='A',
                     help='relax the maximum H...acceptor distance (CDPKit default 2.8 A); '
                          'interaction mode only')
    ph4.add_argument('--formats', default='pml',
                     help='comma separated output formats, e.g. pml,cdf')
    ph4.add_argument('--keep-complexes', action='store_true',
                     help='keep the pseudo-complex PDBs (pdb mode only)')
    ph4.add_argument('--label', default=None, help='prefix for the pharmacophore model name')
    ph4.add_argument('--name', default=None,
                     help='explicit pharmacophore model name (overrides --label)')

    fwd = p.add_argument_group('apo2ph4_generate_ph4.py settings')
    fwd.add_argument('-n', '--num-features', type=int, default=0,
                     help='number of features, 0 = automatic')
    fwd.add_argument('-t', '--distance-threshold', type=int, default=3)
    fwd.add_argument('-H', '--max-hydrophobic', type=int, default=4)
    fwd.add_argument('-D', '--max-HBD', type=int, default=2)
    fwd.add_argument('-P', '--max-PI', type=int, default=1)
    fwd.add_argument('-N', '--max-NI', type=int, default=1)
    fwd.add_argument('--H_energy', type=float, default=-0.45)
    fwd.add_argument('--PI_energy', type=float, default=-0.7)
    fwd.add_argument('--NI_energy', type=float, default=0.9)
    fwd.add_argument('--ph4-args', default=None,
                     help='extra arguments forwarded verbatim to apo2ph4_generate_ph4.py')

    p.add_argument('--dry-run', action='store_true', help='print the steps without running them')
    p.add_argument('-v', '--verbose', action='store_true', help='echo every external command')

    return p.parse_args()


def main():
    args = parse_args()

    if not path.isdir(args.atlas_dir):
        sys.exit('Error: \'%s\' is not a directory' % args.atlas_dir)

    receptor_pdb = find_receptor(args.atlas_dir, args.receptor)
    clusters = discover_clusters(args.atlas_dir)

    print('FTMap/ATLAS directory : %s' % path.abspath(args.atlas_dir))
    print('Receptor              : %s' % path.abspath(receptor_pdb))
    print('Consensus clusters    : %s'
          % ', '.join('CS%d(%d)' % (c.rank, c.size) for c in clusters))

    selected = select_clusters(clusters, args.clusters)
    groups = group_clusters(selected, args.group_within)

    print('Selected              : %s' % ', '.join('CS%d' % c.rank for c in selected))
    print('Sites to process      : %d' % len(groups))

    os.makedirs(args.out_dir, exist_ok=True)
    summaries = []

    for i, group in enumerate(groups):
        summary = process_site(i, group, args.out_dir, path.abspath(receptor_pdb), args,
                               Runner(args.dry_run, args.verbose))

        if summary:
            summaries.append(summary)

    if summaries:
        index_path = path.join(args.out_dir, 'sites.json')

        with open(index_path, 'w') as f:
            json.dump(summaries, f, indent=2)

        print('\n- Done. %d pharmacophore model(s) written.' % len(summaries))
        print('  Index: %s' % path.abspath(index_path))

        for s in summaries:
            print('   %-28s %2d features  %s'
                  % (s['site'], len(s['features']), s['pharmacophore_pml']))


if __name__ == '__main__':
    main()
