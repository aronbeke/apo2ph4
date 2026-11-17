#!/usr/bin/env python3
"""
ftmap_gen_ia_ph4s.py -- interaction pharmacophores for FTMap/ATLAS probes.

Functional equivalent of scripts/gen_ia_ph4s.py, but instead of reading
pre-assembled pseudo co-crystal complex PDB files it reads the receptor once
and combines it in memory with every probe molecule of an SD-file.

Why: the probes keep the bond orders, formal charges and stereo information of
the SD-file instead of going through a PDB round trip, and the receptor is
parsed once rather than once per probe.

Two modes are available:

  interaction  only probe features that have a geometrically valid interaction
               partner in the receptor are kept -- what gen_ia_ph4s.py does.
  ligand       (default) every probe feature is kept, and the exclusion volumes
               are taken from the receptor environment.

FTMap probe poses come from a mapping/clustering procedure, not from a docking
run optimised against interaction geometry, so under the strict interaction
criteria most probe features are discarded (for the bundled 1H24 CS0 hotspot:
12 of 47). The probe positions themselves are the hotspot information, and
apo2ph4_generate_ph4.py re-filters everything by AutoGrid energies anyway,
which is why 'ligand' is the default here.

The emitted pharmacophores are written to a single multi-model file whose format
follows from the file extension (.pml for LigandScout, .cdf for CDPKit).

    python3 ftmap_gen_ia_ph4s.py -r receptor.pdb -p probes.sdf -o pharmacophores.pml
"""

import argparse
import math
import sys

import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Pharm as Pharm
import CDPL.MolProp as MolProp


LIG_RESIDUE_CODE = 'LIG'
LIG_RESIDUE_SEQ_NO = 999
LIG_CHAIN_ID = 'Z'

# maximum distance (A) at which a probe H-bond feature is paired with a
# complementary receptor feature in 'ligand' mode
DEF_HB_PAIR_RADIUS = 4.5


def pair_hbond_feature(ftr, ftr_type, env_donors, env_acceptors, radius):
    """Anchor a probe H-bond feature on the nearest complementary receptor feature.

    apo2ph4_generate_ph4.py clusters H-bond features by their *origin*
    (position +/- length * orientation) with DBSCAN at eps=0.01 A, i.e. it
    expects features of different fragments that address the same receptor atom
    to share an origin exactly.  Interaction pharmacophores satisfy this because
    their orientation points at the interaction partner.  Probe features
    generated without a partner point along the probe's own bonds instead, so
    their origins never coincide and the whole feature type is discarded
    downstream.

    This function restores that invariant using a pure distance criterion rather
    than the strict H-bond geometry that FTMap probe poses rarely satisfy: the
    feature is pointed at the nearest complementary receptor feature within
    `radius`, and its length is set to that distance.  The origin computed by
    apo2ph4_generate_ph4.py then lands exactly on the receptor feature.

    Returns True if a partner was found.
    """
    partners = env_acceptors if ftr_type == Pharm.FeatureType.H_BOND_DONOR else env_donors

    if not partners:
        return False

    pos = Chem.get3DCoordinates(ftr)
    best = None
    best_dist = radius

    for p in partners:
        dx, dy, dz = p[0] - pos[0], p[1] - pos[1], p[2] - pos[2]
        dist = math.sqrt(dx * dx + dy * dy + dz * dz)

        if dist < best_dist and dist > 0.0:
            best, best_dist = (dx, dy, dz), dist

    if best is None:
        return False

    orientation = (best[0] / best_dist, best[1] / best_dist, best[2] / best_dist)

    if ftr_type == Pharm.FeatureType.H_BOND_ACCEPTOR:
        # generate_ph4 computes the acceptor origin as pos - length * orientation
        orientation = (-orientation[0], -orientation[1], -orientation[2])

    Pharm.setOrientation(ftr, orientation)
    Pharm.setLength(ftr, best_dist)
    Pharm.setGeometry(ftr, Pharm.FeatureGeometry.VECTOR)

    return True


def env_feature_positions(env_ph4, ftr_type):
    return [Chem.get3DCoordinates(f) for f in env_ph4 if Pharm.getType(f) == ftr_type]


def read_receptor(file_path: str) -> Chem.BasicMolecule:
    mol = Chem.BasicMolecule()
    reader = Chem.MoleculeReader(file_path)

    try:
        if not reader.read(mol):
            sys.exit(f'Error: reading receptor file \'{file_path}\' failed')

    except Exception as e:
        sys.exit(f'Error: reading receptor file \'{file_path}\' failed: {str(e)}')

    return mol


def tag_as_ligand(mol: Chem.Molecule) -> None:
    """Give every probe atom a residue identity that the protein cannot collide with."""
    for i, atom in enumerate(mol.atoms):
        Biomol.setResidueCode(atom, LIG_RESIDUE_CODE)
        Biomol.setResidueSequenceNumber(atom, LIG_RESIDUE_SEQ_NO)
        Biomol.setChainID(atom, LIG_CHAIN_ID)
        Biomol.setHeteroAtomFlag(atom, True)
        Biomol.setResidueAtomName(atom, f'{Chem.getSymbolForType(atom).upper()}{i + 1}'[:4])


def prepare_complex(compl_mol: Chem.Molecule) -> None:
    """Same preprocessing as scripts/gen_ia_ph4s.py (existing information is kept)."""
    try:
        Chem.perceiveSSSR(compl_mol, False)
        Chem.setRingFlags(compl_mol, False)
        Chem.calcImplicitHydrogenCounts(compl_mol, False)
        Chem.perceiveHybridizationStates(compl_mol, False)
        Chem.setAromaticityFlags(compl_mol, False)

        if Chem.makeHydrogenComplete(compl_mol):                    # make implicit hydrogens explicit
            Chem.calcHydrogen3DCoordinates(compl_mol)
            Biomol.setHydrogenResidueSequenceInfo(compl_mol, False)

        MolProp.calcAtomHydrophobicities(compl_mol, False)          # needed for hydrophobic features

    except Exception as e:
        sys.exit(f'Error: processing of the probe/receptor complex failed: {str(e)}')


def split_complex(compl_mol: Chem.Molecule) -> tuple:
    """Split into (receptor, ligand) by the residue identity set in tag_as_ligand()."""
    rec_mg = Chem.Fragment()
    lig_mg = Chem.Fragment()

    for atom in compl_mol.atoms:
        if (Biomol.getResidueCode(atom) == LIG_RESIDUE_CODE and
                Biomol.getResidueSequenceNumber(atom) == LIG_RESIDUE_SEQ_NO):
            lig_mg.addAtom(atom)
        else:
            rec_mg.addAtom(atom)

    for bond in compl_mol.bonds:
        if lig_mg.containsAtom(bond.getBegin()) and lig_mg.containsAtom(bond.getEnd()):
            lig_mg.addBond(bond)
        else:
            rec_mg.addBond(bond)

    if lig_mg.numAtoms == 0:
        sys.exit('Error: the probe molecule was lost while splitting the complex')

    Chem.extractSSSRSubset(compl_mol, lig_mg, True)
    Chem.extractSSSRSubset(compl_mol, rec_mg, True)

    return (rec_mg, lig_mg)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Generate interaction pharmacophores for FTMap/ATLAS probe molecules.')

    parser.add_argument('-r', '--receptor', required=True,
                        help='receptor structure (PDB)')
    parser.add_argument('-p', '--probes', required=True,
                        help='posed probe molecules (SD-file)')
    parser.add_argument('-o', '--output', required=True,
                        help='output pharmacophore file (.pml, .cdf, ...)')
    parser.add_argument('-m', '--mode', choices=('interaction', 'ligand'), default='ligand',
                        help="'interaction': keep only probe features that have a valid "
                             'interaction partner in the receptor (as in gen_ia_ph4s.py); '
                             "'ligand': keep all probe features and take the exclusion "
                             'volumes from the receptor environment (default -- FTMap probe '
                             'poses are not geometry-optimised against the receptor and lose '
                             'most features under the interaction criteria)')
    parser.add_argument('--hb-max-length', type=float, default=None,
                        metavar='A',
                        help='relax the maximum H...acceptor distance of the H-bond criterion '
                             f'(CDPKit default: {Pharm.HBondingInteractionConstraint.DEF_MAX_HB_LENGTH} A); '
                             'interaction mode only')
    parser.add_argument('--hb-pair-radius', type=float, default=DEF_HB_PAIR_RADIUS,
                        metavar='A',
                        help='maximum distance at which a probe H-bond feature is anchored '
                             'on a complementary receptor feature; ligand mode only')
    parser.add_argument('--no-exclusion-volumes', action='store_true',
                        help='do not generate exclusion volume spheres')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='suppress per-probe progress output')

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    receptor = read_receptor(args.receptor)

    if not args.quiet:
        print(f'- Receptor \'{args.receptor}\': {receptor.numAtoms} atoms')

    probe_reader = Chem.MoleculeReader(args.probes)
    probe = Chem.BasicMolecule()

    ph4_writer = Pharm.FeatureContainerWriter(args.output)
    ia_ph4 = Pharm.BasicPharmacophore()
    out_ph4 = Pharm.BasicPharmacophore()

    ph4_gen = Pharm.InteractionPharmacophoreGenerator()

    ph4_gen.addExclusionVolumes(not args.no_exclusion_volumes)

    # disable unsupported feature types (as in scripts/gen_ia_ph4s.py)
    ph4_gen.envPharmacophoreGenerator.enableFeature(Pharm.FeatureType.HALOGEN_BOND_ACCEPTOR, False)
    ph4_gen.corePharmacophoreGenerator.enableFeature(Pharm.FeatureType.HALOGEN_BOND_DONOR, False)

    if args.hb_max_length is not None:
        if args.mode != 'interaction':
            sys.exit('Error: --hb-max-length is only meaningful in interaction mode')

        hbd = Pharm.FeatureType.H_BOND_DONOR
        hba = Pharm.FeatureType.H_BOND_ACCEPTOR
        min_len = Pharm.HBondingInteractionConstraint.DEF_MIN_HB_LENGTH
        analyzer = ph4_gen.getInteractionAnalyzer()

        # first argument states whether the FIRST feature of the pair is the donor
        analyzer.setConstraintFunction(
            hbd, hba, Pharm.HBondingInteractionConstraint(True, min_len, args.hb_max_length))
        analyzer.setConstraintFunction(
            hba, hbd, Pharm.HBondingInteractionConstraint(False, min_len, args.hb_max_length))

        print(f'- Relaxed H-bond H...acceptor distance to {min_len} - {args.hb_max_length} A')

    n_core_ftrs = 0
    n_out_ftrs = 0
    n_unpaired = [0]
    n_kept_ftrs = [0]
    count = 0

    while True:
        try:
            if not probe_reader.read(probe):
                break

        except Exception as e:
            sys.exit(f'Error: reading probe molecule {count + 1} failed: {str(e)}')

        name = Chem.getName(probe) or f'probe_{count + 1}'

        if not args.quiet:
            print(f'- Processing probe \'{name}\'')

        try:
            compl_mol = Chem.BasicMolecule()

            compl_mol.assign(probe)       # the probe has to come first
            tag_as_ligand(compl_mol)

            compl_mol += receptor         # append the receptor

            prepare_complex(compl_mol)

            rec_mg, lig_mg = split_complex(compl_mol)

            ph4_gen.generate(lig_mg, rec_mg, ia_ph4, True)  # True = extract environment residues

            core_ph4 = ph4_gen.getCorePharmacophore()
            n_core_ftrs += core_ph4.numFeatures

            if args.mode == 'interaction':
                result = ia_ph4
            else:
                # all probe features, plus the exclusion volumes derived from the
                # receptor environment by the interaction pharmacophore generator
                out_ph4.clear()

                env_ph4 = ph4_gen.getEnvironmentPharmacophore()
                env_donors = env_feature_positions(env_ph4, Pharm.FeatureType.H_BOND_DONOR)
                env_acceptors = env_feature_positions(env_ph4, Pharm.FeatureType.H_BOND_ACCEPTOR)

                for ftr in core_ph4:
                    ftr_type = Pharm.getType(ftr)

                    if ftr_type in (Pharm.FeatureType.H_BOND_DONOR,
                                    Pharm.FeatureType.H_BOND_ACCEPTOR):
                        if not pair_hbond_feature(ftr, ftr_type, env_donors, env_acceptors,
                                                  args.hb_pair_radius):
                            n_unpaired[0] += 1
                            continue

                    out_ph4.addFeature().assign(ftr)

                for ftr in ia_ph4:
                    if Pharm.getType(ftr) == Pharm.FeatureType.EXCLUSION_VOLUME:
                        out_ph4.addFeature().assign(ftr)

                result = out_ph4

            n_out_ftrs += result.numFeatures
            n_kept_ftrs[0] += sum(1 for f in result
                                  if Pharm.getType(f) != Pharm.FeatureType.EXCLUSION_VOLUME)

            if not ph4_writer.write(result):
                sys.exit(f'Error: writing the pharmacophore of probe \'{name}\' failed')

        except SystemExit:
            raise

        except Exception as e:
            sys.exit(f'Error: interaction pharmacophore generation for probe \'{name}\' '
                     f'failed: {str(e)}')

        count += 1

    ph4_writer.close()

    if count == 0:
        sys.exit(f'Error: no probe molecules read from \'{args.probes}\'')

    print(f'- Done! {count} pharmacophores ({args.mode} mode) written to \'{args.output}\'')
    print(f'  probe features generated: {n_core_ftrs}, kept: {n_kept_ftrs[0]}, '
          f'written incl. exclusion volumes: {n_out_ftrs}')

    if n_unpaired[0]:
        print(f'  {n_unpaired[0]} H-bond feature(s) dropped: no complementary receptor '
              f'feature within {args.hb_pair_radius} A')

    if args.mode == 'interaction' and n_core_ftrs and n_kept_ftrs[0] < n_core_ftrs * 0.5:
        print('  ! most probe features had no valid interaction partner -- consider '
              '--mode ligand or --hb-max-length')

    sys.exit(0)


if __name__ == '__main__':
    main()
