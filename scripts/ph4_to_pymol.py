#!/usr/bin/env python3
"""
ph4_to_pymol.py -- visualise a LigandScout .pml pharmacophore in PyMOL.

Reads a pharmacophore model written by apo2ph4_generate_ph4.py (or any other
LigandScout-style .pml) and builds a PyMOL scene from it, optionally together
with the receptor and the FTMap probe molecules it was derived from.

Features are rendered as PyMOL objects, not as flat CGO graphics, so they stay
clickable and measurable: every feature is a pseudoatom whose VDW radius is the
feature tolerance, grouped by feature type. H-bond directions are drawn as CGO
arrows pointing at the interaction partner, aromatic features additionally get a
disc in the ring plane.

    # write a session and open it
    python3 ph4_to_pymol.py model.pml -r receptor.pdb -p probes.sdf -o scene.pse
    pymol scene.pse

    # open directly in the PyMOL GUI
    python3 ph4_to_pymol.py model.pml -r receptor.pdb --show

    # just render an image, no GUI needed
    python3 ph4_to_pymol.py model.pml -r receptor.pdb --png view.png
"""

import argparse
import math
import os.path as path
import sys
import xml.etree.ElementTree as ET


# LigandScout-like colours
FEATURE_STYLE = {
    'H':   ((1.00, 0.90, 0.20), 'hydrophobic'),
    'AR':  ((0.40, 0.40, 0.95), 'aromatic'),
    'HBA': ((0.95, 0.25, 0.25), 'H-bond acceptor'),
    'HBD': ((0.25, 0.80, 0.35), 'H-bond donor'),
    'PI':  ((0.20, 0.45, 0.95), 'positive ionizable'),
    'NI':  ((0.95, 0.35, 0.65), 'negative ionizable'),
    'XVOL': ((0.60, 0.60, 0.60), 'exclusion volume'),
}

ARROW_RADIUS = 0.10
ARROW_HEAD_RADIUS = 0.25
ARROW_HEAD_LENGTH = 0.6
DISC_RADIUS = 1.2


def _xyz(elem):
    return (float(elem.get('x3')), float(elem.get('y3')), float(elem.get('z3')))


def _tol(elem, default=1.0):
    try:
        return float(elem.get('tolerance'))
    except (TypeError, ValueError):
        return default


def parse_pml(file_path):
    """Return a list of feature dicts from a LigandScout .pml file."""
    try:
        root = ET.parse(file_path).getroot()
    except ET.ParseError as e:
        sys.exit(f'Error: \'{file_path}\' is not readable as XML: {e}')

    features = []

    for elem in root.iter():
        tag = elem.tag

        if tag == 'point':
            pos = elem.find('position')

            if pos is None:
                continue

            features.append({'type': elem.get('name', '?'), 'pos': _xyz(pos),
                             'tol': _tol(pos)})

        elif tag == 'plane':
            pos = elem.find('position')
            normal = elem.find('normal')

            if pos is None:
                continue

            features.append({'type': elem.get('name', 'AR'), 'pos': _xyz(pos),
                             'tol': _tol(pos),
                             'normal': _xyz(normal) if normal is not None else None})

        elif tag == 'vector':
            origin = elem.find('origin')
            target = elem.find('target')

            if origin is None or target is None:
                continue

            # 'pointsToLigand' marks which end carries the ligand-side feature:
            # true for acceptors (the target is on the ligand), false for donors.
            to_ligand = str(elem.get('pointsToLigand', 'false')).lower() == 'true'
            lig_end, rec_end = (target, origin) if to_ligand else (origin, target)

            features.append({'type': elem.get('name', '?'), 'pos': _xyz(lig_end),
                             'tol': _tol(lig_end), 'partner': _xyz(rec_end)})

        elif tag == 'volume':
            pos = elem.find('position')

            if pos is None:
                continue

            features.append({'type': 'XVOL', 'pos': _xyz(pos), 'tol': _tol(pos)})

    if not features:
        sys.exit(f'Error: no pharmacophore features found in \'{file_path}\'')

    return features


# --------------------------------------------------------------------------
# CGO helpers
# --------------------------------------------------------------------------

def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _norm(v):
    length = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

    return (v[0] / length, v[1] / length, v[2] / length) if length else (0.0, 0.0, 1.0)


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def arrow_cgo(cgo_consts, start, end, color):
    """Cylinder + cone from `start` to `end`."""
    CYLINDER, CONE = cgo_consts

    direction = _norm(_sub(end, start))
    shaft_end = (end[0] - direction[0] * ARROW_HEAD_LENGTH,
                 end[1] - direction[1] * ARROW_HEAD_LENGTH,
                 end[2] - direction[2] * ARROW_HEAD_LENGTH)
    r, g, b = color

    return ([CYLINDER, start[0], start[1], start[2],
             shaft_end[0], shaft_end[1], shaft_end[2], ARROW_RADIUS, r, g, b, r, g, b] +
            [CONE, shaft_end[0], shaft_end[1], shaft_end[2], end[0], end[1], end[2],
             ARROW_HEAD_RADIUS, 0.0, r, g, b, r, g, b, 1.0, 1.0])


def disc_cgo(cgo_consts_disc, center, normal, color, radius=DISC_RADIUS, segments=32):
    """Filled disc centred at `center`, perpendicular to `normal`."""
    BEGIN, END, TRIANGLE_FAN, COLOR, NORMAL, VERTEX = cgo_consts_disc

    n = _norm(normal)
    helper = (1.0, 0.0, 0.0) if abs(n[0]) < 0.9 else (0.0, 1.0, 0.0)
    u = _norm(_cross(n, helper))
    v = _cross(n, u)
    r, g, b = color

    cgo = [BEGIN, TRIANGLE_FAN, COLOR, r, g, b, NORMAL, n[0], n[1], n[2],
           VERTEX, center[0], center[1], center[2]]

    for i in range(segments + 1):
        angle = 2.0 * math.pi * i / segments
        ca, sa = math.cos(angle) * radius, math.sin(angle) * radius
        cgo += [VERTEX,
                center[0] + u[0] * ca + v[0] * sa,
                center[1] + u[1] * ca + v[1] * sa,
                center[2] + u[2] * ca + v[2] * sa]

    cgo.append(END)

    return cgo


# --------------------------------------------------------------------------

def build_scene(cmd, features, args):
    from pymol.cgo import (CYLINDER, CONE, BEGIN, END, TRIANGLE_FAN, COLOR,
                           NORMAL, VERTEX)

    cmd.set('retain_order', 1)
    cmd.set('pse_export_version', 1.8)

    if args.receptor:
        cmd.load(args.receptor, 'receptor')
        cmd.hide('everything', 'receptor')
        cmd.show('cartoon', 'receptor')
        cmd.color('grey80', 'receptor')
        cmd.set('cartoon_transparency', 0.75, 'receptor')
        cmd.set('stick_radius', 0.12, 'receptor')

    if args.probes:
        cmd.load(args.probes, 'probes')
        cmd.hide('everything', 'probes')
        cmd.show('lines', 'probes')
        cmd.color('grey50', 'probes')
        cmd.set('line_width', 1.5)
        # a multi-molecule SD-file loads as one object with one state per molecule,
        # and PyMOL draws only the current state -- show the whole cluster at once
        cmd.set('all_states', 1, 'probes')

        if args.hide_probes:
            # kept in the object panel so it can be toggled back on with one click
            cmd.disable('probes')

    made = []          # every object, for grouping
    sphere_objs = []   # pseudoatom objects only (CGO objects hold no atoms)
    counts = {}

    for i, ftr in enumerate(features):
        ftr_type = ftr['type']
        color, _ = FEATURE_STYLE.get(ftr_type, ((0.8, 0.8, 0.8), 'unknown'))
        counts[ftr_type] = counts.get(ftr_type, 0) + 1

        obj = f'ph4_{ftr_type}'
        color_name = f'ph4col_{ftr_type}'
        cmd.set_color(color_name, list(color))

        cmd.pseudoatom(obj, pos=list(ftr['pos']), vdw=ftr['tol'], elem='PS',
                       name=f'{ftr_type}{counts[ftr_type]}', resn=ftr_type,
                       resi=counts[ftr_type], color=color_name)

        if obj not in sphere_objs:
            sphere_objs.append(obj)
            made.append(obj)

        if ftr.get('partner'):
            arrow = arrow_cgo((CYLINDER, CONE), ftr['pos'], ftr['partner'], color)
            name = f'ph4_{ftr_type}_dir_{i}'
            cmd.load_cgo(arrow, name)
            made.append(name)

        if ftr.get('normal'):
            for sign in (1.0, -1.0):
                normal = tuple(sign * c for c in ftr['normal'])
                name = f'ph4_AR_plane_{i}_{"a" if sign > 0 else "b"}'
                cmd.load_cgo(disc_cgo((BEGIN, END, TRIANGLE_FAN, COLOR, NORMAL, VERTEX),
                                      ftr['pos'], normal, color), name)
                made.append(name)

    for obj in sphere_objs:
        cmd.hide('everything', obj)
        cmd.show('spheres', obj)
        cmd.set('sphere_transparency',
                args.xvol_transparency if obj == 'ph4_XVOL' else args.transparency,
                obj)

    if 'ph4_XVOL' in sphere_objs and args.xvol_as_mesh:
        cmd.hide('spheres', 'ph4_XVOL')
        cmd.show('mesh', 'ph4_XVOL')

    cmd.group('pharmacophore', ' '.join(made))

    # focus on the model, keeping the surrounding residues visible
    cmd.orient('pharmacophore')
    cmd.zoom('pharmacophore', args.zoom_buffer)

    if args.receptor and args.show_sidechains:
        cmd.select('_env', f'byres (receptor within {args.env_radius} of pharmacophore)')
        cmd.show('sticks', '_env and sidechain')
        cmd.color('grey55', '_env and elem C')
        cmd.delete('_env')

    # hydrogens add nothing here and obscure the features
    for obj in ('receptor', 'probes'):
        if obj in cmd.get_object_list():
            cmd.hide('everything', f'{obj} and hydro')

    cmd.set('bg_rgb', [1, 1, 1])
    cmd.set('ray_opaque_background', 0)
    cmd.set('depth_cue', 0)
    cmd.set('specular', 0.2)

    return counts


def parse_args():
    p = argparse.ArgumentParser(
        description='Visualise a LigandScout .pml pharmacophore in PyMOL.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('pml', help='pharmacophore model (.pml)')
    p.add_argument('-r', '--receptor', help='receptor structure to show alongside')
    p.add_argument('-p', '--probes', help='probe/ligand molecules to show alongside')
    p.add_argument('-o', '--output', help='write a PyMOL session (.pse) here')
    p.add_argument('--png', help='render an image to this path and exit')
    p.add_argument('--show', action='store_true', help='open the PyMOL GUI')
    p.add_argument('--width', type=int, default=1400, help='image width')
    p.add_argument('--height', type=int, default=1050, help='image height')
    p.add_argument('--transparency', type=float, default=0.45,
                   help='sphere transparency of the pharmacophore features')
    p.add_argument('--xvol-transparency', type=float, default=0.8,
                   help='sphere transparency of the exclusion volumes')
    p.add_argument('--xvol-as-mesh', action='store_true',
                   help='draw exclusion volumes as a mesh instead of spheres')
    p.add_argument('--hide-probes', action='store_true',
                   help='load the probes but start with them switched off, so the scene '
                        'shows only the pharmacophore model (toggle them back on by '
                        'clicking "probes" in the PyMOL object panel)')
    p.add_argument('--zoom-buffer', type=float, default=7.0,
                   help='padding (A) around the model when framing the view')
    p.add_argument('--env-radius', type=float, default=5.0,
                   help='show receptor side chains within this distance of the model')
    p.add_argument('--show-sidechains', action='store_true', default=True,
                   help='show receptor side chains within 5 A of the model')
    p.add_argument('--no-sidechains', dest='show_sidechains', action='store_false')

    args = p.parse_args()

    if not (args.output or args.png or args.show):
        p.error('nothing to do: pass --show, --png and/or -o')

    for f in (args.pml, args.receptor, args.probes):
        if f and not path.isfile(f):
            sys.exit(f'Error: \'{f}\' not found')

    return args


def main():
    args = parse_args()
    features = parse_pml(args.pml)

    if args.show:
        import pymol
        pymol.finish_launching(['pymol', '-q'])
        from pymol import cmd
    else:
        import pymol2
        session = pymol2.PyMOL()
        session.start()
        cmd = session.cmd

    counts = build_scene(cmd, features, args)

    print(f'- {len(features)} features from \'{args.pml}\':')

    for ftr_type, n in sorted(counts.items()):
        print(f'    {n:3d} x {ftr_type:4s} ({FEATURE_STYLE.get(ftr_type, (None, "?"))[1]})')

    if args.output:
        cmd.save(args.output)
        print(f'- session written to \'{args.output}\'  ->  pymol {args.output}')

    if args.png:
        cmd.set('ray_trace_mode', 0)
        cmd.png(path.abspath(args.png), width=args.width, height=args.height,
                dpi=150, ray=1)
        print(f'- image written to \'{args.png}\'')

    if args.show:
        print('- PyMOL GUI launched')
    else:
        session.stop()


if __name__ == '__main__':
    main()
