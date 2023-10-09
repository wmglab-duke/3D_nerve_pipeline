import json
import os
import sys

import numpy as np
from shapely.geometry import Point


class FascicleConnectivityMap:
    def __init__(self, mapslides):
        self.map_data = None
        self.slides = mapslides
        self.branch_count = 0
        self.branch_inds = []
        self.merge_inds = []
        self.split_inds = []

    def generate_map(self, ignore_dead_end=False, primaries=False):
        # flag: some more work will need to be done using contours overlap rather than centroids
        # flag: currently cannot handle multiple inners
        mapslides = self.slides
        connections = []
        hov_fasc = {}
        active_ids = []
        u_id = 0  # noqa: SIM113
        for fascicle in mapslides[0].fascicles:
            # Initialize the first layer of fascicules
            u_id += 1
            hov_fasc[u_id] = {
                'points': [fascicle.centroid() + (0,)],
                'shapes': [fascicle.inners[0].polygon()],
                'start': True,
            }
            active_ids.append(u_id)
        for i in range(1, len(mapslides)):
            fascicles = mapslides[i].fascicles
            temp_ids = range(0, len(fascicles))
            old_contains_new = {}
            new_contains_old = {}
            for d in active_ids:
                container = [
                    x
                    for x in temp_ids
                    if hov_fasc[d]['shapes'][-1].contains(
                        Point(fascicles[x].inners[0].polygon().representative_point())
                    )  # TODO decide whether to keep change
                ]
                old_contains_new[d] = container
            for x in temp_ids:
                container = [
                    d
                    for d in active_ids
                    if fascicles[x]
                    .inners[0]
                    .polygon()
                    .contains(hov_fasc[d]['shapes'][-1].representative_point())  # TODO decide whether to keep change
                ]
                new_contains_old[x] = container
            for d, cont in old_contains_new.copy().items():
                if len(cont) == 1 and new_contains_old[cont[0]] == [
                    d
                ]:  # checks to see if the fascicle remained between slices
                    hov_fasc[d]['points'].append(fascicles[cont[0]].centroid() + (i,))
                    hov_fasc[d]['shapes'].append(fascicles[cont[0]].inners[0].polygon())
                    old_contains_new.pop(d)
                    new_contains_old.pop(cont[0])
            if len(new_contains_old) > 0 or len(old_contains_new) > 0:
                for d, cont in old_contains_new.copy().items():
                    if len(cont) > 1:  # detects a fascicle split
                        # check for a merge and split at the same time
                        if not np.any([len(vals) > 1 and d in vals for vals in new_contains_old.values()]):
                            self.branch_count += 1
                            self.branch_inds.append(i)
                            self.split_inds.append(i)
                            active_ids.remove(d)
                            fro = [d]
                            to = []
                            for x in cont:
                                u_id += 1
                                hov_fasc[u_id] = {
                                    'points': [fascicles[x].centroid() + (i,)],
                                    'shapes': [fascicles[x].inners[0].polygon()],
                                }
                                active_ids.append(u_id)
                                new_contains_old.pop(x)
                                to.append(u_id)
                            old_contains_new.pop(d)
                            connections.append([fro, to])
                        else:
                            print("newspecialfix")
                            merges = [len(vals) - 1 for vals in new_contains_old.values() if len(vals) > 1]
                            splits = [len(vals) - 1 for vals in old_contains_new.values() if len(vals) > 1]
                            if sum(merges) != 1 or sum(splits) != 1:
                                raise NotImplementedError(
                                    "This code handles a very specific edge case\
                                        and has not been tested for other \
                                            mergsplit multi cases"
                                )
                            self.branch_count += sum(merges)
                            self.branch_count += sum(splits)
                            self.branch_inds.append(i)
                            self.split_inds.append(i)
                            self.merge_inds.append(i)
                            donzo = [value for value in new_contains_old.values() if d in value][0]
                            ocn_matches = {}
                            these_tos = []
                            for x in cont:
                                u_id += 1
                                ocn_matches[x] = u_id
                                hov_fasc[u_id] = {
                                    'points': [fascicles[x].centroid() + (i,)],
                                    'shapes': [fascicles[x].inners[0].polygon()],
                                }
                                active_ids.append(u_id)
                                fro = []
                                to = []
                                fro.extend(new_contains_old[x])
                                fro.extend([key for key, value in old_contains_new.items() if x in value])
                                to.extend([u_id])
                                if len(set(fro)) == 1:
                                    olds = [str(num) + '_old' for num in old_contains_new[fro[0]]]
                                    to.extend(olds)
                                to = list(set(to))
                                fro = list(set(fro))
                                these_tos.append(to)
                                connections.append([fro, to])
                            for fnum in donzo:
                                active_ids.remove(fnum)
                                old_contains_new.pop(fnum)
                            for tonum in cont:
                                new_contains_old.pop(tonum)
                            for to in these_tos:
                                for val in to:
                                    if type(val) is str:
                                        newnum = int(val.split('_')[0])
                                        to[to.index(val)] = ocn_matches[newnum]
                for x, cont in new_contains_old.copy().items():
                    if len(cont) > 1:  # detects a fascicle merge
                        self.branch_count += 1
                        self.merge_inds.append(i)
                        self.branch_inds.append(i)
                        fro = []
                        for d in cont:
                            active_ids.remove(d)
                            fro.append(d)
                        u_id += 1
                        hov_fasc[u_id] = {
                            'points': [fascicles[x].centroid() + (i,)],
                            'shapes': [fascicles[x].inners[0].polygon()],
                        }
                        active_ids.append(u_id)
                        new_contains_old.pop(x)
                        to = [u_id]
                        for d in cont:
                            old_contains_new.pop(d)
                        connections.append([fro, to])
            if len(new_contains_old) > 0 or len(old_contains_new) > 0:
                if len(new_contains_old) > 0:  # detects fascicles that shifted too much but are the same
                    for x, cont in new_contains_old.copy().items():
                        if len(cont) != 1:
                            continue
                        hov_fasc[cont[0]]['points'].append(fascicles[x].centroid() + (i,))
                        hov_fasc[cont[0]]['shapes'].append(fascicles[x].inners[0].polygon())
                        new_contains_old.pop(x)
                        old_contains_new.pop(cont[0])
                if len(old_contains_new) > 0:  # detects fascicles that shifted too much but are the same
                    for d, cont in old_contains_new.copy().items():
                        if len(cont) != 1:
                            continue
                        hov_fasc[d]['points'].append(fascicles[cont[0]].centroid() + (i,))
                        hov_fasc[d]['shapes'].append(fascicles[cont[0]].inners[0].polygon())
                        old_contains_new.pop(d)
                        new_contains_old.pop(cont[0])
                for x, cont in new_contains_old.copy().items():
                    if len(cont) == 0:  # checks for appearing fascicles
                        if not ignore_dead_end:
                            raise RuntimeError("appearing fascicle")
                        u_id += 1
                        active_ids.append(u_id)
                        hov_fasc[u_id] = {
                            'points': [fascicles[x].centroid() + (i,)],
                            'shapes': [fascicles[x].inners[0].polygon()],
                        }
                        new_contains_old.pop(x)
                for d, cont in old_contains_new.copy().items():
                    if len(cont) == 0:  # checks for disappearing fascicles
                        if not ignore_dead_end:
                            # TODO: fix to handle split and merge in same slice
                            raise RuntimeError(f"disappearing fascicle index {i}")
                        active_ids.remove(d)
                        old_contains_new.pop(d)
            if len(new_contains_old) > 0 and len(old_contains_new) > 0:
                sys.exit('Lost a fascicle in connectivity map')
            if len(new_contains_old) > 0 or len(old_contains_new) > 0:
                sys.exit('Lost a fascicle in connectivity map')
        for uid in active_ids:
            hov_fasc[uid]['end'] = True
        # get mean area
        for uid, info in hov_fasc.items():
            hov_fasc[uid]['from'] = []
            hov_fasc[uid]['to'] = []
            hov_fasc[uid]['mean_area'] = np.mean([x.area for x in info['shapes']])
            if 'start' not in info:
                hov_fasc[uid]['start'] = False
            if 'end' not in info:
                hov_fasc[uid]['end'] = False

        # get all connections and all primary connections
        for connection in connections:
            for fro in connection[0]:
                for to in connection[1]:
                    hov_fasc[fro]['to'].append(to)
                    hov_fasc[to]['from'].append(fro)
        if primaries:
            for uid in hov_fasc.keys():
                if len(hov_fasc[uid]['from']) > 1:
                    maxfrom = np.argmax([hov_fasc[fro]['mean_area'] for fro in hov_fasc[uid]['from']])
                    hov_fasc[uid]['primary_from'] = hov_fasc[uid]['from'][maxfrom]
                elif len(hov_fasc[uid]['from']) == 1:
                    hov_fasc[uid]['primary_from'] = hov_fasc[uid]['from'][0]
                else:
                    hov_fasc[uid]['primary_from'] = None
                if len(hov_fasc[uid]['to']) > 1:
                    maxto = np.argmax([hov_fasc[to]['mean_area'] for to in hov_fasc[uid]['to']])
                    hov_fasc[uid]['primary_to'] = hov_fasc[uid]['to'][maxto]
                elif len(hov_fasc[uid]['to']) == 1:
                    hov_fasc[uid]['primary_to'] = hov_fasc[uid]['to'][0]
                else:
                    hov_fasc[uid]['primary_to'] = None

            for uid, info in hov_fasc.items():
                search_id = uid
                allpoint = hov_fasc[search_id]['points'].copy()
                while hov_fasc[search_id]['primary_from'] is not None:
                    search_id = hov_fasc[search_id]['primary_from']
                    allpoint.extend(hov_fasc[search_id]['points'])
                search_id = uid
                while hov_fasc[search_id]['primary_to'] is not None:
                    search_id = hov_fasc[search_id]['primary_to']
                    allpoint.extend(hov_fasc[search_id]['points'])
                superpoints = np.array(allpoint.copy())
                hov_fasc[uid]['superpoints'] = superpoints[superpoints[:, 2].argsort()]
                if len(allpoint) < len(mapslides):
                    hov_fasc[uid]['full_length'] = True
                else:
                    hov_fasc[uid]['full_length'] = False
                if (len(info['to']) == 0 and info['end'] is False) or (
                    len(info['from']) == 0 and info['start'] is False
                ):
                    hov_fasc[uid]['dead_end'] = True
                    if not ignore_dead_end:
                        raise RuntimeError("Dead end fascicle found.")
                else:
                    hov_fasc[uid]['dead_end'] = False

        self.map_data = hov_fasc

    def get_simplemap(self):
        simplemap = {}
        for key, value in self.map_data.items():
            simplemap[key] = {
                'from': len(value['from']),
                'to': len(value['to']),
                'startz': int(np.amin(np.array(value['points'])[:, 2])),
                'endz': int(np.amax(np.array(value['points'])[:, 2])),
            }
        return simplemap

    # TODO add to sample as method
    def save_images(self, path, dims, separate=False, print_ids=False, resize_factor=1, id_font=40):
        # formats points and prints the fascicle ids on images
        for i in range(len(self.slides)):
            sl = self.slides[i]
            ids = []
            if separate:
                paths = {
                    'n': path + f'/n/{i}.tif',
                    'i': path + f'/i/{i}.tif',
                    'p': path + f'/p/{i}.tif',
                }
                for impath in paths.values():
                    os.makedirs(os.path.dirname(impath), exist_ok=True)
                sl.saveimg(
                    paths,
                    dims,
                    separate=True,
                    ids=ids,
                    resize_factor=resize_factor,
                    id_font=id_font,
                )
            else:
                os.makedirs(path, exist_ok=True)
                sl.saveimg(
                    path + f'/{i}.tif',
                    dims,
                    separate=False,
                    ids=ids,
                    resize_factor=resize_factor,
                    id_font=id_font,
                )

    def nervestats(self, outpath):
        areas = [(10**-6) * (s.nerve.area() / np.pi) ** 0.5 for s in self.slides]
        areadict = {'min': min(areas), 'mean': np.mean(areas), 'max': max(areas)}
        with open(outpath, 'w') as f:
            json.dump(areadict, f, indent=2)
