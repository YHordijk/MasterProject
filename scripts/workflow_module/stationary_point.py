import os, paths
import workflow_module.results as wfr



class StationaryPoint:
    '''
    Container class for a stationary point
    Each StationaryPoint object holds a number of wfr.Result objects
    and provides an interface to access the results of the various calculations
    '''
    def __init__(self, results):
        self.results = results
        assert len(set(r['stationary_point'] for r in results)) == 1, 'Not all results belong to the same stationary point'

        self.name = results[0]['stationary_point']
        self.jobs = {r['task']:r for r in results}

    def __repr__(self):
        s = f'{self.name} [{", ".join(self.jobs)}]'
        return s

    def __getitem__(self, key):
        if key in ['reaction', 
                   'functional', 
                   'basis', 
                   'numerical_quality', 
                   'frozen_core', 
                   'phase', 
                   'task', 
                   'stationary_point', 
                   'substituents']:
            return self.results[0][key]

        if key == 'name':
            return self['reaction']

        if key == 'gibbs_COSMORS':
            if not 'COSMORS' in self.jobs:
                return
            return self.jobs['COSMORS']['gibbs']
        if key == 'gibbs':
            if not 'FREQ' in self.jobs:
                return
            return self.jobs['FREQ']['gibbs']



def get_all_stationary_points():
    sorte = []

    for i, r1 in enumerate(wfr.results):
        skip = False
        for col in sorte:
            if all(r1 == c for c in col):
                col.append(r1)
                skip = True
                continue
        if skip: continue
        sorte.append([r1])

    sps = []
    for c in sorte:
        sps.append(StationaryPoint(c))

    return sps

stationary_points = get_all_stationary_points()