import paths
import csv
import os
import utility
import job_results3 as job_results

join = os.path.join

def make_database(res, out_file=paths.results_table):
    def get_data(r):
        d = {}
        d['status'] = r.status
        d['task'] = r.task
        d['reaction'] = r.reaction
        d['enantiomer'] = r.enantiomer
        for R, g in r.substituents:
            d[R] = g
        d['directory'] = r.path
        d['stationary_point'] = r.stationary_point
        d['radical'] = r.radical
        d['inxyz'] = r.data['files'].get('input xyz', None)
        d['outxyz'] = r.data['files'].get('output xyz', None)
        t = int(r.runtime)
        d['runtime'] = f'{t//3600:0>2}:{(t//60)%60:0>2}:{t%60:0>2}'
        d['hash'] = r.hash
        d['step'] = r.step
        return d

    def get_fields():
        return [
          'status',
          'task',
          'reaction',
          'enantiomer',
          'R1', 'R2', 'Rcat', 'Rch',
          'stationary_point',
          'directory',
          'flags',
          'inxyz',
          'outxyz',
          'runtime',
          'step',
          'hash',
          ]


    data = []
    for r in res:
        data.append(get_data(r))

    fields = get_fields()

    #write sep specifier and header
    writer = csv.writer(open(out_file, 'w+', newline=''))
    writer.writerow(['sep=,'])
    writer.writerow([f.upper() for f in fields])

    writer = csv.DictWriter(open(out_file, 'a', newline=''), fields, extrasaction='ignore')
    for d in data:
        writer.writerow(d)





if __name__ == '__main__':

    res = job_results.get_all_results(join(paths.master, 'calculations_test'))
    make_database(res)

    import excel

