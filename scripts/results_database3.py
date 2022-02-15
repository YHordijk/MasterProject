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
        d['basis'] = r.basis
        d['functional'] = r.functional
        d['quality'] = r.numerical_quality
        d['reaction'] = r.reaction
        d['enantiomer'] = r.enantiomer
        for R, g in r.substituents.items():
            d[R] = g
        d['res_dir'] = r.path
        d['calc_dir'] = r.calc_path
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
          'enantiomer',] + list(sorted(all_subs)) + [
          'stationary_point',
          'res_dir',
          'calc_dir',
          'inxyz',
          'outxyz',
          'runtime',
          'step',
          'functional', 'basis', 'quality',
          'hash',
          ]


    data = []
    all_subs = {s for r in res for s in r.substituents.keys()}
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
    res = job_results.get_all_results(paths.calculations, regenerate_all=False)
    make_database(res)

    import excel