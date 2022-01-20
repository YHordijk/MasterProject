import scm.plams as plams
import paths


def load_results(path):
    job = plams.AMSJob()
    job.load_external(path)
    job.path = path
    res = job.results
    res.collect()
    print(dir(res))
    # print(res.files)
    print(job.status)
    print(res.get_main_molecule())

    

    return res



def get_status(res):
    term = result.readrkf('General', 'termination status', 'ams')
    if term == 'NORMAL TERMINATION':
        return True
    elif 'NORMAL TERMINATION' in term:
        return 'WARNING'
    return False



res = load_results(r"D:\Users\Yuman\Desktop\MasterProject\calculations\no_catalyst_H_H\2.no_catalyst_H_H.Rrad")
print(get_status(res))