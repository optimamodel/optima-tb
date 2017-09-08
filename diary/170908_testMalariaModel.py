import os
import pylab as pl
import numpy as np
import csv
from optima_tb.project import Project
from optima_tb.asd import asd
from optima_tb.utils import odict


cascade_path = os.path.abspath('.//cascade-m-170908.xlsx')
databook_path = os.path.abspath('.//databook-m-170908.xlsx')


#%% Function for inspecting result outputs
def plotResults(res_dict, dt, t_start, t_end, labels, den_labels = None, pop_labels = None):
    
    figr, axr = pl.subplots()
    
    for result_label in res_dict.keys():
        
        result = res_dict[result_label]
        try: del vals
        except: pass
        for label in labels:
            try: vals += result.getValuesAt(label, t_start, t_end, pop_labels=pop_labels)[0]
            except: vals = result.getValuesAt(label, t_start, t_end, pop_labels=pop_labels)[0]
            
        try: del den_vals
        except: pass
        if not den_labels is None:
            for label in den_labels:
                try: den_vals += result.getValuesAt(label, t_start, t_end,pop_labels=pop_labels)[0]
                except: den_vals = result.getValuesAt(label, t_start, t_end,pop_labels=pop_labels)[0]
            
        try: axr.plot(np.arange(t_start, t_end, dt), vals/den_vals, '-', label=result_label)
        except: axr.plot(np.arange(t_start, t_end, dt), vals, '-', label=result_label)
    
    axr.legend()
    
    return figr


# wrapper function for plotResults
def visualize(proj, res_dict):
    outputs_to_plot = ['t_hex']

    start = 2001.
    stop = 3825.
    for pop_label in proj.data['pops']['label_names'].keys():
        for output in outputs_to_plot:
            fig = plotResults(res_dict=res_dict, dt=proj.settings.tvec_dt, t_start=start, t_end=stop, labels=[output],
                              pop_labels=[pop_label])
            try:
                fig.suptitle(proj.settings.charac_specs[output]['name'])
            except:
                fig.suptitle(output)

    pl.show()


# write parameters/characteristics specified in stats into a csv-file
def write_results(res, filename='foo.csv'):
    stats = ['t_hinf', 'h_prevx', 't_hdi', 't_htreat', 'm_prevx', 'num_hsus']
    pops = ['gp', 'preg', 'child']

    start = 2000
    stop = 3825
    step = 365
    N = int((stop - start) / step)

    with open(filename, 'wb') as f:
        csvw = csv.writer(f)

        for s in stats:
            csvw.writerow([s])
            csvw.writerow(range(start, stop, step))

            for p in pops:
                vals, _ = res.getValuesAt(s, start - step, stop, pop_labels=p)

                if s.startswith('t_'):
                    tmp = np.zeros(N)
                    for k in range(N):
                        tmp[k] = sum(vals[k*step:(k+1)*step])
                    csvw.writerow(tmp)
                elif s in 'h_prevx' or s in 'm_prevx':
                    tmp = np.zeros(N)
                    for k in range(N):
                        tmp[k] = np.mean(vals[k * step:(k + 1) * step])
                    csvw.writerow(tmp)
                else:
                    csvw.writerow(vals[::step])
            csvw.writerow(['\n'])


# Forces the result of a function evaluation to be interpreted as number not as fraction
def interpret_value_as_number(proj, config='default'):
    # force the transitions from the birth compartments to be interpreted as 'number', not as 'fraction'
    for pop in proj.parsets[config].pop_labels:
        for label in ['num_hb', 'num_mb', 'num_htx']:
            idx = proj.parsets[config].par_ids['cascade'][label]
            proj.parsets[config].pars['cascade'][idx].y_format[pop] = 'number'
    return proj


def setup_project():
    proj = Project(name='malaria', cascade_path=cascade_path, validation_level='warn')
    proj.loadSpreadsheet(databook_path=databook_path)

    dt = proj.settings.tvec_dt = 1. / 1.
    t_start = proj.settings.tvec_start = -1650. # 10 years to eliminate initial transient
    proj.settings.tvec_observed_end = 3825.
    t_end = proj.settings.tvec_end = 3825.

    print('Simulating from ' + str(t_start) + ' to ' + str(t_end) + ' in ' + str(dt) + ' increments')

    proj.makeParset()
    proj = interpret_value_as_number(proj)

    return proj

# convert a dictionary to a list in the order given by pops and pars; inverse of list2dict
def dict2list(D, pops, pars):
    L = []
    for pop in pops:
        for par in pars:
            L.append(D[pop][par])
    return L

# convert list to a dictionary; inverse of dict2list
def list2dict(L,pops, pars):
    d = {}
    l = L.tolist()
    for pop in pops:
        d[pop] = {}
        for par in pars:
            d[pop][par] = l.pop(0)
    return d

# determines the objective value during calibration.
# proj and result contain the data (referenced by the name pname). ref_id and sim_id references the parameter/value
# to be compared. func specifies a function used for comparison, e.g. the mean or sum of all values in a year
def comp_score(proj, pname, result, ref_id, sim_id, func):
    if ref_id in 'tau':
        T = 'cascade'
    else:
        T = 'characs'

    pops = proj.parsets[pname].pop_labels
    data = proj.parsets[pname].pars[T][proj.parsets[pname].par_ids[T][ref_id]]

    start = 2000
    stop = 3825
    step = 365
    N = int((stop - start) / step) + 1

    err = []

    simdata = {}
    for pop in pops:
        simdata[pop] = result.getValuesAt(sim_id, start - step, stop, pop_labels=pop)[0]
        tmp = np.zeros(N)
        for k in range(N):
            # determine the simulation result value
            tmp[k] = func(simdata[pop][k * step:(k + 1) * step])
        simdata[pop] = tmp

        assert (len(simdata[pop]) == len(data.y[pop]) or len(data.y[pop] == 1))

        # compute relative error of the simulation result with the reference value
        err.append(np.mean(abs(simdata[pop] - data.y[pop]) / np.mean(data.y[pop])))

    return err


# objective function for the calibration which is passed to asd
def calibration_objective(vals, proj, pops, pars):
    d = list2dict(vals, pops, pars)

    pname = 'calibration'
    proj.makeManualCalibration(parset_name=pname, rate_dict=d)
    proj = interpret_value_as_number(proj, pname)

    try:
        result = proj.runSim(parset_name=pname)
    except: # if anything goes wrong, discard result. NOTE: requires 'ERROR' logging-mode
        return np.inf

    # compute objective value
    score = []
    score.extend(comp_score(proj, pname, result, ref_id='h_inci', sim_id='t_hinf', func=sum))
    score.extend(comp_score(proj, pname, result, ref_id='h_prevx', sim_id='h_prevx', func=np.mean))
    score.extend(comp_score(proj, pname, result, ref_id='num_hmd', sim_id='t_hdi', func=sum))
    score.extend(comp_score(proj, pname, result, ref_id='tau', sim_id='t_htreat', func=np.mean))

    del proj.parsets[pname]

    return np.mean(score)


def calibrate(proj):
    pops = ['gp', 'preg', 'child']
    pars = ['rate_hdi', 'rho', 'phi', 'f']

    # parameter upper and lower bounds, empirically determined
    ubs = [0.02, 900.0, 0.02, 1.0,      0.01, 350.0, 0.025, 1.0,         0.08, 300.0, 0.028, 1.0]
    lbs = [0.005, 600.0, 0.005, 0.0,      0.001, 100.0, 0.003, 0.0,         0.055, 100.0, 0.010, 0.0]
    vals = (np.array(lbs) + np.array(np.random.rand(len(ubs))) * np.array(np.array(ubs)-np.array(lbs))).tolist()

    assert((len(vals) == len(ubs)) and (len(ubs) == len(lbs)))

    # asd optimisation
    args = {'proj': proj, 'pops': pops, 'pars': pars}
    opt_pars, score, exit_code = asd(calibration_objective, vals, xmin=lbs, xmax=ubs, args=args, maxiters=1)

    print('Optimised score: ' + str(score[-1]))
    print(list2dict(np.array(opt_pars), pops, pars))
    return list2dict(np.array(opt_pars), pops, pars), score[-1] # return configuration and score of best result


# calibrate model in several runs; in each run, a csv-file is written with the simulated values (indexed by the number
# of run) and a txt-file with the configuration whose file name includes the score and index. So sorting the text-files
# alphabetically is equivalent to best to worst ordering
def exec_calib():
    proj = setup_project()

    for n in range(1):
        calibrated_pars, score = calibrate(proj)
        proj.makeManualCalibration(parset_name='calib', rate_dict=calibrated_pars)
        proj = interpret_value_as_number(proj, 'calib')

        res_dict = odict()
        run_cal = proj.runSim(parset_name='calib')
        res_dict['Standard Run'] = run_cal

        path = './/'
        write_results(run_cal, path + 'calib' + format(n, '04') + '.csv')
        with open(path + 'pars_' + format(score, '08') + '_' + format(n, '04') + '.txt', 'w') as f:
            f.write(str(calibrated_pars))

if __name__ == "__main__":
    exec_calib()
