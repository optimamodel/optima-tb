import logging
logger = logging.getLogger(__name__)

import numpy as np

from optima_tb.utils import OptimaException, odict
from optima_tb.cascade import __addCharacteristic
from optima_tb.databook import __addCharacteristicData, getEmptyData
# from optima_tb.plotting import extractFlows



def __saveProgressionResults(filehandle, results, pop, progression_from, progression_to, start_year, year_report, starting_pop):
    """
    Saves the results into a formatted csv file: 
    
    """
    indices = np.in1d(results.t_step - start_year, year_report)

    for prog in progression_to:
        filehandle.write("%s --> %s," % (progression_from, prog))

        filehandle.write(",".join(map(str, results.outputs["entry_%s" % prog][pop][indices] / starting_pop)))
        filehandle.write("\n")


def evaluateDiseaseProgression(proj, specified_progressions, specified_populations=None, starting_pop=1e6,
                               year_track=[1., 2., 3., 4., 5.], birth_transit=None,
                               colormappings=None, output_file="ProgressionData.csv", output_fig=""):
    """
    Determines disease progression from certain compartments to other compartments, 
    by (artifically) setting the population of a compartment to a predefined size
    (default size=1e6) and then letting the disease run. Note that birth rate is set
    to zero, so the total population size (including deaths) should be constant.
    
    Inputs:
        proj                          Project object
        specified_populations         list of population labels
        specified_progressions        dictionary with keys = starting compartment 
                                                      values = list of compartments we want percentage values for
        starting_pop                  int representing population size to initialise starting population to
        year_track                    list of floats, describing timepoints (in years) at which to return percentage 
                                      values for i.e. [1.] = value after 1 yr. Note that partial years can be included
                                      but that they must be multiples of dt
        birth_transit                 string representing parameter transition for birth rate parameter.
                                      Default value: None, indicating no such parameter exists
        output_file                   name of output file (string)
    
    
    Example usage:

    proj = Project(name="MyDisease", cascade_path="mycascade.xlsx")
    specified_populations = ['0-4', '5-14', '15-64', '65+', 'HIV 15+', 'Pris']
    specified_progressions = {'sus': ['ltu','acu'],
                              'ltu': ['ltt','acu','ltu'],
                              'acu': ['act','rec']}
    year_track = [0.5,1.,2.] # times that we report on. Note that these can only be multiple of dt
    outputfile = "MyDiseaseProgression.csv"
    birth_rate_parameter = "b_rate"
    
    evaluateDiseaseProgression(proj,
                           specified_progressions=specified_progressions,
                           specified_populations=specified_populations,
                           year_track=year_track,
                           birth_transit='b_rate',
                           output_file=outputfile)
    
    """
    if specified_populations is None:
        specified_populations = proj.parsets[0].pop_labels
        logger.info("Evaluating disease progression for all populations")

    data = proj.data
    dt = proj.settings.tvec_dt
    start_year = proj.settings.tvec_start


    # turn off aging and other transfers
    for transfer_type in data['transfers']:
        data['transfers'][transfer_type] = odict()

    # data['contacts']['into'] = dict()
    # data['contacts']['from'] = dict()

    # turn off all plotting except for population
    for (charac, val) in proj.settings.charac_specs.iteritems():
        val['plot_characteristic'] = 'n'
        val['databook_order'] = -1

    # set births to 0: here, we set the rate to 0
    if birth_transit is not None:
        for (prog, reporters) in specified_progressions.iteritems():
            for pop in specified_populations:
                data['linkpars'][birth_transit][pop]['t'] = [start_year]
                data['linkpars'][birth_transit][pop]['y'] = [0.]

    # remove all references to entry points in project.data characteristics so that we can add our own ...
    for charac, spec in proj.settings.charac_specs.iteritems():
        if 'entry_point' in spec.keys():
            del spec['entry_point']

    # ... and now we add our own entry points and characteristics for reporting.
    # This is where it gets hacky. We have already loaded settings and data, so we have to add
    # to both.
    full_compartment_list = list(set([x for v in specified_progressions.itervalues() for x in v] + specified_progressions.keys()))
    for prog in full_compartment_list: # specified_progressions.iteritems():
        charac_label = 'entry_%s' % prog
        __addCharacteristic(proj.settings, charac_label=charac_label, full_name=charac_label, entry_point=prog, includes=[prog])
        proj.settings.charac_specs[charac_label]['plot_characteristic'] = 'n' # so that we don't plot these new characteristics
        for pop in specified_populations:
            data = __addCharacteristicData(data, charac_label, pop, [start_year], [0.], 'number')



    # Setup simulations, save output and plot:

    output_file_handle = open(output_file, 'w')

    for pop in specified_populations:

        output_file_handle.write("%s," % pop)
        output_file_handle.write(",".join(map(str, year_track)))
        output_file_handle.write("\n")

        for (prog, reporters) in specified_progressions.iteritems():

            parset_name = "prog%s_pop%s" % (prog, pop)
            charac_label = 'entry_%s' % prog
            proj.makeParset(parset_name)
            parset = proj.parsets[parset_name]

            # set up populations and init values
            par = parset.pars['characs'][parset.par_ids['characs'][charac_label]]
            par.y[pop][0] = starting_pop

            # run simulations
            results = proj.runSim(parset_name=parset_name)

            # get outputs
            __saveProgressionResults(output_file_handle, results, pop, prog, reporters, start_year, year_track, starting_pop)

            # save plots
            proj.plotResults(results, plot_observed_data=False, savePlot=True, figName=output_fig + 'DiseaseProgression_compartment_%s' % parset_name, pop_labels=[pop], colormappings=colormappings)

            # reset for the next loop
            par.y[pop][0] = 0.

        output_file_handle.write("\n"*2)

    output_file_handle.close()



def calculateCumulativeDerivatives(results, settings, from_year, to_year,
                                   comp_labels=None, pop_labels=None,
                                   link_labels=None, include_link_not_exclude=True, link_legend=None, sum_total=False,
                                   plot_inflows=True, plot_outflows=True, exclude_transfers=False):
    """
    Calculate the sum of yearly flows to determine total values over a period.
    
    Parameters:
        results        results
        settings       project settings
        from_year, to_year    period over which to calculate the cumulative 
        comp_labels    list of compartments
        pop_labels     list of populations
            other args as per plotting/plotFlows 
    
    Outputs:
        list of summed flows for each year in the period (from_year,to_year)
    
    
    """
    tvec = results.sim_settings['tvec']

    rates, tvecs, labels = extractDerivatives(results=results,
                                              settings=settings,
                                              tvec=tvec,
                                              plot_inflows=plot_inflows,
                                              plot_outflows=plot_outflows,
                                              comp_labels=comp_labels,
                                              # comp_titles = comp_titles,
                                              pop_labels=pop_labels,
                                              # pop_titles = pop_titles,
                                              # plot_pops = pop_labels,
                                              link_labels=link_labels,
                                              include_link_not_exclude=include_link_not_exclude,
                                              link_legend=link_legend,
                                              sum_total=sum_total,
                                              exclude_transfers=exclude_transfers)

    # sum across populations
    yvals = np.array(rates)
    yvals = yvals.sum(axis=0)[0]
    tvals = np.array(tvecs[0])[0]

    # extract years that we need
    idx = (tvals >= from_year) * (tvals <= to_year)

    summed_derivatives = yvals[idx].sum() * settings.tvec_dt

    return summed_derivatives






def getPIDs(results, poplabels):
    """
    Takes in a list of poplabels and returns the corresponding PIDs
    
    TODO: this can be improved and made more efficient by either a better implementation of this
    look up, OR (better yet) by improving the data structures of mpops.
    """
    pids = []
    for poplabel in poplabels:
        for i, pop in enumerate(results.m_pops):
            if pop.label == poplabel:
                pids.append(i)
    return pids


def extractDerivatives(results, settings, tvec, comp_labels=None, comp_titles=None, pop_labels=None,
              link_labels=None, include_link_not_exclude=True, link_legend=None, sum_total=False,
              plot_inflows=True, plot_outflows=True, exclude_transfers=False):
    """
    
    """

    if link_legend is None: link_legend = dict()

    if pop_labels is not None:
        plot_pids = getPIDs(results, pop_labels)
    else:
        plot_pids = range(len(results.m_pops))
        pop_labels = [pop.label for pop in results.m_pops]

    if comp_labels is None:
        logger.error("No compartments have been selected for flow-plots.")
        comp_labels = []
        return


    if comp_titles is not None and len(comp_titles) != len(comp_labels):
        logger.error("Flow-plot failure due to the number of compartment plot titles not matching the number of compartments to analyse.")

    all_rates = []
    all_tvecs = []
    all_labels = []


    for (i, comp_label) in enumerate(comp_labels):

        for (j, pid) in enumerate(plot_pids):

            rates, tvecs, labels = extractFlows(pop_labels=[pid],
                                                    comp_label=comp_label,
                                                    results=results,
                                                    settings=settings,
                                                    tvec=tvec,
                                                    link_labels=link_labels,
                                                    include_link_not_exclude=include_link_not_exclude,
                                                    link_legend=link_legend,
                                                    plot_inflows=plot_inflows,
                                                    plot_outflows=plot_outflows,
                                                    sum_total=sum_total,
                                                    exclude_transfers=exclude_transfers)

            all_rates.append(rates)
            all_tvecs.append(tvecs)
            all_labels.append(labels)

    return all_rates, all_tvecs, all_labels
