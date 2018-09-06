from collections import defaultdict
import uuid
import optima_tb.settings as project_settings
from optima_tb.utils import OptimaException
import numpy as np
import logging
from optima_tb.parsing import parse_function

logger = logging.getLogger(__name__)

class ModelProgramSet(object):
    # Contains some programs

    def __init__(self,progset,pops):
        # pops is model.pops

        self.uid = uuid.uuid4()
        self.programs = list()
        self.progset = progset # A complete ProgramSet instance
        self.pars = set()
        self.dt = 0 # Timestep set when calling update_cache and same timestep used for scaling in compute_pars
        self.tval_cache = None # Time values are cached and passed to ModelPrograms for program-specific activation/deactivation

        # Construct the ModelProgram objects
        for prog in progset.progs: # For each Program, instantiate a ModelProgram

            # Make a list of all impact Parameter objects
            impact_pars = []
            for pop in pops:
                if pop.label in prog.target_pops:
                    impact_pars += [par for par in pop.pars if par.label in prog.target_pars]

            # Now populate the impact groups with those parameters
            impact_groups = defaultdict(list)
            for par in impact_pars:
                if 'group' in progset.progtype_specs[prog.prog_type]['impact_pars'][par.label]: # If this parameter has an impact group
                    impact_groups[progset.progtype_specs[prog.prog_type]['impact_pars'][par.label]['group']].append(par)
                else:
                    impact_groups[par.label].append(par) # The parameter becomes a member of its own impact group

            is_fraction = prog.cov_format != 'number' # programs leave this as None otherwise (because it's a fraction and not probability?)

            self.programs.append(ModelProgram(prog.label,is_fraction,impact_pars,impact_groups,prog.active_times))
            self.pars.update(impact_pars) # Update the set of all Parameter objects being reached by adding any Parameters that aren't present
            for par in self.pars:
                par.program_overwrite = True # Mark the program as being overwritten by a program

        self.par_by_id = {x.uid:x for x in self.pars} # Map parameter ID to parameter object
        self.constraints = dict() # par_uid:(lower,upper)

    def load_constraints(self,constraints=None):
        self.constraints = dict()
        if constraints is None:
            for par in self.pars:
                self.constraints[par.uid] = [0.0,np.inf]
        else:
            for par_label,d in constraints['impacts'].items():
                for par in self.pars:
                    if par.label == par_label:
                        self.constraints[par.uid] = d['vals']
        return

    def unlink(self):
        self.pars = list(self.pars)
        for i in xrange(0,len(self.pars)):
            self.pars[i] = self.pars[i].uid
        for prog in self.programs:
            prog.unlink()

    def relink(self,objs):
        self.pars = list(self.pars)
        for i in xrange(0,len(self.pars)):
            self.pars[i] = objs[self.pars[i]]
        self.pars = set(self.pars)
        for prog in self.programs:
            prog.relink(objs)

    def get_alloc(self,t,dt,sim_settings):
        # Extract the allocation (spending value at each point in time for all active programs)
        # based on sim_settings
        #
        # This returns alloc (dict with prog_label:[spending vals]) and alloc_is_coverage flag
        # And corresponding time points
        # Note that spending vals could be a single item, or could have the same shape as t
        
        # Take in a set of times and 
        start_year = sim_settings['progs_start']
        init_alloc = sim_settings['init_alloc']
        alloc_is_coverage = sim_settings['alloc_is_coverage']

        alloc = dict()

        for prog in self.progset.progs: # Iterate over classic Programs, not ModelPrograms

            # Store budgets/coverages for programs that are initially allocated.
            if prog.label in init_alloc:
                spending = init_alloc[prog.label]

                # If ramp constraints are active, stored cost and coverage needs to be a fully time-dependent array corresponding to time points.
                if not isinstance(spending,np.ndarray) and 'constraints' in sim_settings and sim_settings['constraints'] is not None and 'max_yearly_change' in sim_settings['constraints'] and prog.label in sim_settings['constraints']['max_yearly_change']:
                    # Todo - alloc is coverage below does not seem to get used?
                    if alloc_is_coverage:
                        default = prog.getCoverage(budget=prog.getDefaultBudget(year=start_year))
                    else:
                        default = prog.getDefaultBudget(year=start_year)
                    default = prog.getDefaultBudget(year=start_year)
                    
                    

                    if np.abs(spending - default) > project_settings.TOLERANCE:
                        if default == 0.0:
                            default = project_settings.TOLERANCE
                        
                        spending_def = t * 0.0 + default
                        spending_new = t * 0.0 + spending
                        try:
                            eps = sim_settings['constraints']['max_yearly_change'][prog.label]['val']
                        except:
                            raise OptimaException('ERROR: A maximum yearly change constraint was passed to the model for "%s" but had no value associated with it.' % prog.label)
                        if 'rel' in sim_settings['constraints']['max_yearly_change'][prog.label] and sim_settings['constraints']['max_yearly_change'][prog.label]['rel'] is True:
                            eps *= default
                        if np.isnan(eps):
                            eps = np.inf  # Eps likely becomes a nan if it was infinity multiplied by zero.
                        if np.abs(eps * dt) < np.abs(spending - default):
                            if np.abs(eps) < project_settings.TOLERANCE:
                                raise OptimaException('ERROR: The change in budget for ramp-constrained "%s" is effectively zero. Model will not continue running; change in program funding would be negligible.' % prog.label)
                            spending_ramp = default + (t - start_year) * eps * np.sign(spending - default)
                            if spending >= default:
                                spending_ramp = np.minimum(spending_ramp, spending_new)
                            else:
                                spending_ramp = np.maximum(spending_ramp, spending_new)
                            spending = spending_def * (t < start_year) + spending_ramp * (t >= start_year)

                # Todo - alloc_is_coverage is mentioned above? Does this statement execute correctly, or should it be an elif?
                if alloc_is_coverage:
                    alloc[prog.label] = prog.getBudget(coverage=spending)
                else:
                    alloc[prog.label] = spending

            # Store default budgets/coverages for all other programs if saturation is selected.
            elif 'saturate_with_default_budgets' in sim_settings and sim_settings['saturate_with_default_budgets'] is True:
                alloc[prog.label] = prog.getDefaultBudget(year=start_year) # Use default budget
            else:
                logger.warn("Program '%s' will not be used because no initial allocation was provided, 'saturate_with_default_budgets' not enabled." % prog.label)

        return alloc

    def update_cache(self,alloc,tvals,dt): # Do the stuff in precalculateprogsetvals
        # alloc must be a spending value - note that if the settings originally had alloc_is_coverage then the alloc would have been
        # passed though getBudget which means that it would go back through getCoverage at this step
        
        self.dt = dt
        self.tval_cache = tvals

        prog_vals = dict()

        for prog in self.progset.progs: # Iterate over classic Programs, not ModelPrograms
            if prog.label in alloc:
                prog.is_active = True
                prog_vals[prog.label] = {'cost':alloc[prog.label], 'cov':prog.getCoverage(budget=alloc[prog.label]), 'impact':{}}
                for par_label in prog.target_pars:
                    prog_vals[prog.label]['impact'][par_label] = prog.getImpact(prog_vals[prog.label]['cov'], impact_label=par_label, years=tvals, budget_is_coverage=True)
            else:
                prog.is_active = False
           
        # Finally, load the coverage and impact into the ModelPrograms
        # Also expand out any scalars
        for prog in self.programs:
            if prog.is_active:

                # Load cost
                prog.cost = prog_vals[prog.label]['cost']
                if not isinstance(prog.cost, np.ndarray) or prog.cost.size == 1:
                    prog.cost = np.full(tvals.shape, prog.cost)

                # Load coverage
                net_cov = prog_vals[prog.label]['cov']
                if not isinstance(net_cov, np.ndarray) or net_cov.size == 1:
                    net_cov = np.full(tvals.shape, net_cov)

                if prog.is_fraction:
                    prog.net_dt_cov = net_cov
                else:
                    prog.net_dt_cov = net_cov * self.dt

                # Load impacts
                prog.net_dt_impact = dict()
                for par in prog.pars:

                    impact = prog_vals[prog.label]['impact'][par.label]
                    # Todo - this might be possible to remove now that we use all years above
                    if not isinstance(impact, np.ndarray) or impact.size == 1:
                        impact = np.full(tvals.shape, impact)

                    if prog.is_fraction:
                        prog.net_dt_impact[par.uid] = impact
                    else:
                        prog.net_dt_impact[par.uid] = impact * self.dt

    def compute_pars(self,ti):
        # This function takes in a timestep and updates the parameter value in place
        contribs = defaultdict(list)

        # STAGE 3.1 Accumulate contributions
        for prog in self.programs:
            if prog.is_active:
                for par_id,contrib in prog.get_contribution(self.tval_cache[ti],ti).items():
                    contribs[par_id].append(contrib) # Append a (frac_cov,value) tuple

        # STAGE 3+4 Normalize and set output value
        impacts = dict()
        frac_dt_cov = dict()

        for par_id,c in contribs.items():
            frac_dt_cov[par_id] = [x[0] for x in c]
            total_cov = sum(frac_dt_cov[par_id] )
            impact = np.array([x[1] for x in c])

            if total_cov > 1:
                impact /= total_cov

            par = self.par_by_id[par_id]

            if par.links:
                impact = np.sum(impact) # Note - sum the rates before doing the conversion
                if par.units == 'fraction':
                    impact = 1.0 - (1.0 - impact) ** (1.0 / self.dt)
                elif par.units == 'number':
                    impact = impact * (1.0 / self.dt)
                elif par.units == 'proportion':
                    impact = impact
                else:
                    raise OptimaException('Unknown units')

            else:
                impact = np.sum(impact)

            if par_id in self.constraints:
                impact = np.clip(impact,self.constraints[par_id][0],self.constraints[par_id][1])

            impacts[par_id] = impact

        return impacts,frac_dt_cov


class ModelProgram(object):

    def __init__(self,label,is_fraction,pars,impact_groups,active_times,cost=None,net_dt_cov=None,net_dt_impact=None):

        self.uid = uuid.uuid4()
        self.label = label
        self.pars = pars # Parameters being reached
        self.is_fraction = is_fraction # True if the coverage and impact are in units of number, assert all Pars have same units
        self.impact_groups = impact_groups # For each impact group, list of parameter objects in that impact group

        self.cost = cost # Spending values at each time point - net_dt_cov and net_dt_impact is defined at same times
        self.net_dt_cov = net_dt_cov # Same for all Parameters, precomputed. If not a transition parameter, then the dt value is the same as the annual value!
        self.net_dt_impact = net_dt_impact # Dict mapping par_uid to impact array
        
        # Map parameter UID to impact group to identify which impact group to use for which each par reached
        # NOTE - this implies a parameter can only belong to one impact group!
        self.pars_to_groups = dict()
        for grp in self.impact_groups:
            for par in self.impact_groups[grp]:
                self.pars_to_groups[par.uid] = grp

        self.is_active = True # Only programs that are active will be used in ModelProgramSet.compute_pars()
        self.active_times = active_times

    def unlink(self):
        for i in xrange(0,len(self.pars)):
            self.pars[i] = self.pars[i].uid
        for grp in self.impact_groups:
            for i in xrange(0,len(self.impact_groups[grp])):
                self.impact_groups[grp][i] = self.impact_groups[grp][i].uid

    def relink(self,objs):
        for i in xrange(0,len(self.pars)):
            self.pars[i] = objs[self.pars[i]]
        for grp in self.impact_groups:
            for i in xrange(0,len(self.impact_groups[grp])):
                self.impact_groups[grp][i] = objs[self.impact_groups[grp][i]]

    def get_contribution(self,t,ti):
        # Return fractional coverage and program value contribution for every Parameter
        # reached by this Program
        # INPUTS
        # t - current simulation time (e.g. 2017.25), for activating/deactivating program
        # ti - current simulation time index (for looking up popsizes etc.)

        par_contribution = dict()

        # Return early if program is not active in this time
        if self.active_times is not None and (t < self.active_times[0] or t >= self.active_times[1]):
            return par_contribution

        # STAGE 2.1 Get impact group sizes
        source_set_size = defaultdict(float)

        for grp in self.impact_groups:
            if self.is_fraction:
                source_set_size[grp] = 1.0
            else:
                for par in self.impact_groups[grp]:
                    source_set_size[grp] += par.source_popsize(ti)

        # STAGE 2.2
        for par in self.pars:

            if not par.links: # If not a transition parameter, return the value directly
                par_contribution[par.uid] = (1.0,self.net_dt_impact[par.uid][ti])
            else:
                # Convert units depending on the *Program* units
                if self.is_fraction:
                    frac_dt_cov = self.net_dt_cov[ti]
                    frac_dt_impact = self.net_dt_impact[par.uid][ti]
                else:
                    grp_size = source_set_size[self.pars_to_groups[par.uid]]
                    if grp_size <= project_settings.TOLERANCE:
                        frac_dt_cov = 0.0
                        frac_dt_impact = 0.0
                    else:
                        frac_dt_cov =  self.net_dt_cov[ti] / grp_size # Note grp_size is 1.0 if self.is_fraction is True
                        frac_dt_impact = self.net_dt_impact[par.uid][ti] / grp_size

                # Convert units depending on the *Parameter* units
                if par.units == 'number':
                    eff_dt_impact = frac_dt_cov*par.source_popsize(ti)
                elif par.units == 'fraction' or par.units == 'proportion':
                    eff_dt_impact = frac_dt_impact
                else:
                    raise OptimaException('Unknown units!')

                par_contribution[par.uid] = (frac_dt_cov,eff_dt_impact) # Return whether or not the par is a fraction here to avoid having to look it up later

        return par_contribution

