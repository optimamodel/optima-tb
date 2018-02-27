from collections import defaultdict
import logging
logger = logging.getLogger(__name__)

class ModelProgramSet(object):
	# Contains some programs

	def __init__(self,progset,pops):
		# pops is model.pops

		self.uid = uuid.uuid4()
		self.programs = list()
		self.alloc = alloc

		self.progset = progset # A complete ProgramSet instance

		self.pars = list() # Array of Parameters that this ModelProgramSet will be responsible for updating
		self.par_id_to_par = list() # Map from UID to parameter object?

		# Construct the ModelProgram objects
		for prog in progset.progs # For each Program, instantiate a ModelProgram

			# Make a list of all impact Parameter objects
			impact_pars = []
			for pop in pops:
				impact_pars += [par for par in pop.pars if par.label in prog.target_pars]

			# Now populate the impact groups with those parameters
			impact_group_index = 0
			impact_groups = defaultdict(list)
			for par in impact_pars:
				if 'group' in settings.progtype_specs[prog.prog_type]['impact_pars'][par_label]:
					for group_label in progset.progtype_specs[prog.prog_type]['impact_par_groups']:
						impact_groups[group_label] += [par for par in impact_pars if par.label in progset.progtype_specs[prog.prog_type]['impact_par_groups'][group_label]]
				else:
					impact_groups[par_label] = [par] # The parameter becomes a member of a new impact group containing only itself

			is_fraction = prog.cov_format != 'number' # programs leave this as None otherwise (because it's a fraction and not probability?)

			self.programs.append(prog.label,is_fraction,impact_pars,impact_groups)

	def unlink(self):
	    for i in xrange(0,len(self.pars)):
	    	self.pars[i] = self.pars[i].uid
	    for prog in self.programs:
	    	prog.unlink()

	def relink(self,objs):
	    for i in xrange(0,len(self.pars)):
	    	self.pars[i] = objs[self.pars[i]]
    	for prog in self.programs:
    		prog.relink(objs)

	def update_cache(self,sim_settings): # Do the stuff in precalculateprogsetvals
		# sim_settings is the dict initialized in Model.build()
		
		# Take in a set of times and 
		start_year = sim_settings['progs_start']
		init_alloc = sim_settings['init_alloc']
		alloc_is_coverage = sim_settings['alloc_is_coverage']

        for prog in self.progset.progs: # Iterate over classic Programs, not ModelPrograms

            # Store budgets/coverages for programs that are initially allocated.
            if prog.label in init_alloc:
                alloc = init_alloc[prog.label]

                # If ramp constraints are active, stored cost and coverage needs to be a fully time-dependent array corresponding to timevec.
                if 'constraints' in sim_settings and 'max_yearly_change' in sim_settings['constraints'] and prog.label in sim_settings['constraints']['max_yearly_change']:
                    if alloc_is_coverage:
                        default = prog.getCoverage(budget=prog.getDefaultBudget(year=start_year))
                    else:
                        default = prog.getDefaultBudget(year=start_year)
                    default = prog.getDefaultBudget(year=start_year)
                    if np.abs(alloc - default) > project_settings.TOLERANCE:
                        alloc_def = sim_settings['tvec'] * 0.0 + default
                        alloc_new = sim_settings['tvec'] * 0.0 + alloc
                        try: eps = sim_settings['constraints']['max_yearly_change'][prog.label]['val']
                        except: raise OptimaException('ERROR: A maximum yearly change constraint was passed to the model for "%s" but had no value associated with it.' % prog.label)
                        if 'rel' in sim_settings['constraints']['max_yearly_change'][prog.label] and sim_settings['constraints']['max_yearly_change'][prog.label]['rel'] is True:
                            eps *= default
                        if np.isnan(eps): eps = np.inf  # Eps likely becomes a nan if it was infinity multiplied by zero.
                        if np.abs(eps * sim_settings['tvec_dt']) < np.abs(alloc - default):
                            if np.abs(eps) < project_settings.TOLERANCE:
                                raise OptimaException('ERROR: The change in budget for ramp-constrained "%s" is effectively zero. Model will not continue running; change in program funding would be negligible.' % prog.label)
                            alloc_ramp = default + (sim_settings['tvec'] - start_year) * eps * np.sign(alloc - default)
                            if alloc >= default: alloc_ramp = np.minimum(alloc_ramp, alloc_new)
                            else: alloc_ramp = np.maximum(alloc_ramp, alloc_new)
                            alloc = alloc_def * (sim_settings['tvec'] < start_year) + alloc_ramp * (sim_settings['tvec'] >= start_year)

                if alloc_is_coverage:
                    prog_vals[prog.label] = {'cost':prog.getBudget(coverage=alloc), 'cov':alloc, 'impact':{}}
                else:
                    prog_vals[prog.label] = {'cost':alloc, 'cov':prog.getCoverage(budget=alloc), 'impact':{}}

            # Store default budgets/coverages for all other programs if saturation is selected.
            elif 'saturate_with_default_budgets' in sim_settings and sim_settings['saturate_with_default_budgets'] is True:
                prog_vals[prog.label] = {'cost':prog.getDefaultBudget(year=start_year), 'cov':prog.getCoverage(budget=prog.getDefaultBudget(year=start_year)), 'impact':{}}

            else:
                logger.warn("Program '%s' was not contained in init_alloc and not saturated, therefore was not created." % prog.label)

            # Convert coverage into impact for programs.
            if prog.label in prog_vals:
                if 'cov' in prog_vals[prog.label]:
                    cov = prog_vals[prog.label]['cov']

                    # Check if program attributes have multiple distinct values across time.
                    do_full_tvec_check = {}
                    for att_label in prog.attributes:
                        att_vals = prog.attributes[att_label]
                        if len(set(att_vals[~np.isnan(att_vals)])) <= 1:
                            do_full_tvec_check[att_label] = False
                        else:
                            do_full_tvec_check[att_label] = True

                    # If attributes change over time and coverage values are greater than zero, impact functions based on them need to be interpolated across all time.
                    # Otherwise interpolating for the very last timestep alone should be sufficient.
                    # This means impact values can be stored as full arrays, single-element arrays or scalars in the case that an impact function has no attributes.
                    for par_label in prog.target_pars:
                        do_full_tvec = False
                        if 'attribs' in prog.target_pars[par_label] and np.sum(cov) > project_settings.TOLERANCE:
                            for att_label in prog.target_pars[par_label]['attribs']:
                                if att_label in do_full_tvec_check and do_full_tvec_check[att_label] is True:
                                    do_full_tvec = True
                        if do_full_tvec is True:
                            years = sim_settings['tvec']
                        else:
                            years = [sim_settings['tvec'][-1]]
                        prog_vals[prog.label]['impact'][par_label] = prog.getImpact(cov, impact_label=par_label, parser=parser, years=years, budget_is_coverage=True)

            # Finally, load the coverage and impact into the ModelPrograms
            for prog in self.programs:

            	# Load coverage
            	if prog.is_fraction:
            		prog.net_dt_cov = prog_vals[prog.label]['cov']
            	else:
            	    prog.net_dt_cov = prog_vals[prog.label]['cov'] * sim_settings['tvec_dt']

            	# Load impacts
            	prog.net_dt_impact = dict()
            	for par in pars:
            		if prog.is_fraction:
            			prog.net_dt_impact[par.uid] = prog_vals[prog.label]['impact'][par.label]
            		else:
            		    prog.net_dt_impact[par.uid] = prog_vals[prog.label]['impact'][par.label] * sim_settings['tvec_dt']

	def update_pars(self,ti):
		# This function takes in a timestep and updates the parameter value in place
		contribs = defaultdict(list)

		# STAGE 3.1 Accumulate contributions
		for prog in progs:
			for par_id,contrib in prog.get_contribution(ti).items():
				contribs[par_id].append(contrib) # Append a (frac_cov,value) tuple

		# STAGE 3+4 Normalize and set output value
		for par_id,c in contribs.items():
			net_cov = sum([x[0] for x in c])
			impacts = np.array([x[1] for x in c])

			if net_cov > 1
				impacts /= net_cov

			self.pars[self.par_id_to_par[par_id]].vals[ti] = sum(impacts)


class ModelProgram(Parameter)

	def __init__(self,label,is_fraction,pars,impact_groups,net_dt_cov=None,net_dt_impact=None):

		self.uid = uuid.UUID4()
		self.label = label
		self.pars = pars # Parameters being reached
		self.is_fraction = is_fraction # True if the coverage and impact are in units of number, assert all Pars have same units
		self.impact_groups = impact_groups # For each impact group, list of parameter objects in that impact group

		self.net_dt_cov = net_dt_cov # Same for all Parameters, precomputed. If not a transition parameter, then the dt value is the same as the annual value!
		self.net_dt_impact = net_dt_impact # Dict mapping par_uid to impact array
		
		# Map parameter UID to impact group to identify which impact group to use for which each par reached
		# NOTE - this implies a parameter can only belong to one impact group!
		self.pars_to_groups = dict()
		for grp in self.impact_groups:
			for par in self.impact_groups[grp]:
				self.pars_to_groups[par.uid] = grp 


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

	def get_contribution(self,ti):
		# Return fractional coverage and program value contribution for every Parameter
		# reached by this Program

		# STAGE 2.0 Get number of people reached by each parameter
		self.popsizes = {}
		for par in self.pars:
			source_element_size[par.uid] = par.source_popsize

		# STAGE 2.1 Get impact group sizes
		source_set_size = dict

		for grp in self.impact_groups:
			if self.is_fraction:
				source_set_size[grp] = 1.0
			else:
				source_set_size[grp] = sum(self.source_element_size[par.uid] for par in self.impact_groups[grp]) # Number of people covered by this impact group

		# STAGE 2.2
		par_contribution = dict
		for par in self.pars:

			if len(par.links) == 0 # If not a transition parameter, return the value directly
				par_contribution[par.uid] = (1.0,self.net_dt_impact[par.uid][ti])


			# Retrieve the fractional values for the parameter's impact group
			frac_dt_cov = grp_frac_dt_cov[par.uid]
			frac_dt_impact = grp_frac_dt_impact[par.uid]

			# Convert units
			# If the *Program* has number units - note grp_size is 1.0 if self.is_fraction
			grp_size = source_set_size[self.pars_to_groups[par.uid]]
			frac_dt_cov =  self.net_dt_cov[ti] / grp_size
			frac_dt_impact = self.net_dt_impact[par.uid][ti] / grp_size

			# If the *Parameter* has number units
			if par.is_fraction:
				eff_dt_cov = frac_dt_cov
				eff_dt_impact = frac_dt_impact
			else:
				eff_dt_cov = frac_dt_cov*source_element_size[par.uid]
				eff_dt_impact = frac_dt_cov*source_element_size[par.uid]

			par_contribution[par.uid] = (frac_dt_cov,eff_dt_cov*eff_dt_impact)

		return par_contribution

