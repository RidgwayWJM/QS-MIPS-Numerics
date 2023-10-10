//Driver for 2D chemically structured model with canonical QS kinetics
//Solves the coupled equations
// n_t = D(u)\nabla^2 n + epsilon*n_{uu} - (f(u,c)n)_u
// c_t = Dc\nabla^2c - beta*c + alpha_0*ubar
// ubar = int_0^\infty u*n du
//
 
#include "generic.h"	// Generic oomph-lib routines
#include "StratifiedHelmholtzElement2D.h" //custom helmholtz element
#include "DensityElement2D.h"	//custom density element
#include "meshes/simple_cubic_mesh.template.h"	//include the meshes
#include "meshes/simple_cubic_mesh.template.cc"
#include <cstdlib>	//needed for srand and rand
#include <ctime>	//needed for time

using namespace std;
using namespace oomph;

//==start_of_namespace================================================
/// Namespace for problem parameters
//====================================================================
namespace GlobalParams
{
	//=====================Basic Problem parameters==============================
	// domain parameters
	double Lu, Lx, Ly, domain_vol; // length of domain in u-, x-, y- directions, respectively (which correspond to x,y,z- directions in oomph-lib)

	// AI parameters
	double D2, beta_decay, alpha;

	//cell parameters
	double a, L, K, lambda, rho_star, epsilon_reg;

	//swimming parameters
	double D1_star, D1_p_star, D_inf_0;	//swimming motility coefficient, it's derivative at u=u*, and a parameter representing either D_inf or D_0
	
	//=======================Numerical parameters=================================================
	double epsilon_t, t_max, dt_out; //target error in adaptive time-stepping, simulation time, and minimum time between writing output files, respectively 

	//number of elements in each coordinate direction
	unsigned num_el_density_u, num_el_helm_u, num_el_x, num_el_y;
	//=======================documentation parameters=================================================
	//name of the restart file (if used)
	std::string Restart_filename = "";

	//=======================Derived Quantities in terms of Basic params===========================
	double u_star;

	//initial condition as a Vector (history values are provided using forward Euler)
	void get_IC(const double& time, const Vector<double>& x, Vector<double>& sol)
	{
		//sets initial condition as a spatially homogeneous uniform distribution
		//Note: x[0]=u, x[1]=x, and sol[0]=n, sol[1]=c, sol[2]=ubar, ,

		// ========== Gaussian distribution, spatially uniform =============
		double epsilon_IC = epsilon_reg; // width of the gaussian to use in the initial condition (cannot use with the unregularized problem, i.e epsilon_reg=0)
		//double epsilon_IC = 1e-4;
		double gaussian = std::exp(-0.5*lambda/epsilon_IC*std::pow(x[0]-u_star,2)); //gaussian factor that shows up in the steady-state and perturbation
		double normalization = std::sqrt(lambda/(2*M_PI*epsilon_IC));
		double ss = rho_star*normalization*gaussian;	//steady-state chemical structure
		// compute steady-state AI concentration	
		double ce_quad_B = K - rho_star*alpha*(a + L)/(beta_decay*lambda);
		double ce_quad_C = -a*alpha*rho_star*K/(beta_decay*lambda);
		double ce = 0.5*(-ce_quad_B + std::sqrt(ce_quad_B*ce_quad_B - 4.0*ce_quad_C)); //homogeneous equilibrium value of AI

		//random perturbation to steady-state
		double amplitude = 1e-2;
		double r = 2*(static_cast <double> (rand()) / static_cast <double> (RAND_MAX) - 0.5);	// random double between -1 and 1
		sol[0] = ss + amplitude*r*ss;
		r = 2*(static_cast <double> (rand()) / static_cast <double> (RAND_MAX) - 0.5);	// new random number
		sol[1] = ce + amplitude*r*ce;
		sol[2] = u_star*rho_star;
	}

	// function for evaluating the effective diffusion coefficient
	void get_D1(const Vector<double>& x, double* D1)	{
		// piecewise linear function
		double u1 = (D_inf_0 - D1_star)/D1_p_star + u_star;
		if(D1_p_star < 0)	{	//decreasing D1
			if(x[0] > u1)	{	
				*D1 = D_inf_0;	//constant
			}
			else	{
				*D1 = D1_p_star*(x[0] - u_star) + D1_star;	//linear
			}
		}
		else	{	//increasing D1
			if(x[0] > u1)	{	
				*D1 = D1_p_star*(x[0] - u_star) + D1_star;	//linear
			}
			else	{
				*D1 = D_inf_0; //constant
			}

		}
	}

} // end of namespace

using namespace GlobalParams;

// Start of Element class for enforcing non-local auxillary equation
// requires a pointer to the mesh containing all the density elements so that the non-local residuals can be calculated
// Row_ind_x,y are the x,y indices of the elements in the density mesh which influence the value of ubar in this element
// we imagine this ubar element occupying the same (x,y) space as the elements in the density mesh with x,y indices row_ind_x, row_ind_y
template<class DENSITY_ELEMENT>
class UbarElement : public virtual GeneralisedElement
{
	public:
	//constructor where data is created internally (within the UbarElement)
	UbarElement(SimpleCubicMesh<DENSITY_ELEMENT>* density_msh_pt, const unsigned row_ind_x, const unsigned row_ind_y) 
		: Density_msh_pt(density_msh_pt), //this line initializes the protected member variable Density_msh_pt equal to the argument density_msh_pt
		  Row_ind_x(row_ind_x),
		  Row_ind_y(row_ind_y)
	{
		Ndata_entries = unsigned(pow(3,Dim) + 0.5);
		//create new data object that will contain the unknowns on this element
		Ubar_data_pt = new Data(Ndata_entries);
		//declare this data as internal data for this element so that the current element is `in charge' of it
		Ubar_data_index = add_internal_data(Ubar_data_pt);
		Ubar_data_created_internally = true;
		
		//We now need to pass (pointers to) the external density data which are required in the residual and jacobian calculations
		//we add all data of all elements in the density mesh with the same fixed x- and y- coordinates as represented by Row_ind_x,y
		//number of elements in u and x direction
		unsigned n_element_u = Density_msh_pt->nx();	
		unsigned n_element_x = Density_msh_pt->ny();
		for(unsigned i=0; i<n_element_u; i++)	{
			//calculate index of element in the density mesh with same x,y coords as this ubar element
			unsigned ext_el_ind = n_element_u*n_element_x*Row_ind_y + n_element_u*Row_ind_x + i;  
			DENSITY_ELEMENT* ext_el_pt = dynamic_cast<DENSITY_ELEMENT*>(Density_msh_pt->element_pt(ext_el_ind));
			unsigned n_node = ext_el_pt->nnode();
			//loop over nodes in the external element and add external data for each
			for(unsigned j=0; j<n_node; j++)	{
				Density_external_data_ind.push_back(add_external_data(ext_el_pt->node_pt(j)));
			}
		}
	}

	//call default routine for eqn numbering as well as the local routine for storing eqn numbers of the ubar eqns
	void assign_all_generic_local_eqn_numbers(const bool& store_local_dof_pt)	{
		GeneralisedElement::assign_internal_and_external_local_eqn_numbers(store_local_dof_pt);
		this->assign_additional_local_eqn_numbers();
	}

	//store local eqn number of the equation for ubar 
	void assign_additional_local_eqn_numbers()	{
		if(Ubar_data_created_internally)  {	//ubar is an internal variable
			for(unsigned j=0; j<Ndata_entries; j++)  {
				Ubar_local_eqn.push_back(internal_local_eqn(Ubar_data_index,j));
			}
		}
		else	{	//ubar is an external variable
			for(unsigned j=0; j<Ndata_entries; j++)	{
				Ubar_local_eqn.push_back(external_local_eqn(Ubar_data_index,j));
			}
		}
	}

	//elemental residual calculation
	void fill_in_contribution_to_residuals(Vector<double>& residuals)	{
		Vector<double> first_moments(Ndata_entries);
		for(unsigned j=0; j<Ndata_entries; j++)	{
			first_moments[j]=0.0;
		}
	        compute_first_moment(first_moments);	//compute the non-local term needed for this (Ubar) element
		for(unsigned j=0; j<Ndata_entries; j++)	{
			if(Ubar_local_eqn[j] >= 0)	{	// should always be true since no BCs for ubar
				residuals[Ubar_local_eqn[j]] = Ubar_data_pt->value(j) - first_moments[j];
			}
		}
	}

	unsigned ndof_types() const	{
		return 1;
	}

	//creates a list of pairs of all unknowns in this element. First entry contains the global equation 
	//number and the second contains the number of the 'dof_type' that this unknown is associated with
	void get_dof_numbers_for_unknowns(
		std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
	{
		if(Ubar_data_created_internally)	{
			// temporary pair (used to store dof lookup prior to being
			// added to list)
			std::pair<unsigned long, unsigned> dof_lookup;
			for(unsigned j=0; j<Ndata_entries; j++)	{
				// determine local eqn number for ubar eqn
				int local_eqn_number = Ubar_local_eqn[j];
				// Is it a dof or is it pinned?
				if (local_eqn_number >= 0)	{
					// store dof lookup in temporary pair: First entry in pair
					// is global equation number; second entry is dof type
					dof_lookup.first = this->eqn_number(local_eqn_number);
					dof_lookup.second = 0;
					// Add to list
					dof_lookup_list.push_front(dof_lookup);
				}
			}
		}
	}

	//computes the non-local elemental contribution to the ubar residuals
	//In this routine we need access to data in the density mesh which we can look-up using knowledge element ordering within oomph-lib (see specific mesh documentation)
	//The density mesh exists in (u,x,y) space (IN THIS ORDER), indexing is: u first, then x, then y
	void compute_first_moment(Vector<double>& first_moments)	{
		// number of elements in u,x  direction for density mesh
		unsigned num_el_density_u = Density_msh_pt->nx();
		unsigned num_el_density_x = Density_msh_pt->ny();
		
		double du = 0.5*Lu/num_el_density_u; // distance between nodes in each element (rectangular quadratic elements only)
		// storage for integrals (one for each node at fixed x in the element)
		Vector<double> s(Dim+1); // storage for local coordinates in each element
		
		//go through all elements in the density mesh with fixed x,y (i.e. a column of elements in u,x,y-space). Note there there are 
		//nodes within each element that have different values of x,y
		for(unsigned j=0; j<num_el_density_u; j++)	{
			//index of the element in the density mesh that influences ubar values in this element
			unsigned elem_ind = num_el_density_u*num_el_density_x*Row_ind_y + Row_ind_x*num_el_density_u + j;	//goes through all elements with fixed x,y in the density mesh
			// upcast from generalisedElement to present type
			DENSITY_ELEMENT* el_pt = dynamic_cast<DENSITY_ELEMENT*>(Density_msh_pt->element_pt(elem_ind));
			unsigned num_node = el_pt->nnode(); 
			// loop through nodes in this element and add the contribution to the appropriate ubar value
			for(unsigned k=0; k<num_node; k++)	{	 
				double nval = el_pt->field_value_density(k); // value of n stored at this node
				el_pt->local_coordinate_of_node(k, s); // local coordinates of kth node
				//global coordinates of node
				double uval = el_pt->node_pt(k)->x(0); 
				// Add the contribution to the appropriate integral at fixed x,y (s[1])
				int integral_ind = 3*int(s[2] + 1.5) + int(s[1] + 1.5); // add an additional 0.1 to avoid truncation errors (range of local coords is [-1,1])
				if(integral_ind > (Dim+1)*(Dim+1)-1 || integral_ind < 0)	{
					throw 
						OomphLibError("Trouble! Evaluation of non-local integral failed. Integral index out of bounds..\n",
			                     	   "FullDensityProblem::update_moment_nodal_values()",
			                     	   OOMPH_EXCEPTION_LOCATION);
				}

				if(std::abs(s[0]) < 1e-2)  	{//middle column where s[0]=0 
					first_moments[integral_ind] += 4.0/3.0*uval*nval*du; // contribution from mid-point of element (in u) according to Simpson's rule
				}
				else	{
					first_moments[integral_ind] += uval*nval*du/3.0; // contribution from end-points of element (in u) according to Simpson's rule
				}
			}
		}
	}// end of compute_first moment
	
	//pointer to data object that contains the three values of ubar in this element
	Data* ubar_data_pt() const	{
		return Ubar_data_pt;
	}

	//access functions to row indices
	unsigned x_ind()	{
		return Row_ind_x;	
	}
	unsigned y_ind()	{
		return Row_ind_y;	
	}

protected:
	
	bool Ubar_data_created_internally;
	//Mesh for density elements
	SimpleCubicMesh<DENSITY_ELEMENT>* Density_msh_pt;

	//Pointer to data item containing the value of ubar in this element
	Data* Ubar_data_pt;

	//index required to obtain Ubar_data_pt from internal_data_pt(...)
	unsigned Ubar_data_index;

	//number of entries in the data for this element
	unsigned Ndata_entries;
	
	//row indices corresponding to the value of x,y with which this element is associated, x = Row_ind_x*dx, y = Row_ind_y*dy
	const unsigned Row_ind_x, Row_ind_y;

	//spatial dimension
	unsigned Dim=2;

	// Local equation number of the Ubar equation
	Vector<int> Ubar_local_eqn;

	Vector<unsigned> Density_external_data_ind;
}; // end of UbarElement class

// =========================  Start of Multi-domain Upgrade of the Density Element ========================
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
class DensityHelmholtzElement2D
	: public virtual DENSITY_ELEMENT,
	public virtual ElementWithExternalElement
{
public:
	//Constructor: call underlying constructors
	DensityHelmholtzElement2D()
		:DENSITY_ELEMENT(), ElementWithExternalElement()
	{
		// one interaction with other elements (in the kinetic term)
		this->set_ninteraction(1);
	}

	//overload the get_kinetics... function from the density element. This provides the coupling between the density and
	// helmholtz equations
	void get_kinetics_density(const unsigned& ipt, const Vector<double>& x, double* kinetics) const		{
		//interaction index
		unsigned interaction = 0;

		//get the AI concentration interpolated from the external element
		const double interpolated_c = dynamic_cast<HELM_ELEMENT*>(
				external_element_pt(interaction,ipt))->interpolated_c_strat_helmholtz(
				external_element_local_coord(interaction,ipt));
		// calculate the kinetics
		*kinetics = a + L*interpolated_c/(K + interpolated_c) - lambda*x[0]; //canonical QS kinetics
	}

	// compute the elements jacobian matrix with finite differencing
	void fill_in_contribution_to_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian)	{
		//fill in by finite differencing
		ElementWithExternalElement::fill_in_contribution_to_jacobian(residuals, jacobian);
	}


};//end of DensityHelmholtzElement

 // =========================  Start of Multi-Domain upgrade of helmholtz element Element ========================
template<class HELM_ELEMENT, class DENSITY_ELEMENT>	//technically the density element is not needed here since the interaction from ubar comes from external data
class HelmholtzDensityElement2D
	: public virtual HELM_ELEMENT,
	public virtual ElementWithExternalElement
{
public:
	//Constructor: call underlying constructors
	HelmholtzDensityElement2D()
		:HELM_ELEMENT(), ElementWithExternalElement()
	{
		this->set_ninteraction(1);
	}
	
	//overload the get_source... function from the stratified helmholtz element which provides the coupling
	// between the helmholtz and moment equations
	void get_source_strat_helmholtz(const unsigned& ipt, const Vector<double>& x, double* source) const
	{
		//spatial component of local coords
		Vector<double> xloc(Dim,0.0);
		for(unsigned k=0; k<Dim; k++)	{
			xloc[k] = integral_pt()->knot(ipt,k+1); // arguments are (i,j) where i is the index of the integration poitn and j is the jth coordinate 
			//note that knot(ipt,0) is the u-coord which we dont need here
		}
		
		// interpolate the external data (i.e. the 3^Dim values of ubar) quadratically
		Vector<double> basis_x(3);
		Vector<double> basis_y(3);
		//calculate basis functions decomposed as products of 1D shape functions
		oneDshape_strat_helmholtz(xloc[0], basis_x);
		oneDshape_strat_helmholtz(xloc[1], basis_y);
		double ubar = 0;
		for(unsigned i=0; i<3; i++)	{ //3 nodes in each direction
			for(unsigned j=0; j<3; j++)	{
				//the indexing of basis functions has to be the same as the ubar data values
				double ubar_dat = external_data_pt(0)->value(3*i + j); 
				ubar+= basis_x[j]*basis_y[i]*ubar_dat;			}
		}
		*source = alpha*ubar;
	}

	// Calculate the element's contribution to the residual vector. We do this by
	// overloading the fill_in_... function for the residuals
	void fill_in_contribution_to_residuals(Vector<double>& residuals)	{
		// fill in the residuals of the helmholtz equation
		StratifiedHelmholtzElement2D::fill_in_contribution_to_residuals(residuals);
	}

	// Calculate the Jacobian by finite differencing, we overload the fill_in_... function
	void fill_in_contribution_to_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian){
		ElementWithExternalElement::fill_in_contribution_to_jacobian(residuals, jacobian); //finite diffeerences for all elements
	}
	//evaluates one-dimensional shape functions, used for interpolation
	void oneDshape_strat_helmholtz(const double& s, Vector<double>& shape)	const {
		shape[0] = s*(s - 1.0)/2.0;	//zero at xloc = 0,1, one at xloc=-1
		shape[1] = (1.0-s)*(1.0 + s);	//zero at xloc = -1,1, one at xloc=0
		shape[2] = s*(s + 1.0)/2.0;	//zero at xloc = -1,0, one at xloc=1

	}
protected:
	//spatial dimension of the problem
	unsigned Dim = 2;

};//end of HelmholtzDensityElement

//=============================== Start of problem class for combined problem===================
// Combined Density, bulk AI problem
// /============================================================================================
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
class FullDensityProblem : public Problem
{
public:
	// Constructor
	FullDensityProblem(const unsigned& num_el_u_fine, const unsigned& num_el_u_coarse,
		       	const unsigned& num_el_x, const unsigned& num_el_y); 

	//Destructor 
	~FullDensityProblem();

	//actions taken before Newton iteration at each timestep (empty)
	void actions_before_newton_solve() {}
	//actions taken after the Newton iteration at each timestep
	void actions_after_newton_solve() {}

	//actions taken before each newton step within one newton solve
	void actions_before_newton_step() {}

	//actions taken after each newton step within one newton solve	
	void actions_after_newton_step() {}

	//actions taken before the residuals are checked in the newton iteration
	void actions_before_newton_convergence_check() {}

	//actions taken before distribution (but after partitioning of the mesh)
	void actions_before_distribute()	{}
	//actions taken after the problem has been distributed over multiple processors
	void actions_after_distribute()	{
		oomph_info << "Problem succesfully distributed. In actions_after_distribute(...)\n";
		Multi_domain_functions::setup_multi_domain_interactions
			<DENSITY_ELEMENT, HELM_ELEMENT>(this, density_msh_pt(), helm_msh_pt());
		oomph_info <<"Succesffully re-setup multi-domain interactions\n";
	}

	//create the elements and mesh responsible for the nonlocal ubar auxillary equation
	void create_ubar_elements();

	//assigns initial conditions 
	void set_initial_condition();

	// Global error norm for adative timestepping: RMS error, based on difference betweeen 
	// predicted and actual value at all nodes
	double global_temporal_error_norm();

	//computes the total population of the density solution
	double compute_solution_mass(int itime);	
	// document the solution
	void doc_solution();
	
	//write restart file
	void doc_restart();

	//dump problem to disk to allow for restart
	void dump_it(ofstream& dump_file);

	//sets up restart from restart file specified via command line
	void restart();

	// read solution from disk for restart
	void restart(ifstream& restart_file);

	// function that outputs important problem parameters
	void output_parameters_to_file(ofstream& param_file);
	
	//access function for mesh used by the density equation
	SimpleCubicMesh<DENSITY_ELEMENT>* density_msh_pt()	{
		return dynamic_cast<SimpleCubicMesh<DENSITY_ELEMENT>*>(Density_msh_pt);
	}

	//access function for mesh used by the helmholtz equaiton
	SimpleCubicMesh<HELM_ELEMENT>* helm_msh_pt()	{
		return dynamic_cast<SimpleCubicMesh<HELM_ELEMENT>*>(Helm_msh_pt);
	}

	//access function to Next_dt
	double& next_dt() { return Next_dt;}
	unsigned doc_number() { return unsigned(Doc_info.number()); }
protected:
	//Mesh for density elements
	SimpleCubicMesh<DENSITY_ELEMENT>* Density_msh_pt;
	//mesh for stratified helmholtz elements
	SimpleCubicMesh<HELM_ELEMENT>* Helm_msh_pt;
	//mesh for the set of ubar elements
	Mesh* Ubar_msh_pt;
	//spatial dimension of the problem
	unsigned Dim = 2;
	//indices of elements in the ubar mesh in the order in which they are added
	
private:
	//DocInfo objects for outputing solution and restart files to disk 
	DocInfo Doc_info;
	DocInfo Restart_info;
	//suggested next timestep which is stored here so that it can be read/written to/from restart files
	double Next_dt;
}; //end of problem class

///================ Start of  Problem constructor==============================================
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::FullDensityProblem(const unsigned& num_el_u_fine, const unsigned& num_el_u_coarse, 
		const unsigned& num_el_x, const unsigned& num_el_y)
{
	//setup labels for ouptut
	Doc_info.set_directory(CommandLineArgs::Argv[2]);
	Restart_info.set_directory(CommandLineArgs::Argv[2]);	//use same directory for restart files as for output files

	Doc_info.number() = 0;
	Restart_info.number() = 0;

	//area of the (spatial) domain
	domain_vol = Lx*Ly;

	// allocate timestepper with adaptive timestepping (true) or not (false)
	add_time_stepper_pt(new BDF<2>(true));
	
	//specify that the problem is nonlinear and requires more than one newton iteration at each timestep	
	Problem_is_nonlinear=true;

	//specify that we don't need to keep error below tolerance as long as the decrease in the timestep is small
	Keep_temporal_error_below_tolerance=false;
	DTSF_min_decrease=0.9;		//minimum allowed relative decrease in the time step

	Max_newton_iterations = 8;	// maximum allowable iterations for non-linear solver before giving up	

	//Create the mesh objects for each equation
	//parameter to the mesh. The argument to the constructor indicates
	//the number of elements in the mesh.
	double u_min_mesh = max(0.0, u_star - Lu/2.0); //centre the mesh around u_star which is (roughly) where we expect the solution to be localized
	Density_msh_pt = new 
		SimpleCubicMesh<DENSITY_ELEMENT>(num_el_u_fine, num_el_x, num_el_y, 
				u_min_mesh, u_min_mesh + Lu, 0.0, Lx, 0.0, Ly, time_stepper_pt());
	Helm_msh_pt = new
		SimpleCubicMesh<HELM_ELEMENT>(num_el_u_coarse, num_el_x, num_el_y,
				u_min_mesh, u_min_mesh + Lu, 0.0, Lx, 0.0, Ly, time_stepper_pt());
	//store copies of all elements on each processor as `halo elements'
	Density_msh_pt->set_keep_all_elements_as_halos();
	Helm_msh_pt->set_keep_all_elements_as_halos();
	
	oomph_info << "Created sub meshes for density and helmholtz-moment elements.\n";
	        
	// Loop over the elements in the density mesh to setup element specific things that cannot be done by element constructor
	unsigned n_element_density = density_msh_pt()->nelement();	
	for(unsigned i=0; i<n_element_density;i++)	{
		//upcast from generalised element to present type
		DENSITY_ELEMENT *el_pt = dynamic_cast<DENSITY_ELEMENT*>(density_msh_pt()->element_pt(i));
	
		//set the diffusion-in-u coefficient and motility coefficient in the density eqution
		el_pt->epsilon_diffusion_pt() = &epsilon_reg;
		el_pt->motility_fct_pt() = get_D1;

		//ignore external geometric data in 'external' helmholtz element when computing jacobian matrix
		//because the interaction does not involve derivatives
		el_pt->ignore_external_geometric_data();
	}

	//create sub mesh of ubar elements
	create_ubar_elements();
	oomph_info <<"Created submesh for auxillary field ubar(x).\n";

	// Loop over the elements in the helmholtz_moment mesh to setup element specific things that cannot be done by mesh contrustor
	unsigned n_element_helm = helm_msh_pt()->nelement();
	for(unsigned i=0; i<n_element_helm; i++)	{
		//upcast from generalised element to present type
		HELM_ELEMENT *el_pt = dynamic_cast<HELM_ELEMENT*>(helm_msh_pt()->element_pt(i));
		
		//set the decay and diffusion coefficient pointers in the helmholtz equation
		el_pt->decay_coeff_pt() = &beta_decay;
		el_pt->diff_coeff_pt() = &D2;
		unsigned ubar_el_ind = i/num_el_u_coarse; //find the index of the ubar element that is needed for the ith helmholtz element
		//upcast this element from GeneralisedElement to present type	
		UbarElement<DENSITY_ELEMENT>* ubar_el_pt 
			= dynamic_cast<UbarElement<DENSITY_ELEMENT>*>(Ubar_msh_pt->element_pt(ubar_el_ind));
		//add the ubar data in this ubar element to the list of external data in the current helmholtz element (this
		//tells the element to include it in such things as jacobian calculation)
		el_pt->add_external_data(ubar_el_pt->ubar_data_pt());
	}

	oomph_info << "Element-specific quantities initialized.\n";

	// combine the submeshes and build the combined mesh
	add_sub_mesh(Density_msh_pt);
	add_sub_mesh(Helm_msh_pt);
	add_sub_mesh(Ubar_msh_pt);
	build_global_mesh();	
	oomph_info << "Global Mesh built.\n";

	//setup the interactions between the density and helmholtz submeshes
	Multi_domain_functions::setup_multi_domain_interactions<DENSITY_ELEMENT,HELM_ELEMENT>(this, density_msh_pt(), helm_msh_pt());	

	oomph_info << "Interactions between submeshes setup\n";
	//Assign the global and local equations numbers for the problem
	oomph_info << "Number of equations is " << assign_eqn_numbers(false) << std::endl;	//need to set to false for MPI problems	

}  // end of constructor

//=================Start of destructor=====================================
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::~FullDensityProblem()	{
	//cleanup
	delete Density_msh_pt;
	delete Helm_msh_pt;
	delete Ubar_msh_pt;
}

//==================Start of create_ubar_elements()==========================
// creates the elements and meshes required for the nonlocal ubar auxillary equation
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::create_ubar_elements()	{
	Ubar_msh_pt = new Mesh;
	for(unsigned j=0; j<num_el_y; j++)	{
		for(unsigned i=0; i<num_el_x; i++)	{
			//create the new element
			UbarElement<DENSITY_ELEMENT>* ubar_el_pt =
				new UbarElement<DENSITY_ELEMENT>(density_msh_pt(), i, j);
			Ubar_msh_pt->add_element_pt(ubar_el_pt);	//add to mesh
		}
	}
	Ubar_msh_pt->set_keep_all_elements_as_halos();
}

//==================Start of compute_solution_mass()==========================
// computes the total population in the solution, not used in the actual numerical solution of the problem, 
// rather to ensure that the IC has the correct population. Also just a useful helper function for processing.
// Argument is the time at which the solution is evaluated. itime =0 is the current solution
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
double FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::compute_solution_mass(int itime)	{
	double mass = 0.0;
	unsigned num_el = density_msh_pt()->nelement();
	for(unsigned i=0; i<num_el; i++)	{
		//pointer to ith element
		DENSITY_ELEMENT* el_pt = dynamic_cast<DENSITY_ELEMENT*>(density_msh_pt()->element_pt(i));
		// number of nodes in the element
		const unsigned n_node = el_pt->nnode();
		//allocate memeory for shape functions
		Shape psi(n_node);
		DShape dpsi(n_node,Dim + 1); // dummy storage for calling dshape_eulerian

		// memory for local coords
		Vector<double> s(Dim + 1);
				
		const unsigned n_intpt = el_pt->integral_pt()->nweight();
	
		// loop over integration points
		for(unsigned ipt=0; ipt<n_intpt; ipt++)	{
			//local coordinates of integration point
			for(unsigned k=0; k<Dim + 1; k++)	{
				s[k] = el_pt->integral_pt()->knot(ipt,k); // arguments are (i,j) where i is the index of the integration poitn and j is the jth coordinate
			}
				
			//get the integral weight
			double w = el_pt->integral_pt()->weight(ipt);
			//call the derivatives of the shape and test
			double J = el_pt->dshape_eulerian(s, psi, dpsi);
			double W = w*J;
			// Calculate the interpolated density 
			double interpolated_n = 0.0;
			for(unsigned j=0; j<n_node; j++)	{
				//interpolated_n += el_pt->field_value_density(j)*psi[j];
				interpolated_n += el_pt->node_pt(j)->value(itime,0)*psi[j];
			}
			// add contribution to total mass 
			mass += interpolated_n*W;
		}//end of loop over integration points
	} //end of loop over elements
	return mass;
}//end of compute_solution_mass(...)

//==================Start of restart()==========================
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::restart()	{
	// Pointer to restart file
	ifstream* restart_file_pt=0;
	//get restart filename
	char filename[100];
	sprintf(filename,"%s_on_proc%i.dat",GlobalParams::Restart_filename.c_str(),
			this->communicator_pt()->my_rank()); 
	// Open restart file
	restart_file_pt= new ifstream(filename,ios_base::in);
	if (restart_file_pt!=0)	{
		oomph_info << "Opened " << filename << 
			" for restart. " << std::endl;
	}
	else	{
   		std::ostringstream error_stream;
   		error_stream 
	  		<< "ERROR while trying to open " << filename << 
		    	" for restart." << std::endl;
	     
   		throw OomphLibError(
		  		error_stream.str(),
		  		OOMPH_CURRENT_FUNCTION,
		  		OOMPH_EXCEPTION_LOCATION);
	}
	// Read restart data:
	if (restart_file_pt!=0)	{
		// Read the data from restart file
     		restart(*restart_file_pt);
	}
}

// read solution from disk for restart
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::restart(ifstream& restart_file)	{
	
	//get the predicted next dt and solution indices for documentation from the restart file
	double local_next_dt=0.0;
	unsigned local_doc_info_number=0, local_restart_info_number=0;

	if(restart_file.is_open())	{
		oomph_info << "restart file exists" << endl;

		//read in the predicted next value of dt
		string input_string;
		getline(restart_file, input_string, '#');
		restart_file.ignore(80,'\n');	//ignore rest of line
		local_next_dt = double(atof(input_string.c_str()));

		//Read in the index of the last output file from the restart file and have the soln output files indexed from here
		getline(restart_file, input_string, '#');	//get up to terminating comment
		restart_file.ignore(80, '\n');			//ignore rest of line
		local_doc_info_number = unsigned(atof(input_string.c_str()));	//read in step number

		//Read in the index of the restart file and index future outputs from there
		getline(restart_file, input_string, '#');	//get up to terminating comment
		restart_file.ignore(80, '\n');			//ignore rest of line
		local_restart_info_number = unsigned(atof(input_string.c_str()));	//read in step number

	}
	else	{
		oomph_info << "restart file does not exist" << endl;
	}
	//determine next dt from all local (processor) values
	double next_dt = 0.0;
	MPI_Allreduce(&local_next_dt, &next_dt, 1, MPI_DOUBLE, MPI_MAX, 
			communicator_pt()->mpi_comm());
	Next_dt = next_dt;
	oomph_info << "Next dt: " << Next_dt << endl;

	//now determine the doc_info number for docementing solution at future times
	unsigned doc_info_number = 0;
	MPI_Allreduce(&local_doc_info_number, &doc_info_number, 1, MPI_UNSIGNED, MPI_MAX,
			communicator_pt()->mpi_comm());
	Doc_info.number() = doc_info_number;
	oomph_info << "Restarted Doc_info.number(): " << doc_info_number << endl;
	
	//now determine the doc_info number of the restart file
	unsigned restart_info_number = 0;
	MPI_Allreduce(&local_restart_info_number, &restart_info_number, 1, MPI_UNSIGNED, MPI_MAX,
			communicator_pt()->mpi_comm());
	Restart_info.number() = restart_info_number + 1;	//needs to be incremented so that the next restart file doesn't overwrite the current one
	oomph_info << "Restarted Restart_info.number(): " << restart_info_number << endl;

	//read in data from restart file
	Problem::read(restart_file);
}

//==================Start of set_initial_conditions()==========================
// sets the initila condtion
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::set_initial_condition()	{
	 // Assign initial condition from specified IC
	 //---------------------------------------------
	#ifdef OOMPH_HAS_MPI
		unsigned nproc = MPI_Helpers::communicator_pt()->nproc();	//number of processors
	#endif

	 // Choose initial timestep	
     	double dt0=1e-3;

	// Initialise timestep -- also sets the weights for all timesteppers
	// in the problem.
	initialise_dt(dt0);
	
	//use this as the first suggested timestep
	Next_dt = dt0;

	// Backup time in global Time object
	double backed_up_time=time_pt()->time();
		   
	// Past history fo needs to be established for t=time0-deltat, ...
	// Then provide current values (at t=time0) which will also form
	// the initial guess for the first solve at t=time0+deltat
		   
	// Vector of solution values
	Vector<double> soln_local(3);	//soln in current processor
	Vector<double> soln(3);	//soln synced across processors
	Vector<double> x(Dim+1);
	//element spacing
	double dx = Lx/num_el_x, dy = Ly/num_el_y;

	//Find number of nodes in density mesh
	unsigned num_nod_density = density_msh_pt()->nnode();
	unsigned num_nod_helm = helm_msh_pt()->nnode();

	// Set continuous times at previous timesteps
	int nprev_steps=time_stepper_pt()->nprev_values();
	Vector<double> prev_time(nprev_steps+1);
	for (int itime=nprev_steps;itime>=0;itime--)	{
   		prev_time[itime]=time_pt()->time(unsigned(itime));
	} 
	// Loop over current & previous timesteps
	for (int itime=nprev_steps;itime>=0;itime--)	{
		// Continuous time
   		double time=prev_time[itime];
   		oomph_info << "setting IC at time =" << time << std::endl;
		     
   		// Loop over the nodes in density mesh to set initial guess everywhere
		for (unsigned i=0;i<num_nod_density;i++)     {
  			// Get nodal coordinates
			for(unsigned k=0; k<Dim+1; k++)	{
  				x[k]=density_msh_pt()->node_pt(i)->x(k);
  			}
			// Get intial solution
 			GlobalParams::get_IC(time,x,soln_local); 
			//sync IC across processors (required due to RNG in perturbation)
			for(unsigned ii=0; ii<3; ii++)	{
				#ifdef OOMPH_HAS_MPI
					double *soln_local_ii_vec = (double*)malloc(sizeof(double) * nproc);
					MPI_Allgather(&soln_local[ii], 1, MPI_DOUBLE, 
							soln_local_ii_vec, 1, MPI_DOUBLE, 
							communicator_pt()->mpi_comm());
					soln[ii] = soln_local_ii_vec[0];	//just pick the IC on one of the procs
				#else
					soln[ii] = soln_local[ii];
				#endif
			}	

  			// Assign initial condition for n
  			density_msh_pt()->node_pt(i)->set_value(itime,0,soln[0]);	//n is at index 0 at each node

  			// Loop over coordinate directions: Previous position = present position
  			for (unsigned j=0;j<Dim+1;j++)	{
  				density_msh_pt()->node_pt(i)->x(itime,j)=x[j];
  			}
   		}
		//correct the IC so that the total population is equal to the value specified (value initially differs to to RNG
		//in IC)
       		double total_population = compute_solution_mass(itime);
   		oomph_info << "total population is " << total_population <<"...adding a constant to the IC to correct... \n";
   		for(unsigned i=0; i<num_nod_density; i++)	{
   			double nval = density_msh_pt()->node_pt(i)->value(itime,0);
   			density_msh_pt()->node_pt(i)->set_value(itime,0,nval*rho_star*domain_vol/total_population);
   		}
		oomph_info << "total population is now " << compute_solution_mass(itime) << endl;

   		// Loop over the nodes in helmholtz mesh to set initial guess everywhere
		for(unsigned i=0; i<num_nod_helm; i++)	{
		   	//calculate nodal coords
		     	for(unsigned k=0; k<Dim+1; k++)	{
   				x[k]=helm_msh_pt()->node_pt(i)->x(k);
   			}
  			// Get initial condition at current nodal coords
   			GlobalParams::get_IC(time,x,soln_local);
			//sync IC across processors (required due to RNG in perturbation)
			for(unsigned ii=0; ii<3; ii++)	{
				#ifdef OOMPH_HAS_MPI
					double *soln_local_ii_vec = (double*)malloc(sizeof(double) * nproc);
					MPI_Allgather(&soln_local[ii], 1, MPI_DOUBLE,
							soln_local_ii_vec, 1, MPI_DOUBLE, 
							communicator_pt()->mpi_comm());
					soln[ii] = soln_local_ii_vec[0];	//just pick the IC on one of the procs
				#else
					soln[ii] = soln_local[ii];
				#endif
			}

   			//assign this initial condition to all nodes in the mesh with the same values of x and y
   			//so that c is independent of u
   			double x_i = x[1], y_i = x[2];
   			for(unsigned i2=0; i2<num_nod_helm; i2++)	{
   				double xcurr = helm_msh_pt()->node_pt(i2)->x(1);
   				double ycurr = helm_msh_pt()->node_pt(i2)->x(2);
   				if(abs(xcurr - x_i) < 0.05*dx && abs(ycurr - y_i) < 0.05*dy)	{
   					helm_msh_pt()->node_pt(i2)->set_value(itime, 0, soln[1]);
   				}
   			}
   			//loop over coordinate directions and set previous nodal position = present nodal position
  			for (unsigned k=0;k<Dim+1;k++)	{
   				helm_msh_pt()->node_pt(i)->x(itime,k)=x[k];
   			}
   		}
   		//assign values to ubar 
   		unsigned n_ubardata = unsigned(pow(3,Dim) + 0.5);
   		for(unsigned j=0; j<num_el_y; j++)	{
   			for(unsigned i=0; i<num_el_x; i++)	{
   				for(unsigned k=0; k<n_ubardata; k++)	{// loop over the values stored in each ubar element
   					//calculate the coordinates with which the data values in the ubar element are associated
   					unsigned data_ind_x = k%3;
   					unsigned data_ind_y = k/3;
   					x[1] = i*dx + data_ind_x*dx/2.0;	//value of x[0] doesn't matter here
   					x[2] = j*dy + data_ind_y*dy/2.0;
   					GlobalParams::get_IC(time,x,soln_local);
					//sync IC across processors (required due to RNG in perturbation)
					for(unsigned ii=0; ii<3; ii++)	{
						#ifdef OOMPH_HAS_MPI
							double *soln_local_ii_vec = (double*)malloc(sizeof(double) * nproc);
							MPI_Allgather(&soln_local[ii], 1, MPI_DOUBLE,
								soln_local_ii_vec, 1, MPI_DOUBLE, 
								communicator_pt()->mpi_comm());
							soln[ii] = soln_local_ii_vec[0];	//just pick the IC on one of the procs
						#else
							soln[ii] = soln_local[ii];
						#endif
					}
   					//get index of the ubar element associated with these coordinates
   					unsigned ubar_el_ind = num_el_x*j + i;
   					//upcast from generalised element to present type
   					UbarElement<DENSITY_ELEMENT>* ubar_el_pt = 
   						dynamic_cast<UbarElement<DENSITY_ELEMENT>*>(Ubar_msh_pt->element_pt(ubar_el_ind));
   					ubar_el_pt->ubar_data_pt()->set_value(k,soln[2]);
   				}
   			}
   		}
	} //end of loop over previous times   
	// Reset backed up time for global timestepper
	time_pt()->time()=backed_up_time;
}// end of set_initial_condition()	

// Global error norm for adative timestepping: RMS error, based on difference betweeen 
// predicted and actual value at all nodes
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
double FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::global_temporal_error_norm()
{
	double global_error = 0.0;	   
	 //number of nodes in the density and helmholtz meshes
	 unsigned n_node_density = density_msh_pt()->nnode();
	 unsigned n_node_helm = helm_msh_pt()->nnode();
		
	 //Loop over the nodes in the density mesh and calculate the estimated error in the values
	 for(unsigned i=0;i<n_node_density;i++)
	  {
	   // Get error in solution: Difference between predicted and actual
	   // value for each solution stored at ith node
	   double error = density_msh_pt()->node_pt(i)->time_stepper_pt()->
		   temporal_error_in_value(density_msh_pt()->node_pt(i),0);
	   //Add the square of the individual error to the global error	
	   global_error += error*error;
	  }
	//Loop over the nodes in the helmholtz mesh and calculate the estimated error in the values
	 for(unsigned i=0;i<n_node_helm;i++)
	  {
	   // Get error in solution: Difference between predicted and actual
		 double error = helm_msh_pt()->node_pt(i)->time_stepper_pt()->
      			 temporal_error_in_value(helm_msh_pt()->node_pt(i),0);
		 //Add the square of the individual error to the global error
		 global_error += error*error;
	  } 
	 // Divide by the number of nodes in both meshes
	 global_error /= double(n_node_density + n_node_helm);
	
	 // Return square root...
	 return sqrt(global_error);
		
}// end of global_temporal_error_norm

//======================================Start of doc_solution =============================================
//document the solution
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::doc_solution()
{
	 ofstream some_file;
	 char filename[100];	
	 oomph_info << std::endl;
	 oomph_info << "=================================================" << std::endl;
	 oomph_info << "Docing solution for t=" << time_pt()->time() << std::endl;
	 oomph_info << "=================================================" << std::endl;
		
	 oomph_info << " Timestep: " << Doc_info.number() << std::endl;
		
	 // Output the density to file 
	 //-----------------
	 sprintf(filename,"%s/soln_density%i_on_proc%i.dat",Doc_info.directory().c_str(),
	         Doc_info.number(), this->communicator_pt()->my_rank());
	 some_file.open(filename);
	 density_msh_pt()->output(some_file);
 	 some_file.close();
	 oomph_info <<"Density written to file\n";
 	 // Output the AI and 1st moment to file
	 //-----------------
	 sprintf(filename,"%s/soln_helm%i_on_proc%i.dat",Doc_info.directory().c_str(),
	         Doc_info.number(), this->communicator_pt()->my_rank());
	 some_file.open(filename);
	 helm_msh_pt()->output(some_file);
 	 some_file.close();
	 oomph_info <<"AI and first moments written to file\n";
	 
	 //increment the doc counter
	 Doc_info.number()++;
	 }// end of doc_solution

//======================================Start of doc_restart =============================================
//write restart file
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::doc_restart(){

	ofstream some_file;
       	char filename[100];		
	oomph_info << std::endl;
	oomph_info << "=================================================" << std::endl;
	oomph_info << "Writing restart file." <<  std::endl;
	oomph_info << "=================================================" << std::endl;

	//output restart data to file	
	sprintf(filename,"%s/restart%i_on_proc%i.dat",Restart_info.directory().c_str(),
			 Restart_info.number(), this->communicator_pt()->my_rank()); 
	some_file.open(filename);
	dump_it(some_file);
	some_file.close();
	oomph_info << "Restart data written\n";

	//increment the doc counter
	Restart_info.number()++;
}

// dump the solution to disk to allow for restart
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::dump_it(ofstream& dump_file)	{
	//write the predicted next time step
	dump_file << Next_dt << " # suggested next timestep" << endl;
	//write the index of the documented solutions for future restart
	dump_file << Doc_info.number() << " # output soln index" << endl;
	//write the index of the restart file for future restarts
	dump_file << Restart_info.number() << " # current restart index" << endl;

	// call generic dump()
	Problem::dump(dump_file);
}

// function that outputs important problem parameters
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::output_parameters_to_file(ofstream& param_file)	{
	param_file << "a = "<< a <<"\nL = " << L <<"\nK = " << K <<"\nlambda = " << lambda
		<< "\nbeta = " << beta_decay << "\nalpha = " << alpha << "\nD2 = " << D2 <<"\nD1_star = " 
		<< D1_star << "\nD1_p_star = " << D1_p_star << "\nD_inf_0 = " <<
		D_inf_0 << "\nrho_star = " << rho_star << "\nepsilon = " << epsilon_reg << endl;
}

//===========================End of problem class definitions ==============================

// Check that input is valid and extract parameters from file
// there should be minimum two command line arguments, the first specifies the parameter 
// values in a file,the second specifies the output file directory
int read_input_params_from_file()	{
	ifstream infile; //input filestream
	if(CommandLineArgs::Argc<=2)	{// invalid number of inputs
		std::ostringstream error_stream;
		error_stream 
			<< "ERROR, incorrect number of command line arguments, usage is ./unstead_density_md ";
		error_stream << "param_filename result_directory " << endl << endl
			<< "Optional flags are: --restart-file [restart_filename_prefix] " << endl;
	     
	     throw OomphLibError(
	      error_stream.str(),
	      OOMPH_CURRENT_FUNCTION,
	      OOMPH_EXCEPTION_LOCATION);
	}
	// read in parameters from file and assign them appropriately
	infile.open(CommandLineArgs::Argv[1]);
	std::string buffer;
	unsigned n_params = 22, num_params_assigned=0; // number of parameters that need to be assigned
	oomph_info << "Assigning parameter values...\n";
	while(std::getline(infile,buffer))	{//read in one line
		std::string param_name, param_value;
		//ignore lines that don't contain `=`
		size_t pos = buffer.find('=');
		if(pos != std::string::npos)	{ // found '=' in line
			param_name = buffer.substr(0,pos); // left-hand side 
			param_value = buffer.substr(pos + 1, std::string::npos); //right-hand side
			double val = atof(param_value.c_str());

			//remove whitespaces from parameter name
			param_name.erase(remove(param_name.begin(), param_name.end(),' '), param_name.end());
			if(param_name=="a")		{a = val; 	 	}
			else if(param_name=="L")	{L = val;  		}
			else if(param_name=="K")	{K = val;  		}
			else if(param_name=="lambda")	{lambda = val;  	}
			else if(param_name=="D2")	{D2 = val;  		}
			else if(param_name=="beta")	{beta_decay = val;  	}
			else if(param_name=="alpha")	{alpha = val;  		}
			else if(param_name=="rho_star")	{rho_star = val;  	}
			else if(param_name=="D1_star")	{D1_star = val;		}
			else if(param_name=="D1_p_star"){D1_p_star = val;	}
			else if(param_name=="D_inf_0")	{D_inf_0 = val;		}
			else if(param_name=="epsilon")	{epsilon_reg = val;  	}
			else if(param_name=="Lu")	{Lu = val;		}
			else if(param_name=="Lx")	{Lx = val;		}
			else if(param_name=="Ly")	{Ly = val;		}
			else if(param_name=="epsilon_t"){epsilon_t = val;	}
			else if(param_name=="t_max")	{t_max = val;		}
			else if(param_name=="dt_out")	{dt_out = val;		}
			else if(param_name=="Nu_density"){num_el_density_u = val;}
			else if(param_name=="Nu_helm")	{num_el_helm_u = val;	}
			else if(param_name=="Nx")	{num_el_x = val;	}
			else if(param_name=="Ny")	{num_el_y = val;	}
			else{
				oomph_info << "Error, could not assign parameter '" << param_name <<"'\n";
				return 0;// should never get here
			}
			oomph_info << param_name << " = " << val << endl;
			num_params_assigned ++;
		}
	
	} // end of file reached
	infile.close();
	if(num_params_assigned != n_params) {
		oomph_info << "Error, some parameters left unspecified in input file '" << CommandLineArgs::Argv[1] <<"'\n";
		oomph_info << "Number of parameters set:" << num_params_assigned << "\nNumber of parameters needed:" << n_params << endl;
		return 0;
	}
	oomph_info <<"===============================\n";
	return 1;  
}

// computes the intracellular concentration at which the kinetics is in equilibrium with spatially homogeneous AI c=ce. i.e. f(u_star,ce)=0
void compute_ustar_homogeneous()
{
	double ce_quad_B = K - rho_star*alpha*(a + L)/(beta_decay*lambda);
	double ce_quad_C = -a*alpha*rho_star*K/(beta_decay*lambda);
	double ce = 0.5*(-ce_quad_B + std::sqrt(ce_quad_B*ce_quad_B - 4.0*ce_quad_C)); //homogeneous equilibrium value of AI
	u_star = (a + L*ce/(K+ce))/lambda; 
	oomph_info << "Equilibrium AI concentration in the uniform state: ce = " << ce << endl;
}

////////////////////////////////////////////////////////////////////////

//======start_of_main==================================================
/// Driver for Moment Equations Problem
//=====================================================================
int main(int argc, char* argv[])
 {
	#ifdef OOMPH_HAS_MPI
	 MPI_Helpers::init(argc, argv);
	 unsigned nproc = MPI_Helpers::communicator_pt()->nproc();
	 cout << "Code has been compiled with mpi support \n" 
		 << "Running on " << nproc <<" processors." << endl;
	#else
	 cout << "Code compiled WITHOUT mpi support." << endl;
	#endif
	
	//Store command line arguments	
	CommandLineArgs::setup(argc, argv);

	//Define possible command line arguments and parse the ones that were specified
	
	//name of restart file
	CommandLineArgs::specify_command_line_flag("--restart-file", &GlobalParams::Restart_filename);
	
	//name of file specifying problem distribution (for MPI)
	std::string partition_filename = "";
	CommandLineArgs::specify_command_line_flag("--partition-file", &partition_filename);

	//parse the input
	CommandLineArgs::parse_and_assign();

	//switch off output modifier
	oomph_info.output_modifier_pt() = &default_output_modifier;

	// Define processor-labeled output file for all on-screen output
	std::ofstream output_stream;
	char proc_out_filename[100];
	sprintf(proc_out_filename, "%s//output.%i", CommandLineArgs::Argv[2], MPI_Helpers::communicator_pt()->my_rank());
	output_stream.open(proc_out_filename);
	oomph_info.stream_pt() = &output_stream;
	OomphLibWarning::set_stream_pt(&output_stream);
	OomphLibError::set_stream_pt(&output_stream);
		
	//read parameters from input file and exit upon any errors
	if(!read_input_params_from_file())	{
		std::ostringstream error_stream;
		throw OomphLibError(
				error_stream.str(),OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
	}

	//compute the intracellular concentration at which the kinetics are in equilibrium with the spatially homogeneous AI concentration
	compute_ustar_homogeneous();
	//compute_linear_kinetic_terms();
 	oomph_info << "Critical intracellular concentration in the homogeneous state u_*=" << u_star << endl;	 

	//Build the problem 
  	FullDensityProblem<
		DensityHelmholtzElement2D<
					DensityElement2D,
	       				StratifiedHelmholtzElement2D>,
		HelmholtzDensityElement2D<
					StratifiedHelmholtzElement2D,
					DensityElement2D>
		> 
		problem(num_el_density_u,num_el_helm_u, num_el_x, num_el_y);
	
	// Set IC
	srand(static_cast<unsigned>(time(NULL)));	//seed for random number generator, this is needed for initial condition
	problem.set_initial_condition();	//assign IC

#ifdef OOMPH_HAS_MPI
	// ============= distribute the problem ===============
	Vector<unsigned> used_element_partition;	//vector containing problem distribution
	unsigned n_used = problem.mesh_pt()->nelement();
	used_element_partition.resize(n_used);
	bool report_stats = false;	//indicates whether to document partitioning stats
		
	// if partition file has been specified, read in distribution from file and then distribute the problem
	if(CommandLineArgs::command_line_flag_has_been_set("--partition-file"))	{
		oomph_info << "Partition file specified, opening for restart...";
		//read in partition from file
		std::ifstream partition_file;
		partition_file.open(partition_filename);
		if(partition_file.is_open())	{
			string input_string;
			oomph_info <<"... done.\n";
			//read  first lines which contains number of elements in partition.... we just ignore this line
			getline(partition_file, input_string, '\n');
			for(unsigned e=0; e<n_used; e++)	{
				getline(partition_file, input_string, '\n');	//get up to terminating comment
				used_element_partition[e] = unsigned(atof(input_string.c_str()));	//read in step number
			}
			partition_file.close();
			//distribute
			problem.distribute(used_element_partition, report_stats);
		}
		else	{
			oomph_info << "Could not open" << partition_filename << "for restart"<< endl;
		}
	}
	//if no partition file specified, distribute using default partition and then document the partition
	else	{
		used_element_partition = problem.distribute(report_stats);	//distribute

		//doc partition
		char new_partition_filename[100];
		sprintf(new_partition_filename, "%s/partition.dat",CommandLineArgs::Argv[2]);
		std::ofstream partition_file;	//partition
		partition_file.open(new_partition_filename);
		//write to file
		partition_file << n_used << std::endl;
		for(unsigned e=0; e<n_used; e++)	{
			partition_file << used_element_partition[e] << std::endl;
		}
		partition_file.close();
	}
	
	oomph_info <<"Finished distribution of the problem\n";

	// ====================================================
#endif 
	// open file for writing output times to
	ofstream t_out_file;
	char t_out_filename[100];
	sprintf(t_out_filename, "%s/t_out.dat",CommandLineArgs::Argv[2]);
	t_out_file.open(t_out_filename, ios::app);
		 
	//restart from file if restart flag has been specified 
	if(CommandLineArgs::command_line_flag_has_been_set("--restart-file"))	{
			problem.restart();
			oomph_info << "Finished restart\n";
	}	
	else	{	//if not restarting, then doc the initial condition
		//Output initial condition
		problem.doc_solution();
		if(MPI_Helpers::communicator_pt()->my_rank() == 0)	{
			t_out_file << problem.time_pt()->time() << endl; // output initial time to file
		}
	}
	double dt=problem.next_dt();
	
	// times to record solution at
	double t_out = std::max(problem.time_pt()->time(), dt_out);

	// write problem parameters to file
	ofstream param_file;
	char param_filename[100];
	sprintf(param_filename, "%s/params.dat",CommandLineArgs::Argv[2]);
	param_file.open(param_filename);
	problem.output_parameters_to_file(param_file);	//write equation-specific quantities

	//now write numerical quantities
	param_file << "dt_output = " << dt_out << "\nnum_el_density_u = " << num_el_density_u 
	 <<"\n num_el_helm_u = " << num_el_helm_u << "\nnum_el_x = " << num_el_x 
	 <<"\nnum_el_y = " << num_el_y << "\nepsilon_t = " << epsilon_t << endl;
	param_file.close();

	//setup counters for timestepping loop
	unsigned steps_since_write = 0; //number of steps since the last time the restart file was written
	unsigned steps_btw_restarts= 100;	//write restart file every steps_btw_restarts timesteps 
	unsigned steps = 0;	//number of steps taken thus far
	unsigned max_steps = 5000;	//maximum number of steps before terminating
	// Timestepping loop: Don't know how many steps we're going to take in advance
	while (problem.time_pt()->time()<t_max && steps < max_steps)
	 {
	  
	   // Take an adaptive timestep -- the input dt is the suggested timestep.
	   // The adaptive timestepper will adjust dt until the required error
	   // tolerance is satisfied. The function returns a suggestion
	   // for the timestep that should be taken for the next step. This
	   // is either the actual timestep taken this time or a larger
	   // value if the solution was found to be "too accurate". 
	   double dt_new=problem.adaptive_unsteady_newton_solve(dt,epsilon_t);
	
	   // Use dt_new as suggestion for the next timestep
	   dt=dt_new;

	   //store for restart
	   problem.next_dt() = dt_new;
	   oomph_info << "Current time = " << problem.time_pt()->time() << std::endl;
	       
	   //output solution at specified times
	   if(problem.time_pt()->time()>t_out)
	   {
		   problem.doc_solution();
		   t_out = problem.time_pt()->time() + dt_out;
		   //record the time in a separate output file
		   if(MPI_Helpers::communicator_pt()->my_rank() == 0)	{ //record only once
			   t_out_file << problem.time_pt()->time() << endl;
		   }

		//write a restart file if enough timesteps have been taken
		   if(steps_since_write>steps_btw_restarts)	{
			   problem.doc_restart();
			   steps_since_write=0;	//update counters
		   }
		   oomph_info << " Current Doc_info.number(): " << problem.doc_number() << endl;
	   }  
	   //increment number of steps since last restart file was written
	   steps_since_write++;
	   
	   //update the count of steps taken
	   steps++;
	   if(steps > max_steps)	{
		   oomph_info << "Reached maximum number of steps... terminating" << endl;
	   }
	  } // end of timestepping loop

	//write one last restart file
	problem.doc_restart();

	#ifdef OOMPH_HAS_MPI	
	 MPI_Helpers::finalize();
	#endif
	 
} // end of main
