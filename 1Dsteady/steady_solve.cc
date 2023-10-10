//Driver for steady 1D chemically structured model with canonical QS kinetics
//Solves the coupled equations
// 0 = D(u)n_{xx} n + epsilon*n_{uu} - (f(u,c)n)_u
// 0 = Dc c_{xx} - beta*c + alpha_0*ubar
// ubar = int_0^\infty u*n du
// \int_0^\infy n du = rho_*

// Generic oomph-lib routines
#include "generic.h"
#include "StratifiedHelmholtzElement.h"
#include "DensityElement1D.h"
#include "meshes/rectangular_quadmesh.template.h" //include the mesh

using namespace std;
using namespace oomph;

//==start_of_namespace================================================
/// Namespace for problem parameters
//====================================================================
namespace GlobalParams
{
	//=====================Basic Problem parameters==============================
	// domain parameters
	double Lu, Lx; // length of domain in u- and x-directions, respectively (which correspond to x,y- directions in oomph-lib)
	double u_min_mesh;	//left boundary of the domain

	// bulk parameters
	double D2, beta_decay, alpha;

	//kinetic parameters
	double a, L, K, lambda;
	
	//target cell density and population to achieve in steady-state calculation, actual solution will have a different value, we use this only to calculate u_star
	double rho_target;
	double total_population;

	// regularization parameter
	double epsilon_reg;

	//swimming parameters
	double D1_star, D1_p_star;	//effective diffusion coefficient and it's derivative at u=u*
	double D_inf_0; //diffusion coefficient at u=0 or u->inf
	
	//=======================Numerical parameters=================================================
	//number of elements in each coordinate direction
	unsigned num_el_density_u, num_el_helm_u, num_el_x;

	double rho_eval; 	//value of x at which rho(x) is evaluated for the bifurcation diagram

	//name of the restart file (if used)
	std::string Restart_filename = "";
	//=======================Derived Quantities in terms of Basic params===========================
	double u_star;
	//=============================================================================================

	// pointer to the single-valued data that gets replaced with the mass constraint
	Data* traded_density_data_pt;

	//initial condition 
	void get_IC(const double& time, const Vector<double>& x, Vector<double>& sol)	{
		// ========== Gaussian distribution, spatially uniform =============

		double epsilon_IC = epsilon_reg; // width of the gaussian to use in the initial condition (cannot use with the unregularized problem, i.e epsilon_reg=0)
		double gaussian = std::exp(-0.5*lambda/epsilon_IC*std::pow(x[0]-u_star,2)); //gaussian factor that shows up in the steady-state and perturbation
		double normalization = std::sqrt(lambda/(2*M_PI*epsilon_IC));
		double ss = rho_target*normalization*gaussian;	//steady-state concentration
		
		//We need an initial starting point on the MIPS branch in order to calculate this branch. We don't know
		//analytically what the solutions on this branch are,but we can make a reasonable guess using a sigmoid function.
		//A more systematic approach is to use the eigenfunctions from the linear stabiltiy analysis, but in practice this can require very 
		//fine tuning of the bifurcation parameter when the bifurcating branch is almost supercritical

		double sigmoid_x = exp(2.0*(x[1]-5.0))/(1.0 + 1.0*exp(2.0*(x[1] - 5.0)));	//sigmoid function of the class e^x/(1 + e^x).

		double ustar_mode = 2*u_star*sigmoid_x;
		gaussian = std::exp(-0.5*lambda/epsilon_IC*std::pow(x[0] - ustar_mode,2)); //gaussian factor that shows up in the steady-state and perturbation
		normalization = std::sqrt(lambda/(2*M_PI*epsilon_IC));
		ss = rho_target*normalization*gaussian;	//steady-state concentration

		sol[0] = (ss/rho_target)*sigmoid_x*62.5/3.0;
		sol[1] = 35.0*sigmoid_x;
		sol[2] = u_star*rho_target*sigmoid_x/3.0;
	}

	// function for evaluating the effective diffusion coefficient
	void get_D1(const Vector<double>& x, double* D1)	{
		// linear function
		double u1 = (D_inf_0 - D1_star)/D1_p_star + u_star;
		if(D1_p_star < 0)	{ 	//decreasing D1
			if(x[0] > u1)	{
				*D1 = D_inf_0;	//constant
			}
			else	{
				*D1 = D1_p_star*(x[0]-u_star) + D1_star;	//linear
			}
		}
		else	{	//increasing D1
			if(x[0] > u1)	{
				*D1 = D1_p_star*(x[0] - u_star) + D1_star;	//linear
			}
			else	{
				*D1 = D_inf_0;	//constant
			}
		}
	}
} // end of namespace

using namespace GlobalParams;

// A class that implements an element responsible for enforcing the mass 
// constraint that the integral of the density solution must be equal to a specified value
class MassConstraintElement : public GeneralisedElement	{
private:
	//pointer to target mass of the solution
	double* Prescribed_mass_pt;
	//the data containing the traded density dof is stored as external data for 
	//this element. store it's index for accessing via external_data_pt(...)
	unsigned External_data_index_of_traded_density;
	//index of the value in the traded density data that corresponds to the traded density value
	unsigned Index_of_traded_density_value;
	
	//local eqn number for the traded density
	inline int density_traded_local_eqn()	{
		return this->external_local_eqn(External_data_index_of_traded_density,
				Index_of_traded_density_value);
	}	
	//function that fills in the residual corresponding to the volume constraint (only part of it
	//is taken into account here, the rest are in the corrsponding MassComputationElements)
	//contribution to jacobian is calculated if flag is set
	void fill_in_generic_contribution_to_residuals_mass_constraint(Vector<double>& residuals)	{
		const int local_eqn = this->density_traded_local_eqn();
		if(local_eqn>=0)	{
			residuals[local_eqn] -= *Prescribed_mass_pt;
		}
	}
public:
	//Constructor
	MassConstraintElement(double *prescribed_mass_pt, 
			Data* density_traded_data_pt, const unsigned& index_of_traded_density)	{
		//pointer to the prescribed mass in the solution
		Prescribed_mass_pt = prescribed_mass_pt;
		
		External_data_index_of_traded_density = add_external_data(density_traded_data_pt);
		Index_of_traded_density_value = index_of_traded_density;
	}

	//Destructor
	~MassConstraintElement()  {}
	//access function to data that contains the traded density dof
	inline Data* density_traded_data_pt()	{
		return external_data_pt(External_data_index_of_traded_density);
	}
	//access function to the traded density value
	inline double density_traded()	{
		return density_traded_data_pt()->value(Index_of_traded_density_value);
	}
	//access function to the index at which the traded density is stored
	inline unsigned index_of_traded_density()	{
		return Index_of_traded_density_value;
	}

	//fill in residuals for the mass constraint
	void fill_in_contribution_to_residuals(Vector<double>& residuals)	{
		this->fill_in_generic_contribution_to_residuals_mass_constraint(residuals);
	}
	// fill in contribution to jacobian (no contribution from the master element) and residuals
	void fill_in_contribution_to_jacobian(Vector<double>& residuals,
			DenseMatrix<double>& jacobian)	{
		//no contribution to jacobian from this element as it only subtracts the target population
		this->fill_in_generic_contribution_to_residuals_mass_constraint(residuals);
	}
};

/* Class that implements the elements needed for computing the solution mass for enforcing the mass 
 * constraint in the density equation
 */
template<class DENSITY_ELEMENT>
class MassComputationElement : public virtual GeneralisedElement	{
private:
	//index required to access the traded density data from external_data_pt(...)
	unsigned Data_number_traded_density;
	//index of the value in the traded density data that corresponds to the actual
	//traded density dof
	unsigned Index_of_traded_density_value;

	//vector for storing indices of external density data attached to this element (this external data is added in the constructor
	//and is not the traded density data)
	Vector<unsigned> External_local_density_dof;

	//density element that this element is responsible for computing the mass of
	DENSITY_ELEMENT* Density_el_pt;

	//local eqn number for the traded density
	inline int density_traded_local_eqn()	{
		return this->external_local_eqn(Data_number_traded_density, 
				Index_of_traded_density_value);
	}
	//calculates elemental contribution to the residuals and, if flag is set, the jacobian
	void fill_in_generic_residual_contribution_mass_compute(
		Vector<double>& residuals, DenseMatrix<double>& jacobian, bool flag)	{

		// number of nodes in the element
		const unsigned n_node = Density_el_pt->nnode();
		//index of density values in the density data
		const unsigned density_ind = Density_el_pt->field_index_density();
		//allocate memeory for shape functions
		Shape psi(n_node);
		DShape dpsi(n_node,2); // dummy storage for calling dshape_eulerian

		// memory for local coords
		Vector<double> s(2);
				
		const unsigned n_intpt = Density_el_pt->integral_pt()->nweight();
		// equation number of traded density
		const int local_eqn = this->density_traded_local_eqn();
		//storage for local dof in the density element 
		int local_dof;
	
		// loop over integration points
		for(unsigned ipt=0; ipt<n_intpt; ipt++)	{
			//local coordinates of integration point
			s[0] = Density_el_pt->integral_pt()->knot(ipt,0); // arguments are (i,j) where i is the index of the integration poitn and j is the jth coordinate
			s[1] = Density_el_pt->integral_pt()->knot(ipt,1);
	
			//get the integral weight
			double w = Density_el_pt->integral_pt()->weight(ipt);
			//call the derivatives of the shape and test
			double J = Density_el_pt->dshape_eulerian(s, psi, dpsi);
			double W = w*J;
			// if not a Dirichlet boundary, then add contribution
			if(local_eqn >=0)	{	
				// Calculate the interpolated density 
				double interpolated_n = 0.0;
				for(unsigned j=0; j<n_node; j++)	{
					interpolated_n += Density_el_pt->field_value_density(j)*psi[j];
				}
				// add contribution to total mass to the residual
				residuals[local_eqn] += interpolated_n*W;
				if(flag)	{	//calculate jacobian 
					for(unsigned j2=0; j2<n_node; j2++)	{
						local_dof = external_local_eqn(External_local_density_dof[j2], density_ind);
						if(local_dof >=0)	{
							//add elemental contribution to jacobian
							jacobian(local_eqn, local_dof) += psi[j2]*W;
						}
					}
				}
			}
		}//end of outer loop over integration points

		//The above jacobian calculation works for all but one of the entries in the global jacobian. This occurs if
		//the traded density data is actually contained in the attached density element. This will happen for exactly
		//one element in the mesh and for exactly one data value in that element. 
		//Find the corresponding entry in the jacobian and put in correct value
		if(flag && local_eqn>=0 )	{
			//go through all external density data
			for(unsigned i=0; i<n_node; i++)	{
				unsigned ith_density_dof = External_local_density_dof[i];
				//check if the external density data and the traded density data point to the same value
				if(external_data_pt(ith_density_dof)->value_pt(density_ind) 
						== external_data_pt(Data_number_traded_density)->value_pt(density_ind))	
				{
					local_dof = external_local_eqn(External_local_density_dof[i], density_ind);
					if(local_dof >=0)	{
						//since the two data pointers point to the same value, their entries in the jacobian are equal
						jacobian(local_eqn, local_eqn) = jacobian(local_eqn, local_dof);	
					}
				}
			}
		}

	}//end of fill_in_generic_residual_contribution(...)

	//computes the mass density on this element, each element has three values corresponding to 
	//the three different nodal values of x in the density element
	//Provides access to the density for post-processing purposes
	void compute_elemental_mass_density(Vector<double>& mass_density)	{
		// number of nodes in the element
		const unsigned n_node = Density_el_pt->nnode();
		
		//distance between nodes in u-direction
		double du = 0.5*Lu/num_el_density_u; // distance between nodes in each element (rectangular quadratic elements only)
		//memory for local coords
		Vector<double> s(2);

		//loop over nodes and calculate the mass with 1D simpsons rule
		for(unsigned j=0; j< n_node; j++)	{
			//get the nodal coordinates
			Density_el_pt->local_coordinate_of_node(j,s);
			//get value of n stored at this node
			double nval = Density_el_pt->field_value_density(j);
			//which mass value do we add onto
			unsigned mass_index = int(s[1] + 1.5);
			if(mass_index > 2 || mass_index < 0)	{
				throw 
					OomphLibError("Trouble! Evaluation of local mass density contribution failed. Integral index out of bounds..\n",
							"FullDensityProblem::compute_mass_density()",
							OOMPH_EXCEPTION_LOCATION);
			}
			//add contribution to integral with Simpson's rule
			unsigned u_ind =  unsigned(s[0] + 1.5);
			if(u_ind == 0 || u_ind == 2)	{
				mass_density[mass_index]+=nval*du/3.0;
			}
			else if(u_ind==1)	{
				mass_density[mass_index]+=nval*du*4.0/3.0;
			}
			else	{
				throw 
					OomphLibError("Trouble! Evaluation of local mass density contribution failed. Simpons rule index out of bounds..\n",
							"FullDensityProblem::compute_mass_density()",
							OOMPH_EXCEPTION_LOCATION);
			}
		}
	}
public:
	//constructor: the density data passed gets added as external data 
	MassComputationElement(DENSITY_ELEMENT* const& density_el_pt) {
		//add all density element data in the external element as external data
		unsigned dof;
		unsigned n_node = density_el_pt->nnode();
		for(unsigned j=0; j<n_node; j++)	{
			dof = this->add_external_data(density_el_pt->node_pt(j));
			External_local_density_dof.push_back(dof);
		}
		//store a pointer to this element
		Density_el_pt = dynamic_cast<DENSITY_ELEMENT*>(density_el_pt);
	}
	
	//fill in residuals
	void fill_in_contribution_to_residuals(Vector<double>& residuals)	{
		this->fill_in_generic_residual_contribution_mass_compute(residuals, GeneralisedElement::Dummy_matrix, 0);
	}
	//fill in jacobian and residuals
	void fill_in_contribution_to_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian)	{
		//compute analytical jacobian
		this->fill_in_generic_residual_contribution_mass_compute(residuals, jacobian, 1);
		//compute jacobian with fd 
		//GeneralisedElement::fill_in_contribution_to_jacobian(residuals, jacobian);
	}

	//sets the "master" mass control element
	void set_mass_constraint_element(
			MassConstraintElement* const& mass_constraint_el_pt)	{
		//added traded density pointer as external data
		Data_number_traded_density = 
			this->add_external_data(mass_constraint_el_pt->density_traded_data_pt());
		Index_of_traded_density_value = mass_constraint_el_pt->index_of_traded_density();
	}

	//interpolates the mass density at a given value of x  using 1d shape functions
	double interpolate_elemental_mass_density(double x)	{
		//number of nodes in density element
		unsigned n_node = Density_el_pt->nnode();
		//leftmost value of x in the density element
		double xmin = 1e10;	//(dummy)
		for(unsigned j=0; j<n_node; j++)	{
			xmin = min(Density_el_pt->node_pt(j)->x(1), xmin);
		}
		//calculate local coordinates  of the given x
		double dx = Lx/num_el_x;
		double xlocal = 2.0*(x - xmin)/dx - 1.0;	//local coorindate between -1 and 1
		
		//evaluate 1D shape functions at local coordinates
		Vector<double> one_d_shape(3);
		one_d_shape[0] = xlocal*(xlocal - 1.0)/2.0;
		one_d_shape[1] = -(xlocal - 1.0)*(xlocal + 1.0);
		one_d_shape[2] = xlocal*(xlocal + 1.0)/2.0;

		//allocate and initialize the mass densities on this element
		Vector<double> mass_density(3);
		for(unsigned j=0; j<3; j++)	{mass_density[j]=0.0;	}
		//compute the mass densities
		compute_elemental_mass_density(mass_density);
		//interpolate the mass densities using quadratic interpolation (same as shape functions)
		double interpolated_mass_density = 0.0;
		for(unsigned j=0; j<3; j++)	{
			interpolated_mass_density += mass_density[j]*one_d_shape[j];
		}
		return interpolated_mass_density;	
	}
};

// Start of Element class for enforcing non-local auxillary equation
// requires a pointer to the mesh containing all the density elements so that the non-local residuals can be calculated
template<class DENSITY_ELEMENT>
class UbarElement : public virtual GeneralisedElement
{
public:
	//constructor for data created internally (within the UbarElement)
	UbarElement(RectangularQuadMesh<DENSITY_ELEMENT>* density_msh_pt, const unsigned row_ind_x) 
		: Density_msh_pt(density_msh_pt), //this line initializes the protected member variable Density_msh_pt equal to the argument density_msh_pt
		  Row_ind_x(row_ind_x)
	{
		//create new data object that will contain the unknowns on this element
		Ubar_data_pt = new Data(3);
		//declare this data as internal data for this element so that the current element is `in charge' of it
		Ubar_data_index = add_internal_data(Ubar_data_pt);
		Ubar_data_created_internally = true;
		
		//We now need to pass (pointers to) the external density data which are required in the residual and jacobian calculations
		unsigned ext_el_ind, n_node;	//index of element containing data to be made external data of this element, number of nodes in the element
		unsigned n_element_u = Density_msh_pt->nx();	//number of elements in u direction
		for(unsigned i=0; i<n_element_u; i++)	{
			ext_el_ind = n_element_u*Row_ind_x + i;	//calculate index of external density element
			DENSITY_ELEMENT* ext_el_pt = dynamic_cast<DENSITY_ELEMENT*>(Density_msh_pt->element_pt(ext_el_ind));
			n_node = ext_el_pt->nnode();
			//loop over nodes in the external element and add external data for each
			for(unsigned j=0; j<n_node; j++)	{
				this->add_external_data(ext_el_pt->node_pt(j));
			}
		}
		//store dimensions of current element and index of density value in the density element
		ext_el_ind = n_element_u*Row_ind_x;
		DENSITY_ELEMENT* ext_el_pt = dynamic_cast<DENSITY_ELEMENT*>(Density_msh_pt->element_pt(ext_el_ind));
		n_node = ext_el_pt->nnode();
		Density_external_ind = ext_el_pt->field_index_density();
		for(unsigned j=0; j<n_node; j++)	{
			Umin = min(Umin, ext_el_pt->node_pt(j)->x(0));
			Xmin = min(Xmin, ext_el_pt->node_pt(j)->x(1));
		}
	}

	//constructor for data created externally (outside the UbarElement)
	UbarElement(RectangularQuadMesh<DENSITY_ELEMENT>* density_msh_pt, const unsigned row_ind_x, Data* ubar_data_pt)	
		: Density_msh_pt(density_msh_pt), 
		  Ubar_data_pt(ubar_data_pt),
		  Row_ind_x(row_ind_x)
	{
		//declare data passed to constructor as external data so that the problem is `in charge' of it
		Ubar_data_index = add_external_data(Ubar_data_pt);
		Ubar_data_created_internally = false;
		//We now need to pass (pointers to) the external density data which are required in the residual and jacobian calculations
		//number of elements in u direction
		unsigned ext_el_ind, n_node;	//index of element containing data to be made external data of this element, number of nodes in the element
		unsigned n_element_u = Density_msh_pt->nx();	
		for(unsigned i=0; i<n_element_u; i++)	{
			unsigned ext_el_ind = n_element_u*Row_ind_x + i;
			DENSITY_ELEMENT* ext_el_pt = dynamic_cast<DENSITY_ELEMENT*>(Density_msh_pt->element_pt(ext_el_ind));
			unsigned n_node = ext_el_pt->nnode();
			//loop over nodes in the external element and add external data for each. Also store the index required to 
			//obtain the data from external_data_pt(...)
			for(unsigned j=0; j<n_node; j++)	{
				this->add_external_data(ext_el_pt->node_pt(j));
			}
		}
		//store dimensions of current element and index of density value in the density element
		ext_el_ind = n_element_u*Row_ind_x;
		DENSITY_ELEMENT* ext_el_pt = dynamic_cast<DENSITY_ELEMENT*>(Density_msh_pt->element_pt(ext_el_ind));
		n_node = ext_el_pt->nnode();
		Density_external_ind = ext_el_pt->field_index_density();
		for(unsigned j=0; j<n_node; j++)	{
			Umin = min(Umin, ext_el_pt->node_pt(j)->x(0));
			Xmin = min(Xmin, ext_el_pt->node_pt(j)->x(1));
		}

	}

	unsigned ndof_types() const	{
		return 1;
	}

	//call default routine for eqn numbering as well as the local routine for storing eqn numbers of the ubar eqns
	void assign_all_generic_local_eqn_numbers(const bool& store_local_dof_pt)	{
		GeneralisedElement::assign_internal_and_external_local_eqn_numbers(store_local_dof_pt);
		this->assign_additional_local_eqn_numbers();
	}
	
	//store local eqn number of the equation for ubar 
	void assign_additional_local_eqn_numbers()	{
		if(Ubar_data_created_internally)  {	//ubar is an internal variable
			for(unsigned j=0; j<3; j++)  {
				Ubar_local_eqn.push_back(internal_local_eqn(Ubar_data_index,j));
			}
		}
		else	{	//ubar is an external variable
			for(unsigned j=0; j<3; j++)	{
				Ubar_local_eqn.push_back(external_local_eqn(Ubar_data_index,j));
			}
		}
	}

	//wrapper function for elemental residual calculation
	void fill_in_contribution_to_residuals(Vector<double>& residuals)	{
		this->fill_in_generic_residual_contribution_ubar(residuals, GeneralisedElement::Dummy_matrix, 0);
	}

	//wrapper function for computing elemental jacobian 
	void fill_in_contribution_to_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian)	{
		//analytical jacobian 
		this->fill_in_generic_residual_contribution_ubar(residuals, jacobian, 1);
		//compute jacobian by fd
		//GeneralisedElement::fill_in_contribution_to_jacobian(residuals, jacobian);
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
			for(unsigned j=0; j<3; j++)	{
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

	//pointer to data object that contains the three values of ubar in this element
	Data* ubar_data_pt() const	{
		return Ubar_data_pt;
	}
private:
	//calculates elemental contribution to residuals and jacobian if flag is set
	void fill_in_generic_residual_contribution_ubar(Vector<double>& residuals, DenseMatrix<double>& jacobian, bool flag)	{
		Vector<double> first_moments(3);
		for(unsigned j=0; j<3; j++)	{
			first_moments[j]=0.0;
		}
	        compute_first_moment(first_moments);

		for(unsigned j=0; j<3; j++)	{
			if(Ubar_local_eqn[j] >= 0)	{	// should always be true, <0 corresponds to dirichlet boundary condition
				residuals[Ubar_local_eqn[j]] = Ubar_data_pt->value(j) - first_moments[j];
			}
		}
		//calculate jacobian if flag is set
		if(flag)	{
			//diagonal blocks
			this->fill_in_diagonal_blocks_jacobian_analytic_ubar(jacobian);
			//off-diagonal blocks
			this->fill_in_off_diagonal_blocks_jacobian_analytic_ubar(jacobian);
		}

	}
	//helper function for evaluating off-diagonal blocks of the elemental jacobian matrix. No integration loop as the integral
	//has been evaluated analytically assuming quadratic shape functions (i.e. use Simpson's rule)
	void fill_in_off_diagonal_blocks_jacobian_analytic_ubar(DenseMatrix<double>& jacobian)	{
		//storage for local eqn number and local unknown
		int local_eqn = 0, local_dof = 0;	
		// number of elements in u directions for density mesh
		unsigned n_element_density_u = Density_msh_pt->nx();
		unsigned n_element_x = Density_msh_pt->ny();
		
		double du = Lu/n_element_density_u; // length of each density element in u-direction
		double dx = Lx/n_element_x;
		double nodal_distance_u = du/2.0;  //distance between nodes in the u-direction
		double nodal_distance_x = dx/2.0;
		
		//loop over entries in the single internal data value
		for(unsigned i=0; i<3; i++)	{
			local_eqn = Ubar_local_eqn[i];
			double x_i = Xmin + nodal_distance_x*i;
			if(local_eqn >=0)	{
			//loop over external density data
			unsigned n_ext_data = this->nexternal_data();
			for(unsigned j=0; j<n_ext_data; j++)	{
				//cast external data to node so as to be able to access it's position in the mesh
				Node* ext_data = dynamic_cast<Node*>(external_data_pt(j));
				//global coordinates of node (value of u which is at the centre of the current shape function's support)
				double uc = ext_data->x(0);
				double xc = ext_data->x(1);
				//determine the horizontal node number of the current node
				unsigned current_node = unsigned((uc - u_min_mesh)/nodal_distance_u + 0.5 );	//extra 0.5 is to prevent truncation errors

			       	//determine the local degree of freedom from the external data
				local_dof = this->external_local_eqn(j, Density_external_ind);
				//add contribution to jacobian  
				//only nonzero entries are for density values at the same value of x as the current (ith) ubar value
				if(local_dof >=0 && abs(xc-x_i)<1e-3)	{
				// Add the contribution to the appropriate jacobian entry
					if(current_node%2==1) {	//middle shape function	
						jacobian(local_eqn, local_dof) = -2.0*uc*du/3.0;
					}
					else	{	//left or right shape function
						jacobian(local_eqn, local_dof) = - uc*du/3.0;
						//end points contribute only half
						if(current_node==0 || current_node==2*n_element_density_u)	{
							jacobian(local_eqn, local_dof) = 0.5*jacobian(local_eqn, local_dof);
						}
					}
				}
					
			}
			}
		}//end of loop over 3 internal data values
	}// end of fill_in_off_diagonal...

	//helper function for evaluating diagonal blocks of the elemental jacobian matrix, this block is simply an identity matrix with permuted rows
	void fill_in_diagonal_blocks_jacobian_analytic_ubar(DenseMatrix<double>& jacobian)	{
		int local_eqn = 0, local_unknown = 0;	//storage for local eqn number and local unknown number
		//loop over the 3 entries in the ubar data value of this element
		for(unsigned j=0; j<3; j++)	{
			local_eqn = Ubar_local_eqn[j];
			if(local_eqn >= 0)	{
				local_unknown = local_eqn;
				jacobian(local_eqn,local_unknown) = 1.0;
			}
		}
	}

	//helper function for residual calucation. computes the non-local part of the residuals
	void compute_first_moment(Vector<double>& first_moments)	{
		// number of elements in u directions for density mesh
		unsigned n_element_density_u = Density_msh_pt->nx();
		
		double du = 0.5*Lu/n_element_density_u; // distance between nodes in each element (rectangular quadratic elements only)
		// storage for integrals (one for each node at fixed x in the element)
		Vector<double> s(2); // storage for local coordinates in each element
		
		//Go through all elements with fixed x and sum contribution to integral from each
		//Note there there are nodes within each element that have different values of x
		for(unsigned j=0; j<n_element_density_u; j++)	{
			unsigned elem_ind = Row_ind_x*n_element_density_u + j;
			// upcast from generalisedElement to present type
			DENSITY_ELEMENT* el_pt = dynamic_cast<DENSITY_ELEMENT*>(Density_msh_pt->element_pt(elem_ind));
			unsigned num_node = el_pt->nnode();
			// loop through nodes at fixed x in each element
			for(unsigned k=0; k<num_node; k++)	{	 
				double nval = el_pt->field_value_density(k); // value of n stored at this node
				el_pt->local_coordinate_of_node(k, s); // local coordinates of kth node
				//global coordinates of node
				double uval = el_pt->node_pt(k)->x(0); 
				// Add the contribution to the appropriate integral at fixed x (s[1])
				int integral_ind = int(s[1] + 1.0 + 0.5); // add an additional 0.5 to avoid truncation errors (range of local coords is [-1,1])
				if(integral_ind > 2 || integral_ind < 0)	{
					throw 
						OomphLibError("Trouble! Evaluation of non-local integral failed. Integral index out of bounds..\n",
			                     	   "FullDensityProblem::update_moment_nodal_values()",
			                     	   OOMPH_EXCEPTION_LOCATION);
				}
				if(std::abs(s[0]) < 1e-2) {//middle column where s[0]=0  	
					first_moments[integral_ind] += 4.0/3.0*uval*nval*du; // contribution from mid-point of element according to Simpson's rule
				}
				else	{
					first_moments[integral_ind] += uval*nval*du/3.0; // contribution from end-points of element according to Simpson's rule
				}
			}
		}
	}// end of compute_first moment

	//evaluates the oneD shape functions for quadratic elements at local coordinate s. index is 0,1,2 corresponding
	//to left middle right shape functions
	void oneD_shape_fn(double* phi_1D, double s, int index )	{
		if(index ==0)	{
			*phi_1D = 0.5*s*(s - 1.0);
		}
		else if(index == 1)	{
			*phi_1D = 1.0 - s*s;
		}
		else	{
			*phi_1D = 0.5*s*(s + 1.0);
		}
	}


protected:
	
	bool Ubar_data_created_internally;
	//Mesh for density elements
	RectangularQuadMesh<DENSITY_ELEMENT>* Density_msh_pt;

	//Pointer to data item containing the value of ubar in this element
	Data* Ubar_data_pt;

	//index required to obtain Ubar_data_pt from internal_data_pt(...)
	unsigned Ubar_data_index;
	
	//row index corresponding to the value of x with which this element is associated, x = Row_ind_x*dx
	unsigned Row_ind_x;

	//index of density value stored in external data
	unsigned Density_external_ind;

	// Local equation number of the Ubar equation
	Vector<int> Ubar_local_eqn;

	//minimum values of the global coordinates in this element (dummy values)
	double Umin=1e12, Xmin=1e12;

}; // end of MassElement class


// =========================  Start of Multi-domain Upgrade of the Density Element ========================
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
class DensityHelmholtzElement1D
	: public virtual DENSITY_ELEMENT,
	public virtual ElementWithExternalElement
{
public:
	//Constructor: call underlying constructors
	DensityHelmholtzElement1D()
		:DENSITY_ELEMENT(), ElementWithExternalElement()	{
		// one interaction with other elements (in the kinetic term)
		this->set_ninteraction(1);
	}

	//overload the get_kinetics... function from the density element. This provides the coupling between the density and
	// helmholtz equations
	void get_kinetics_density(const unsigned& ipt, const Vector<double>& x, double* kinetics) const	{
		//interaction index
		unsigned interaction = 0;

		//get the AI concentration interpolated from the external element
		const double interpolated_c = dynamic_cast<HELM_ELEMENT*>(
				external_element_pt(interaction,ipt))->interpolated_c_strat_helmholtz(
				external_element_local_coord(interaction,ipt));
		// calculate the kinetics
		*kinetics = a + L*interpolated_c/(K + interpolated_c) - lambda*x[0]; 
	}

	// compute the elements jacobian matrix with finite differencing
	void fill_in_contribution_to_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian)		{
		//fill in by finite differencing
		//ElementWithExternalElement::fill_in_contribution_to_jacobian(residuals, jacobian);
		//fill in analytically
		//diagonal blocks
		DENSITY_ELEMENT::fill_in_contribution_to_jacobian(residuals, jacobian);
		//off-diagonal blocks
		this->fill_in_off_diagonal_jacobian_blocks_analytic_density(residuals, jacobian);

	}
	// helper function to compute off-diagonal jacobian blocks analytically
	void fill_in_off_diagonal_jacobian_blocks_analytic_density(Vector<double>& residuals, DenseMatrix<double>& jacobian)	{
	
		//solution index for n
		unsigned n_nodal_ind = this->field_index_density();

		//storage for local unknown and local equation
		int local_unknown = 0, local_eqn = 0;

		//number of nodes in this element
		unsigned n_node = this->nnode();

		//memory allocation for shape functions and derivatives
		Shape psi(n_node);
		DShape dpsi(n_node,2);

		//number of integration points
		const unsigned n_intpt = this->integral_pt()->nweight();
		
		//loop over integration points
		for(unsigned ipt=0; ipt<n_intpt; ipt++)	{
			
			//weight associated with this integration point
			double w = this->integral_pt()->weight(ipt);
			//derivatives of shape and test functions
			double J = this->dshape_eulerian_at_knot(ipt, psi, dpsi);
			// pre-multiply weights and jacobian
			double W = w*J;
			
			//the external elements' local eqn numbering map to global eqn numbers is different than for this element, 
			//global eqn numbers are used instead
			// vector of global eqn numbers corresponding to the external element's data
			Vector<unsigned> global_eqn_number_of_external_element_data;
			// vector containing all derivatives of c with respect to shape functions in it's own element (here c is in an external element)
			Vector<double> dinterpolated_c_ddata_strat_helmholtz;
			this->get_derivatives_wrt_external_element_data(ipt,
					dinterpolated_c_ddata_strat_helmholtz,
					global_eqn_number_of_external_element_data);

			//find out how many external element data there are
			const unsigned n_external_element_data = global_eqn_number_of_external_element_data.size();
			//Loop over shape functions and calculate n at the integration point 
			double interpolated_n = 0.0;
			for(unsigned j=0; j<n_node; j++)	{
				interpolated_n += this->field_value_density(j)*psi[j];
			}

			//interaction index
			unsigned interaction = 0;

			//get the AI concentration interpolated from the external element
			const double interpolated_c = dynamic_cast<HELM_ELEMENT*>(
				external_element_pt(interaction,ipt))->interpolated_c_strat_helmholtz(
				external_element_local_coord(interaction,ipt));

			// evaluate derivative of kinetics w.r.t c
			double fc = L*K/((K + interpolated_c)*(K + interpolated_c));

			//Assemble jacobian terms
			//loop over test functions
			for(unsigned j=0; j<n_node; j++)	{
				// Calculate derivatives of n residual w.r.t c, other off-diagonal blocks are zero
				local_eqn = this->nodal_local_eqn(j, n_nodal_ind);
				if(local_eqn>=0)	{ // if not a point on dirichlet boundary
					for(unsigned j2=0; j2<n_external_element_data; j2++)	{ // loop over external data
						//find local eqn number corresponding to global unknown
						local_unknown = this->local_eqn_number(
							global_eqn_number_of_external_element_data[j2]);
						if(local_unknown>=0)	{ // if locla unknown doesn't correspond to point on dirichlet boundary
							jacobian(local_eqn, local_unknown) += 
								-interpolated_n*dpsi(j,0)*fc*dinterpolated_c_ddata_strat_helmholtz[j2]*W;
						}
					}
				}
			}// end of outer loop over nodes
		}//end of loop over integration points

	}
//helper function for determining global equation numbers corresponding to data in external element
	void get_derivatives_wrt_external_element_data(const unsigned &ipt, Vector<double>& derivatives, 
			Vector<unsigned>& global_eqn_numbers)	{

		//dynamic cast external element to correct type
		unsigned interaction = 0;
		HELM_ELEMENT* source_el_pt =
			dynamic_cast<HELM_ELEMENT*>(external_element_pt(interaction, ipt));
		source_el_pt->dinterpolated_c_ddata_strat_helmholtz(
				external_element_local_coord(interaction, ipt), derivatives, global_eqn_numbers);
	}



};//end of DensityHelmholtzElement1D


 // =========================  Start of Multi-Domain upgrade of helmholtz Element ========================
template<class HELM_ELEMENT, class DENSITY_ELEMENT>	//technically the density element is not needed here since the interaction from ubar comes from external data
class HelmholtzDensityElement1D
	: public virtual HELM_ELEMENT,
	public virtual ElementWithExternalElement
{
private:
	//evaluates oneD shape functions at (1D) knot point ipt
	void oneD_shape_at_knot_helm(const unsigned& ipt, Vector<double>& psi_1D) const	{
		//storage for shape function
		unsigned ndof = 3;
		psi_1D.resize(ndof,0.0);
		
		double xloc;
		xloc = integral_pt()->knot(ipt,1); // arguments are (i,j) where i is the index of the integration poitn and j is the jth coordinate
		
		// interpolate the external data (i.e. the three values of ubar) quadratically
		Vector<double> basis_fn(3);
		psi_1D[0] = xloc*(xloc-1.0)/2.0;	//zero at xloc = 0,1, one at xloc=-1
		psi_1D[1] = (1.0-xloc)*(1.0+xloc);	//zero at xloc = -1,1, one at xloc=0
		psi_1D[2] = xloc*(xloc+1.0)/2.0;	//zero at xloc = -1,0, one at xloc=1
	}

	// helper-function for computing off-diagonal jacobian blocks analytically
	void fill_in_off_diagonal_jacobian_blocks_analytic(Vector<double>& residuals, DenseMatrix<double>& jacobian)	{
	
		//solution indices
		unsigned c_nodal_ind = this->field_index_strat_helmholtz();
		unsigned ubar_external_data_ind = 0;	//the index of the (ONLY) external data for this element

		//storage for local unknown and local equation
		int local_unknown = 0, local_eqn = 0;

		//number of nodes in this element
		unsigned  n_node = this->nnode();

		//memory allocation for shape functions and derivatives
		Shape psi(n_node);
		DShape dpsi(n_node,2);

		//number of integration points
		const unsigned n_intpt = this->integral_pt()->nweight();
		
		//storage for nodal_local_coords
		Vector<double> ext_local_dof_coords(2);

		//loop over integration points
		for(unsigned ipt=0; ipt<n_intpt; ipt++)	{
			
			//weight associated with this integration point
			double w = this->integral_pt()->weight(ipt);
			//derivatives of shape and test functions
			double J = this->dshape_eulerian_at_knot(ipt, psi, dpsi);
			// pre-multiply weights and jacobian
			double W = w*J;

			//compute oneD shape functions
			Vector<double> psi_1D;
			oneD_shape_at_knot_helm(ipt, psi_1D);
				
			//Assemble jacobian terms
			//loop over test functions
			for(unsigned j=0; j<n_node; j++)	{
				// calculate derivatives of the c-equation residuals w.r.t to ubar
				local_eqn = this->nodal_local_eqn(j, c_nodal_ind);
				//if not a dirichlet boundary
				if(local_eqn>=0)	{
					//number of entries in the (ONLY) external data 
					unsigned nexternal_dofs = external_data_pt(0)->nvalue(); 
					// loop over shape functions
					for(unsigned j2=0; j2<nexternal_dofs; j2++)
					{
						local_unknown = this->
							external_local_eqn(ubar_external_data_ind, j2);
						//if not a dirichlet boundary
						if(local_unknown>=0)	{
							jacobian(local_eqn, local_unknown) += 
								-alpha*psi[j]*psi_1D[j2]*W;
						}
					}
				}
			}// end of outer loop over nodes
		}//end of loop over integration points
	}//end of fill_in_off_diagonal_....

public:
	//Constructor: call underlying constructors
	HelmholtzDensityElement1D()
		:HELM_ELEMENT(), ElementWithExternalElement()	{
		this->set_ninteraction(1);
	}
	
	//overload the get_source... function from the stratified helmholtz element which provides the coupling
	// between the helmholtz and moment equations
	void get_source_strat_helmholtz(const unsigned& ipt, const Vector<double>& x, double* source) const	{
		// interpolate the external data (i.e. the three values of ubar) using 1D shape functions
		Vector<double> psi_1D;
		oneD_shape_at_knot_helm(ipt, psi_1D);
		unsigned n_shape = psi_1D.size();
		double ubar=0.0;
		for(unsigned j=0; j<n_shape; j++)	{
			ubar += psi_1D[j]*external_data_pt(0)->value(j);
		}
		*source = alpha*ubar;
	}

	// Calculate the element's contribution to the residual vector. We do this by
	// overloading the fill_in_... function for the residuals
	void fill_in_contribution_to_residuals(Vector<double>& residuals)	{
		// fill in the residuals of the helmholtz equation
		HELM_ELEMENT::fill_in_contribution_to_residuals(residuals);
	}

	// Calculate the Jacobian by finite differencing, we overload the fill_in_... function
	void fill_in_contribution_to_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian)	{
		//finite diffeerences for all elements
		//ElementWithExternalElement::fill_in_contribution_to_jacobian(residuals, jacobian); 
		//analytical jacobian	
		//diagonal blocks
		HELM_ELEMENT::fill_in_contribution_to_jacobian(residuals, jacobian);
		//off-diagonal blocks
		this->fill_in_off_diagonal_jacobian_blocks_analytic(residuals, jacobian);

	}
};//end of HelmholtzDensityElement1D

//=============================== Start of problem class for combined problem===================
// Combined Density, bulk AI problem
// /============================================================================================
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
class FullDensityProblem : public Problem
{
public:
	// Constructor
	FullDensityProblem(const unsigned& num_el_u_fine, const unsigned& num_el_u_coarse, const unsigned& num_el_x); 

	//Destructor (empty)
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
	void actions_before_distribute() {}
	//actions taken after the problem has been distributed over multiple processors
	void actions_after_distribute()	{
		oomph_info << "Problem successfully distributed. In actions_after_distribute(...)\n";
		Multi_domain_functions::setup_multi_domain_interactions 
			<DENSITY_ELEMENT, HELM_ELEMENT>(this, density_msh_pt(), helm_msh_pt());
		oomph_info << "Successfully re-setup multi-domain interactions\n";
	}
	//create elements and mesh responsible for enforcing population constraint
	void create_mass_constraint_elements();

	//create the elements and mesh responsible for the nonlocal ubar auxillary equation
	void create_ubar_elements();

	//computes the population density at a given value of x, uses simpsons rule to evaluate the 
	// integral in the u-direction and then 1d quadratic interpolation(same as shape functions) to interpolate in x
	double compute_population_density(double x);

	//assigns initial conditions 
	void set_initial_condition();

	// document the solution
	void doc_solution();

	//write a restart file
	void doc_restart();

	//dump problem to disk to allow for restart
	void dump_it(ofstream& dump_file);

	//setup restart from restart file specified via command line
	void restart();

	// read solution from disk for restart
	void restart(ifstream& restart_file);

	// function that outputs important problem parameters
	void output_parameters_to_file(ofstream& param_file);
	
	//access function for mesh used by the density equation
	RectangularQuadMesh<DENSITY_ELEMENT>* density_msh_pt()	{
		return dynamic_cast<RectangularQuadMesh<DENSITY_ELEMENT>*>(Density_msh_pt);
	}

	//access function for mesh used by the helmholtz equaiton
	RectangularQuadMesh<HELM_ELEMENT>* helm_msh_pt()	{
		return dynamic_cast<RectangularQuadMesh<HELM_ELEMENT>*>(Helm_msh_pt);
	}

	//Data item for traded density
	Data* Traded_density_data_pt; 

	//mesh for computing the mass in the density and enforcing the mass constraint
	Mesh* Mass_constraint_msh_pt;

	//access function to solution file counter
	unsigned doc_number() { return unsigned(Doc_info.number());  }
	
protected:
	//Mesh for density elements
	RectangularQuadMesh<DENSITY_ELEMENT>* Density_msh_pt;
	//mesh for stratified helmholtz elements
	RectangularQuadMesh<HELM_ELEMENT>* Helm_msh_pt;

	//mesh for the set of ubar elements
	Mesh* Ubar_msh_pt;	
private:
	//DocInfo objects for outputing solution and restart files to disk
	DocInfo Doc_info;
	DocInfo Restart_info;
}; //end of problem class

///================ Start of  Problem constructor==============================================
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::FullDensityProblem(
		const unsigned& num_el_u_fine, const unsigned& num_el_u_coarse, const unsigned& num_el_x) 
{
	//setup labels for output
	Doc_info.set_directory(CommandLineArgs::Argv[2]);
	Restart_info.set_directory(CommandLineArgs::Argv[2]); //use same directory for restart and output files
	//initialise file counters
	Doc_info.number() = 0;	
	Restart_info.number() = 0;

	// allocate timestepper (steady = no timestepping)
	add_time_stepper_pt(new Steady<2>);
	//add_time_stepper_pt(new BDF<2>(true));
	Newton_solver_tolerance = 1e-8;	
	Max_newton_iterations = 15;	// maximum allowable newton iterations for non-linear solver before giving up	
	Max_residuals = 1e6;	//maximum allowed residuals before giving up, set high because even very rough initial guesses can converge

	//total (fixed) population in the domain. This is the total mass in the solution which we will prescribe through the use of MassConstraintElements
	total_population = rho_target*Lx;

	//Create the mesh objects for each equation
	//parameter to the mesh. The arguments to the constructors indicates
	//the number of elements in each direction in each mesh as well as the size of the domain.
	u_min_mesh = max(u_star - Lu/2.0,0.0);
	//u_min_mesh = 0.0;
	Density_msh_pt = new 
		RectangularQuadMesh<DENSITY_ELEMENT>(num_el_u_fine, num_el_x, u_min_mesh, u_min_mesh + Lu, 0.0, Lx, time_stepper_pt());
	Helm_msh_pt = new
		RectangularQuadMesh<HELM_ELEMENT>(num_el_u_coarse, num_el_x, u_min_mesh, u_min_mesh + Lu, 0.0, Lx, time_stepper_pt());
	oomph_info << "Created sub meshes for density and helmholtz elements.\n";
	//store copies of all elements on each processor as 'halo elements' 
	Density_msh_pt->set_keep_all_elements_as_halos();
	Helm_msh_pt->set_keep_all_elements_as_halos();

	//hijack the first data value in the first density element (any element would work) so that the residual can
	//be replaced with a mass-enforcing constraint
	Traded_density_data_pt = dynamic_cast<DENSITY_ELEMENT*>(
			Density_msh_pt->element_pt(0))->hijack_nodal_value(0,0); //pointer to traded data value	
	//create the mesh of elements required for imposing the mass constraint
	create_mass_constraint_elements();
	oomph_info << "Created sub mesh for enforcing population constraint.\n";
	// Loop over the elements in the density mesh to setup element specific things that cannot be done by element constructor
	unsigned n_element_density = density_msh_pt()->nelement();	
	for(unsigned i=0; i<n_element_density;i++)	{
		//upcast from generalised element to present type
		DENSITY_ELEMENT *el_pt = dynamic_cast<DENSITY_ELEMENT*>(density_msh_pt()->element_pt(i));
	
		//set the effective (spatial) diffusion coefficient and diffusion-in-u coefficient in the density eqution
		el_pt->epsilon_diffusion_pt() = &epsilon_reg;
		el_pt->motility_fct_pt() = get_D1;

		//ignore external geometric data in 'external' helmholtz element when computing jacobian matrix
		//because the interaction does not involve derivatives
		el_pt->ignore_external_geometric_data();

	}
	//create sub mesh of ubar elements 
	create_ubar_elements();
	oomph_info <<"Created submesh for auxillary field ubar(x).\n";

	// Loop over the elements in the helmholtz mesh to setup element specific things that cannot be done by mesh contrustor
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

	// combine the submeshes and build the global mesh
	add_sub_mesh(Density_msh_pt);
	add_sub_mesh(Helm_msh_pt);
	add_sub_mesh(Mass_constraint_msh_pt);
	add_sub_mesh(Ubar_msh_pt);
	build_global_mesh();	//this needs to be done before the ubar elements are added to the mesh, otherwise they get ignored
	oomph_info << "Global Mesh built.\n";

	//setup the interactions between the density and helmholtz submeshes
	Multi_domain_functions::setup_multi_domain_interactions<DENSITY_ELEMENT,HELM_ELEMENT>(this, density_msh_pt(), helm_msh_pt());			
	oomph_info << "Interactions between submeshes setup\n";

	//Assign the global and local equations numbers for the problem 
	oomph_info << "Number of equations is " << assign_eqn_numbers(false) << std::endl;	//must use false for MPI, switch to true for serial
	
}  // end of constructor
 
//=================Start of destructor=====================================
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::~FullDensityProblem()	{
	//delete the meshes and global data
	delete Traded_density_data_pt;
	delete Density_msh_pt;
	delete Helm_msh_pt;
	delete Mass_constraint_msh_pt;
	delete Ubar_msh_pt;
}

//==================Start of create_mass_constraint_elements()==========================
// creates the elements and meshes required to enforce the population constraint
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::create_mass_constraint_elements()	{
	//`master element' that stores the traded data and prescribed population
	MassConstraintElement* mass_constraint_element = 
		new MassConstraintElement(&total_population, Traded_density_data_pt, 0);
	//add this element to a new mesh that enforces the population constraint
	Mass_constraint_msh_pt = new Mesh;
	Mass_constraint_msh_pt->add_element_pt(mass_constraint_element);
	
	//loop over all elements in the density mesh and create a mass_computation element for each density element
	unsigned n_element = density_msh_pt()->nelement();
	for(unsigned j=0; j<n_element; j++)	{
		//get the current density element 
		DENSITY_ELEMENT* density_el_pt = 
			dynamic_cast<DENSITY_ELEMENT*>(density_msh_pt()->element_pt(j));
		//create the new mass computation element and pass 
		MassComputationElement<DENSITY_ELEMENT>* el_pt = 
			new MassComputationElement<DENSITY_ELEMENT>(density_el_pt);
		//specify the mass constraint element as the `master element'
		el_pt->set_mass_constraint_element(mass_constraint_element);
		//add the new element to the mesh
		Mass_constraint_msh_pt->add_element_pt(el_pt);
	}
	Mass_constraint_msh_pt->set_keep_all_elements_as_halos();	
}

//==================Start of create_ubar_elements()==========================
// creates the elements and meshes required for the nonlocal ubar auxillary equation
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::create_ubar_elements()	{
	Ubar_msh_pt = new Mesh;	//create mesh
	for(unsigned i=0; i<num_el_x; i++)	{	//add the elements
		UbarElement<DENSITY_ELEMENT>* ubar_el_pt =
			new UbarElement<DENSITY_ELEMENT>(density_msh_pt(), i);
		Ubar_msh_pt->add_element_pt(ubar_el_pt);
	}
	Ubar_msh_pt->set_keep_all_elements_as_halos();	//retain halo copies on all procs
}

//==========================Start of restart()=====================================
//opens restart file and calls the restart routine
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::restart()	{
	//Pointer to restart file
	ifstream* restart_file_pt=0;
	//get restart filename
	char filename[100];
	sprintf(filename, "%s_on_proc%i.dat", GlobalParams::Restart_filename.c_str(),
			this->communicator_pt()->my_rank());
	//Open restart file
	restart_file_pt = new ifstream(filename, ios_base::in);
	if(restart_file_pt !=0)	{
		oomph_info << "Opened " << filename << 
			"for restart. " << std::endl;
	}
	else	{
		std::ostringstream error_stream;
		error_stream 
			<< "ERROR while trying to open " << filename <<
			" for restart " << std::endl;
	}
	//Read restart data:
	if(restart_file_pt !=0)	{
		restart(*restart_file_pt);
	}
}

//==================Start of compute_population_density()==========================
//calculates the mass density of the density solution and interpolates it at the specified value of x
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
double FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::compute_population_density(double x)	{

	//figure out which elements we need for computing the mass density for the given x
	double dx = Lx/num_el_x;
	unsigned row_ind = x/dx;
	if(row_ind == num_el_x)	{ // only happens when x==Lx. just use the top row of elements
		row_ind--;
	}
	double mass_density = 0.0;
	//number of elements in the u-direction in mass constraint mesh mesh
	unsigned n_element_u = Density_msh_pt->nx();
	//loop over masscomputationelements, note that the first element in this mesh is a mass CONSTRAINT element, so skip this one
	for(unsigned j=0; j<n_element_u; j++)	{
		unsigned elem_ind = row_ind*n_element_u + j + 1;	//+1 skips the mss constraint element
		//pointer to jth element in mass compuatation, upcast from generalised element to present type
		MassComputationElement<DENSITY_ELEMENT>* mass_el_pt 
			= dynamic_cast<MassComputationElement<DENSITY_ELEMENT>*>(Mass_constraint_msh_pt->element_pt(elem_ind));
		mass_density += mass_el_pt->interpolate_elemental_mass_density(x);
	}	
	return mass_density;
}

// read solution from disk for restart
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::restart(ifstream& restart_file)	{
	
	//get the predicted next dt (not needed but read anyway) and solution indices for documentation from the restart file
	unsigned local_doc_info_number=0, local_restart_info_number=0;

	if(restart_file.is_open())	{
		oomph_info << "restart file exists" << endl;

		//read in the predicted and next value of dt
		string input_string;
		getline(restart_file, input_string, '#');
		restart_file.ignore(80,'\n');  //ignore rest of line
		//Read in the index of the last output file from the restart file and have the soln output files indexed from here
		getline(restart_file, input_string, '#');		//get up to terminating comment
		restart_file.ignore(80,'\n');	//ignore rest of line
		local_doc_info_number = unsigned(atof(input_string.c_str())); //read in step number
		//Read in the index of the restart file and index future outputs from there
		getline(restart_file, input_string, '#');	//get up to terminating comment
		restart_file.ignore(80,'\n');	//ignore rest of line
		local_restart_info_number = unsigned(atof(input_string.c_str()));
	}
	else	{
		oomph_info << "restart file does not exist" << endl;
	}
	//determine the doc_info number for documenting solution at future times
	unsigned doc_info_number = 0;
	MPI_Allreduce(&local_doc_info_number, &doc_info_number, 1, MPI_UNSIGNED, MPI_MAX,
			communicator_pt()->mpi_comm());
	Doc_info.number() = doc_info_number;
	oomph_info << "Restarted Doc_info.number(): " << doc_info_number << endl;

	//now determine the doc_info number of the restart file
	unsigned restart_info_number = 0;
	MPI_Allreduce(&local_restart_info_number, &restart_info_number, 1, MPI_UNSIGNED, MPI_MAX,
			communicator_pt()->mpi_comm());
	Restart_info.number() = restart_info_number + 1; //needs to be incremented so that the next restart file doesn't overwrite the current one
	oomph_info << "Restarted Restart_info.number(): " << restart_info_number << endl;

	Problem::read(restart_file);
}

//==================Start of set_initial_conditions()==========================
// sets the initila condtion
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::set_initial_condition()	{
	   
	// Choose initial timestep
     	double dt0=1e-3;

	// Initialise timestep -- also sets the weights for all timesteppers
	// in the problem.
	initialise_dt(dt0);
	
	// Backup time in global Time object
	double backed_up_time=time_pt()->time();
		   
	// Past history fo needs to be established for t=time0-deltat, ...
	// Then provide current values (at t=time0) which will also form
	// the initial guess for the first solve at t=time0+deltat
		   
	// Vector of solution values
	Vector<double> soln(3);
	Vector<double> x(2);
		   
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
		for (unsigned i=0;i<num_nod_density;i++)	{
			// Get nodal coordinates
			x[0]=density_msh_pt()->node_pt(i)->x(0);
			x[1]=density_msh_pt()->node_pt(i)->x(1);

			// Get intial solution
			GlobalParams::get_IC(time,x,soln); 

			// Assign initial condition for n
			density_msh_pt()->node_pt(i)->set_value(itime,0,soln[0]);

			// Loop over coordinate directions: Previous position = present position
       		      for (unsigned j=0;j<2;j++)	{
      			      density_msh_pt()->node_pt(i)->x(itime,j)=x[j];
		      }
	     }
	     //Loop over the nodes in the helmholtz-moment mesh to set initial guess everywhere
	     for(unsigned i=0;i<num_nod_helm;i++)	{
		     //get nodal coordinates
		     x[0]=helm_msh_pt()->node_pt(i)->x(0);
		     x[1]=helm_msh_pt()->node_pt(i)->x(1);

		     // Get initila condition
		     GlobalParams::get_IC(time,x,soln);

		     //Assign initial condition for c and ubar
		     helm_msh_pt()->node_pt(i)->set_value(itime, 0, soln[1]); // c is at index 0 at each node

		     //no dirichlet boundary conditions for ubar or c 
		     //Loop over coordinate directiosn: previous position = present position
		     for (unsigned j=0;j<2;j++)	{
      			      helm_msh_pt()->node_pt(i)->x(itime,j)=x[j];
      		      }
	     }
	     unsigned n_ubar_data = Ubar_msh_pt->nelement();
    	     //assign values in the ubar mesh
	     double dx = Lx/num_el_x;
    	     for(unsigned i=0; i<n_ubar_data; i++)	{
    		     for(unsigned j=0; j<3; j++)	{// loop over the three values stored in each ubar element
     			     x[1] = i*dx + j*dx/2.0;	//value of x[0] doesn't matter here
     			     GlobalParams::get_IC(time,x,soln);
			     //upcast from generalised element to present type
			     UbarElement<DENSITY_ELEMENT>* ubar_el_pt = 
				     dynamic_cast<UbarElement<DENSITY_ELEMENT>*>(Ubar_msh_pt->element_pt(i));
			     ubar_el_pt->ubar_data_pt()->set_value(j,soln[2]);
    		     }
	     }
	} //end of loop over previous times   
     	// Reset backed up time for global timestepper
	time_pt()->time()=backed_up_time;
}// end of set_initial_condition()	


//======================================Start of doc_info =============================================
//document the solution
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::doc_solution()	
{
	ofstream some_file;
	char filename[100];	
	oomph_info << std::endl;
	oomph_info << "=================================================" << std::endl;
	oomph_info << "Docing solution" << std::endl;
	oomph_info << "=================================================" << std::endl;

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

	//incremenet the doc counter
	Doc_info.number()++;	
}// end of doc_solution

//======================================Start of doc_restart =============================================
//write restart file
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::doc_restart() {

	ofstream some_file;
	char filename[100];
	oomph_info << std::endl;
	oomph_info << "=================================================" << std::endl;
        oomph_info << "Writing restart file." <<  std::endl;
        oomph_info << "=================================================" << std::endl;

	//output restart data to file
	sprintf(filename, "%s/restart%i_on_proc%i.dat", Restart_info.directory().c_str(), 
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
	//write dummy next time step in case we want to run a time-dependent simulation starting from last output
	dump_file << "0.0001" << " # suggested next dt (dummy)" << endl;
	//write the index of the documented solutions for future restart
	dump_file << Doc_info.number() << " # output soln index" << endl;
	//write the index of the restart file for future restarts
	dump_file << Restart_info.number() << " # current restart index " << endl;
	
	// call generic dump()
	Problem::dump(dump_file);
}

// function that outputs important problem parameters
template<class DENSITY_ELEMENT, class HELM_ELEMENT>
void FullDensityProblem<DENSITY_ELEMENT, HELM_ELEMENT>::output_parameters_to_file(ofstream& param_file)	{
	param_file << "a = "<< a <<"\nL = " << L <<"\nK = " << K <<"\nlambda = " << lambda
		<< "\nbeta = " << beta_decay << "\nalpha = " << alpha << "\nD2 = " << D2 <<"\nD1_star = " 
		<< D1_star << "\nD1_p_star = " << D1_p_star << "\nDinf_0 = " << D_inf_0 << "\nm = " << rho_target << "\nepsilon = " 
		<< epsilon_reg << endl;
}

//===========================End of problem class definitions ==============================

// Check that input is valid and extract parameters from file
int read_input_params_from_file()
{
	ifstream infile; //input filestream
	if(CommandLineArgs::Argc<=2)	{// invalid number of inputs
		std::ostringstream error_stream;
		error_stream 
			<< "ERROR, too few command line arguments, usage is ./unstead_density_md ";
		error_stream << "param_filename result_directory restart_filename" << std::endl;
	     
	     throw OomphLibError(
	      error_stream.str(),
	      OOMPH_CURRENT_FUNCTION,
	      OOMPH_EXCEPTION_LOCATION);
	}
	infile.open(CommandLineArgs::Argv[1]);

	// read in parameters from file and assign them appropriately
	std::string buffer;
	// number of parameters that need to be assigned and number assigned so far
	unsigned n_params = 18, num_params_assigned=0; 
	while(std::getline(infile,buffer))	{	//read in one line at a time
		std::string param_name, param_value;
		//ignore lines that don't contain `=`
		size_t pos = buffer.find('=');
		if(pos != std::string::npos)	{ // found '=' in line
			param_name = buffer.substr(0,pos); // left-hand side 
			param_value = buffer.substr(pos + 1, std::string::npos); //right-hand side
			double val = atof(param_value.c_str()); // convert RHS to double from string

			//remove whitespaces from parameter name
			param_name.erase(remove(param_name.begin(), param_name.end(),' '), param_name.end());
			if(param_name=="a")		{a = val; 	 	}
			else if(param_name=="L")	{L = val;  		}
			else if(param_name=="K")	{K = val;  		}
			else if(param_name=="lambda")	{lambda = val;  	}
			else if(param_name=="D2")	{D2 = val;  		}
			else if(param_name=="beta")	{beta_decay = val;  	}
			else if(param_name=="alpha")	{alpha = val;  		}
			else if(param_name=="m")	{rho_target = val;  	}
			else if(param_name=="D1_star")	{D1_star = val;		}
			else if(param_name=="D1_p_star"){D1_p_star = val;	}
			else if(param_name=="D_inf_0")	{D_inf_0 = val;		}
			else if(param_name=="epsilon")	{epsilon_reg = val;  	}
			else if(param_name=="Lu")	{Lu = val;		}
			else if(param_name=="Lx")	{Lx = val;		}
			else if(param_name=="Nu_density"){num_el_density_u = val;}
			else if(param_name=="Nu_helm")	{num_el_helm_u = val;	}
			else if(param_name=="Nx")	{num_el_x = val;	}
			else if(param_name=="rho_eval")	{rho_eval = val;		}
			else	{
				oomph_info << "Error, could not assign parameter '" << param_name <<"'\n";
				return 0;// should never get here unless there is a problem assigning a parameter
			}
			//print out parameters that got assigned and their values
			oomph_info << param_name << " = " << val << endl;
			num_params_assigned ++;
		}
	
	} // end of file reached
	if(num_params_assigned != n_params) { //not enough parameters got assigned, perhaps some missing?
		oomph_info << "Error, some parameters left unspecified in input file '" << CommandLineArgs::Argv[1] <<"'\n";
		oomph_info << "Number of parameters set:" << num_params_assigned << "\nNumber of parameters needed:" << n_params << endl;
		return 0;
	}
	return 1;  

}

// computes the intracellular concentration at which the kinetics is in equilibrium with spatially homogeneous AI c=ce. i.e. f(u_star,ce)=0
void compute_ustar_homogeneous()
{
	double ce_quad_B = K - rho_target*alpha*(a + L)/(beta_decay*lambda);
	double ce_quad_C = -a*alpha*rho_target*K/(beta_decay*lambda);
	double ce = 0.5*(-ce_quad_B + std::sqrt(ce_quad_B*ce_quad_B - 4.0*ce_quad_C)); //homogeneous equilibrium value of AI
	u_star = (a + L*ce/(K+ce))/lambda; 
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//======start_of_main==================================================
/// Driver for Problem
//=====================================================================


int main(int argc, char* argv[])
 {
	#ifdef OOMPH_HAS_MPI
	 MPI_Helpers::init(argc, argv);
	 unsigned nproc = MPI_Helpers::communicator_pt()->nproc();
	 cout << "Code has been compiled with mpi support \n"
		 << "Running on " << nproc << " processors." << endl;
	#else
	 cout << "Code compiled WITHOUT mpi support." << endl;
	#endif

	 //Store command line arguments	
	CommandLineArgs::setup(argc, argv);

	//Define possible command line arguments and parse the ones that were specified 

	//name of restart file
	CommandLineArgs::specify_command_line_flag("--restart-file",&GlobalParams::Restart_filename);

	//parse the input
	CommandLineArgs::parse_and_assign();

	//switch  off output modified
	oomph_info.output_modifier_pt() = &default_output_modifier;

	//Define processor-labelled output file for all on-screen otuput
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
 	oomph_info << "Critical intracellular concentration in the homogeneous state u_*=" << u_star << endl;	 

	//Build the problem, pass the density and helmholtz element types to problem
	//the 'Hijacked' element allows us to `hijack' an element and replace a residual with a mass-enforcing constraint
  	FullDensityProblem<
		Hijacked<DensityHelmholtzElement1D<	
					DensityElement1D,
	       				StratifiedHelmholtzElement
					>>,
		HelmholtzDensityElement1D<
					StratifiedHelmholtzElement,
					DensityElement1D
					>
		> 
		problem(num_el_density_u,num_el_helm_u, num_el_x);

	// setup labels for output
	DocInfo doc_info;
	doc_info.set_directory(CommandLineArgs::Argv[2]);
	doc_info.number()=0;	
	
	// Set IC
	problem.set_initial_condition();
	
#ifdef OOMPH_HAS_MPI
	// ====================== distribute the problem ==================
	Vector<unsigned> used_element_partition;
	bool report_stats = false;

	//partition the problem using default oomph-lib partitioning
	Vector<unsigned> used_partition;
	used_partition.resize(problem.mesh_pt()->nelement());
	used_element_partition = problem.distribute(report_stats);

	oomph_info << "Finished distribution of the problem\n";
	// =================================================================
#endif

	//restart form file if restart flag has been specified 
	if(CommandLineArgs::command_line_flag_has_been_set("--restart-file"))	{
		problem.restart();
		oomph_info << "Finished restart\n";
	}
	else	{ // if not restarting, then doc the initial condition
		problem.doc_solution();
	}
		 
	// write initial problem parameters to file
	ofstream param_file;
	char param_filename[100];
	sprintf(param_filename, "%s/params.dat",doc_info.directory().c_str());
	param_file.open(param_filename);
	problem.output_parameters_to_file(param_file);	//write equation-specific quantities
	//now write numerical quantities
	param_file << "\nnum_el_density_u = " << num_el_density_u  
		<<"\n num_el_helm_u = " << num_el_helm_u 
		<< "\nnum_el_x = " << num_el_x << endl;
	param_file.close();
	
	//create output file for bifurcation diagram
	ofstream bifurcation_file;
	char bif_filename[100];
	sprintf(bif_filename,"%s/bifurcation.dat",doc_info.directory().c_str());
#ifdef OOMPH_HAS_MPI
	if(MPI_Helpers::communicator_pt()->my_rank()==0)	
#endif
	{
		bifurcation_file.open(bif_filename, std::ios_base::app);
	}

	//perform one steady-solve to get close to steady-state for current D1p_star
	problem.steady_newton_solve();
	//document bifurcation data
	double m0 = problem.compute_population_density(0.0);
	double m1 = problem.compute_population_density(rho_eval);
#ifdef OOMPH_HAS_MPI
	if(MPI_Helpers::communicator_pt()->my_rank()==0)	
#endif
	{
		bifurcation_file << std::setprecision(10) << D1_p_star << " " << std::setprecision(10) << m0 << " " << m1 << endl;
	}

	//============Construct bifurcation diagram using pseudo-arclength continuation================
	//start continuation run
	double ds_max = 5e-4;	//maximum allowed arc-length step (approximate)
	double ds = 5e-3;	//initial arc-length step (approximate), sign determines whether we continue left or right
	for(unsigned i=0; i<500; i++)	{
		oomph_info << "Current D1_p_star = " << std::setprecision(10) << D1_p_star << endl;
		double sign = 0.0;
		if(ds>0)	{
			sign = 1.0;
		}
		else	{
			sign = -1.0;
		}
		//set next arc length step less than or equal to maximal step length, accounting for sign
		ds = sign*min(ds_max,abs(ds));
		oomph_info << "next arclength step = " << ds << endl;
		//do a continuation step
		ds = problem.arc_length_step_solve(&D1_p_star,ds);	

		//compute bifurcation information (max and min steady-state density always at the edges)
		m0 = problem.compute_population_density(0.0);
		m1 = problem.compute_population_density(rho_eval);
	#ifdef OOMPH_HAS_MPI
		if(MPI_Helpers::communicator_pt()->my_rank()==0)	
	#endif
		{
			bifurcation_file << std::setprecision(10) << D1_p_star << " " << std::setprecision(10) << m0 << " " << m1 << endl;
		}

		//document solution 
		if(i%10==0)	{
			problem.doc_solution();
			problem.doc_restart();
		}		
	}
	       
#ifdef OOMPH_HAS_MPI
	MPI_Helpers::finalize();	
	if(MPI_Helpers::communicator_pt()->my_rank()==0)
#endif
	{
		bifurcation_file.close();		
	}

} // end of main
