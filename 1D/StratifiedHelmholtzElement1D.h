/*
//A class that implements a modified helmholtz equation element.
//Importantly, the Laplacian does not contain derivates with respect to the first variable, only
//the last one, which represent the spatial coordinate

// the modified helmholtz equation is
// c_t = Dc c_{xx} - beta*c + f(u,x,t),
// where beta is the decay coefficient, Dc the diffusion coeff, and f is a source function

*/

#include "generic/Qelements.h"
#include "math.h"

using namespace oomph;
using namespace std;

class StratifiedHelmholtzElement1D : public virtual QElement<2,3>
{
  public:
  //Constructor
  //initialises Source function pointer and pointer to decay coeff to null
  StratifiedHelmholtzElement1D() : QElement<2,3>() 
  {
	  Source_fct_pt = 0; // pointer to source function
	  Decay_coeff_pt = 0; // pointer to decay coefficient
	  Diff_coeff_pt = 0; // pointer to diffusion coefficient
  }
  
  //funciton pointer to source function
  typedef void (*StratifiedHelmholtzSourceFctPt) (const Vector<double>& x, double* f);
  
  // Access function to source function
  StratifiedHelmholtzSourceFctPt& source_fct_pt()
  {
	  return Source_fct_pt;
  }
  // Access function to source function (const version)
  StratifiedHelmholtzSourceFctPt source_fct_pt() const
  {
	  return Source_fct_pt;
  }

  /// Get pointer to decay coefficient
  double*& decay_coeff_pt()
  {
    return Decay_coeff_pt;
  }
  
  
  /// Get the decay coefficient 
  double decay_coeff()
  {
    return *Decay_coeff_pt;
  }

  double*& diff_coeff_pt()
  {
	  return Diff_coeff_pt;
  }
  double diff_coeff()
  {
	  return *Diff_coeff_pt;
  }
  
  // amount of memory needed for allocation at each node
  inline unsigned required_nvalue(const unsigned &j) const {return 1;}

  //index at which the solution is stored at each node
  virtual unsigned field_index_strat_helmholtz() const {return 0;}

  //Access function to the first (and in single-physics problems only) data value stored at each node 
  double field_value_strat_helmholtz(const unsigned &j) {return node_pt(j)->value(field_index_strat_helmholtz());}


  //interpolates the field variable c at local coordinate s from local nodal values
  inline double interpolated_c_strat_helmholtz(const Vector<double>& s) const
  {
	  // number of nodes in element
	  const unsigned n_node = nnode();

	  // local shape function
	  Shape psi(n_node);

	  // find values of shape function
	  shape(s, psi);

	  //initialise value of c
	  double interpolated_c = 0.0;

	  //loop over element nodes and sum
	  for(unsigned j=0; j<n_node; j++)
	  {
		  interpolated_c += this->nodal_value(j, field_index_strat_helmholtz())*psi[j];
	  }
	  return interpolated_c;

  }

  double dcdt_strat_helmholtz(const unsigned& j) const
  {
	  TimeStepper* time_stepper_pt = this->node_pt(j)->time_stepper_pt();

	  //initialise dcdt
	  double dcdt = 0.0;

	  //Loop over the timesteps if the timestepper is non-steady
	  if(!time_stepper_pt->is_steady())
	  {
		//index at which the field variable is stored at local node j
		const unsigned c_nodal_index = field_index_strat_helmholtz();
		//number of timesteps (past and present)
		const unsigned n_time = time_stepper_pt->ntstorage();

		//add contributions to the time derivative
		for(unsigned t=0; t<n_time; t++)
		{
			dcdt += time_stepper_pt->weight(1,t)*nodal_value(t, j, c_nodal_index);
		}
	  }
	  return dcdt;
  }

  //evaluates the user-defined source function (if present) at global coordinate x, returns zero if not specified
  inline virtual void get_source_strat_helmholtz(const unsigned& ipt, const Vector<double>& x, double* source) const
  {
	  // return zero if no source function has been specified
	  if(Source_fct_pt == 0)
	  {
		  *source = 0.0;
	  }
	  else
	  {
		  (*Source_fct_pt)(x,source);
	  }
  }

  // wrapper for function that evaluates residuals
  void fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
	   fill_in_generic_residual_contribution_strat_helmholtz(residuals, GeneralisedElement::Dummy_matrix, 0);
  }
  
  // wrapper for function that evaluates jacobian
  void fill_in_contribution_to_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian)
  {
	   fill_in_generic_residual_contribution_strat_helmholtz(residuals, jacobian, 1);
  }

  void dinterpolated_c_ddata_strat_helmholtz(Vector<double>& s, 
		  Vector<double>& dc_ddata, Vector<unsigned>& global_eqn_numbers); 

  //self-test: return 0 for ok
  unsigned self_test();

  void output(ostream& outfile);

 protected:
 
  virtual void fill_in_generic_residual_contribution_strat_helmholtz(Vector<double>& residuals, DenseMatrix<double>& jacobian, const unsigned& flag);

  //pointer to source function
  StratifiedHelmholtzSourceFctPt Source_fct_pt;	

  // pointer to decay coefficient
  double* Decay_coeff_pt;

  // pointer to diffusion coefficient
  double* Diff_coeff_pt;
	

}; 
