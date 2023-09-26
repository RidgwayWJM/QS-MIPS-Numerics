/* A class that implements the density equation in 2 spatial dimensions for our chemically stratified model. The equation is
 * n_t = epsilon*n_uu + D1(u)\nabla^2 n  -[n*f(u,x)]_u,
 * where f is a user-specified kinetic function 
 */


#include "generic/Qelements.h"
#include "math.h"

using namespace oomph;
using namespace std;

class DensityElement2D : public virtual QElement<3,3>
{
  public:
  //Constructor
  //initialises Source function pointer and pointer to decay coeff to null
  DensityElement2D() : QElement<3,3>() 
  {
	  Motility_fct_pt = 0; // pointer to diffusion coefficient
	  Kinetic_fct_pt = 0; //pointer to motility function
	  Epsilon_diffusion_pt = 0;
  }
  
  //funciton pointer to kinetic function
  typedef void (*DensityKineticFctPt) (const Vector<double>& x, double* f);
  
  // Access function to kinetic function
  DensityKineticFctPt& kinetic_fct_pt()
  {
	  return Kinetic_fct_pt;
  }
  // Access function to kinetic function (const version)
  DensityKineticFctPt kinetic_fct_pt() const
  {
	  return Kinetic_fct_pt;
  }

  //function pointer to motility function
  typedef void (*DensityMotilityFctPt) (const Vector<double>& x, double* D1);

  // Access function to motility function
  DensityMotilityFctPt& motility_fct_pt()
  {
	  return Motility_fct_pt;
  }

  //Access function to motility function (const version)
  DensityMotilityFctPt motility_fct_pt() const
  {
	  return Motility_fct_pt;
  }

  /// Get pointer to diffusion-in-u coefficient
  double*& epsilon_diffusion_pt()
  {
    return Epsilon_diffusion_pt;
  }
  
  
  /// Get the diffusion-in-u coefficient
  double epsilon_diffusion()
  {
    return *Epsilon_diffusion_pt;
  }


  
  // amount of memory needed for allocation at each node
  inline unsigned required_nvalue(const unsigned &j) const {return 1;}

  //index at which the solution is stored at each node
  virtual unsigned field_index_density() const {return 0;}

  //Access function to the first (and in single-physics problems only) data value stored at each node 
  double field_value_density(const unsigned &j) {return node_pt(j)->value(field_index_density());}


  //interpolates the field variable c at local coordinate s from local nodal values
  inline double interpolated_n_density(const Vector<double>& s) const
  {
	  // number of nodes in element
	  const unsigned n_node = nnode();

	  // local shape function
	  Shape psi(n_node);

	  // find values of shape function
	  shape(s, psi);

	  //initialise value of c
	  double interpolated_n = 0.0;

	  //loop over element nodes and sum
	  for(unsigned j=0; j<n_node; j++)
	  {
		  interpolated_n += nodal_value(j, field_index_density())*psi[j];
	  }
	  return interpolated_n;

  }

  double dndt_density(const unsigned& j) const
  {
	  TimeStepper* time_stepper_pt = node_pt(j)->time_stepper_pt();

	  //initialise dcdt
	  double dndt = 0.0;

	  //Loop over the timesteps if the timestepper is non-steady
	  if(!time_stepper_pt->is_steady())
	  {
		//index at which the field variable is stored at local node j
		const unsigned density_nodal_index = field_index_density();
		//number of timesteps (past and present)
		const unsigned n_time = time_stepper_pt->ntstorage();

		//add contributions to the time derivative
		for(unsigned t=0; t<n_time; t++)
		{
			dndt += time_stepper_pt->weight(1,t)*nodal_value(t, j, density_nodal_index);
		}
	  }
	  return dndt;
  }

  //evaluates the user-defined kinetic function (if present) at global coordinate x, returns zero if not specified
  inline virtual void get_kinetics_density(const unsigned& ipt, const Vector<double>& x, double* kinetics) const
  {
	  // return zero if no kinetic function has been specified
	  if(Kinetic_fct_pt == 0)
	  {
		  *kinetics = 0.0;
	  }
	  else
	  {
		  (*Kinetic_fct_pt)(x,kinetics);
	  }
  }

  //evaluates the user-defined motility coefficient (if present) at global coordinate x, returns 1 if not specified
  inline virtual void get_motility_coeff(const unsigned& ipt, const Vector<double>& x, double* D1) const
  {
	  //return 1 if no motility coefficient has been specified
	  if(Motility_fct_pt == 0)
	  {
		  *D1 = 1.0;
	  }
	  else
	  {
		  (*Motility_fct_pt)(x,D1);
	  }
  }

  // wrapper for function that evaluates residuals
  void fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
	   fill_in_generic_residual_contribution_density(residuals, GeneralisedElement::Dummy_matrix, 0);
  }
  
  // wrapper for function that evaluates jacobian
  void fill_in_contribution_to_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian)
  {
	   fill_in_generic_residual_contribution_density(residuals, jacobian, 1);
  }

  void get_global_eqn_numbers(Vector<unsigned>& global_eqn_numbers);

  //self-test: return 0 for ok
  unsigned self_test();

  void output(ostream& outfile);

 protected:
 
  virtual void fill_in_generic_residual_contribution_density(Vector<double>& residuals, DenseMatrix<double>& jacobian, const unsigned& flag);

  //pointer to kinetic function
  DensityKineticFctPt Kinetic_fct_pt;
  //pointer to motility function
  DensityMotilityFctPt Motility_fct_pt;
  //pointer to diffusion-in-u coefficient 
  double* Epsilon_diffusion_pt;
  //spatial dimension of the problem
  unsigned Dim = 2;

}; 
