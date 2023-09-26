/* A class that implements the structured density equation in 2 spatial dimensions. The equation is
 * n_t = epsilon*n_uu + D1(u)\nabla^2 n  -[n*f(u,x)]_u,
 * where f is a user-specified kinetic function 
 */

#include "DensityElement2D.h"



//compute element residuals and/or jacobian 
// flag=1 compute both
// flag=0 compute residuals only
// NOTE: fill_in_... functions should not initialise residuals or jacobians to zero as otherwise in multi-physics problems it will wipe contributions from
// other elements
void DensityElement2D::fill_in_generic_residual_contribution_density(Vector<double>& residuals, DenseMatrix<double>& jacobian, const unsigned& flag)
{
	// number of nodes in the element
	const unsigned n_node = nnode();
	//allocate memeory for shape functions and derivatives
	Shape psi(n_node);
	DShape dpsi(n_node, Dim + 1); // contains derivatives both w.r.t x and u (u is the first variable)
	
	// memory for local coords
	Vector<double> s(Dim + 1);
			
	const unsigned n_intpt = integral_pt()->nweight();
	// integers to store local equation numbers and unknowns
	int local_eqn_number = 0;
	int local_dof_number = 0;

	// loop over integration points
	for(unsigned ipt=0; ipt<n_intpt; ipt++)	{
		
		for(unsigned k=0; k<Dim+1; k++)	{
			s[k] = integral_pt()->knot(ipt,k); // arguments are (i,j) where i is the index of the integration poitn and j is the jth coordinate
		}
		//get the integral weight
		double w = integral_pt()->weight(ipt);
		//call the derivatives of the shape and test
		double J = dshape_eulerian(s, psi, dpsi);
		double W = w*J;

		// Allocate memory for the solution, their derivatives, and global
		// poisition at the knot point
		Vector<double> interpolated_x(Dim + 1,0.0);
		Vector<double> interpolated_dndx(Dim + 1, 0.0);
		double interpolated_n = 0.0, interpolated_dndt=0.0;

		// Calculate the interpolated values by looping over shape functions
		for(unsigned j=0; j<n_node; j++)
		{
			//calculate interpolated coordinates
			for(unsigned k=0; k<Dim + 1; k++)	{
				interpolated_x[k] += nodal_position(j,k)*psi[j];
				interpolated_dndx[k] += field_value_density(j)*dpsi(j,k);
			}
			//interpolated field and time-derivative
			interpolated_n += field_value_density(j)*psi[j];
			interpolated_dndt += dndt_density(j)*psi[j];
		}
		// evaluate kinetics and motility coefficient
		double kinetics = 0.0;
		double D1 = 0.0;

		get_kinetics_density(ipt, interpolated_x, &kinetics);
		get_motility_coeff(ipt, interpolated_x, &D1);
					
		//ASSEMBLE RESIDUALS AND JACOBIAN
		//loop over test functions (same as shape functions here)
		for (unsigned j=0; j<n_node; j++)
		{
			// get the local equation number
			local_eqn_number = nodal_local_eqn(j, field_index_density());

			// if not a Dirichlet boundary, then add contribution
			if(local_eqn_number >=0)
			{
				// Steady-State Density Equation Residuals (Weak Form) //
				// diffusion terms
				for(unsigned k=0; k<Dim; k++)	{
					residuals[local_eqn_number] += D1*interpolated_dndx[k+1]*dpsi(j,k+1)*W;
				}
				//diffusion in epsilon
				residuals[local_eqn_number] += epsilon_diffusion()*interpolated_dndx[0]*dpsi(j,0)*W;

				// nonlinear kinetics term
				residuals[local_eqn_number] += -kinetics*interpolated_n*dpsi(j,0)*W;	
				//time-derivative term
				residuals[local_eqn_number] += interpolated_dndt*psi[j]*W;
				//calculate Jacobian if the flag is set
				if(flag)
				{
					// loop over shape functions
					for(unsigned i=0; i<n_node; i++)
					{
						//get local degree of freedom number
						local_dof_number = nodal_local_eqn(i, field_index_density());
						// add contribution if not on a Dirichlet boundary
						if(local_dof_number>=0)
						{
							// diffusion terms
							for(unsigned k=0; k<Dim; k++)	{
								jacobian(local_eqn_number, local_dof_number) += D1*dpsi(i,k+1)*dpsi(j,k+1)*W;
							}
							//diffusion in epsilon
							jacobian(local_eqn_number, local_dof_number) += epsilon_diffusion()*dpsi(i,0)*dpsi(j,0)*W;
							// kinetics term
							jacobian(local_eqn_number, local_dof_number) +=
								-kinetics*psi[i]*dpsi(j,0)*W;
							//time-derivative (mass-matrix) term	
							jacobian(local_eqn_number, local_dof_number) += 
								psi[i]*psi[j]*node_pt(i)->time_stepper_pt()->weight(1,0)*W; //zeroth weight of first time-derivative
						}
					}//end of inner loop over test functions
				}
			}
		}//end of outer loop over test functions
	}//end of outer loop over integration points

}//end of fill_in_generic_residual_contribution(...)

  //determines the global equation numbers for the current element
  void DensityElement2D::get_global_eqn_numbers(Vector<unsigned>& global_eqn_numbers)	{
  
  	const unsigned n_node = nnode();	  
	const unsigned nodal_ind = field_index_density();
	unsigned n_dof = 0;
	// find number of dofs associated with interpolated n
	for(unsigned j=0; j<n_node; j++)	{
		int global_eqn = this->node_pt(j)->eqn_number(nodal_ind);
		if(global_eqn>=0)	{ // if not a diirchlet boundary point, add to the count
			n_dof++;
		}
	}
	global_eqn_numbers.resize(n_dof,0);

	//loop over nodes again and set eqn numbers
	unsigned count = 0;
	for(unsigned j=0;j<n_node;j++)	{
		int global_eqn = this->node_pt(j)->eqn_number(nodal_ind);
		if(global_eqn>=0)	{ // if not a diirchlet boundary point, add to the count
			global_eqn_numbers[count] = global_eqn;
			count++;
		}

	}
  }

  void DensityElement2D::output(ostream &outfile)	{
	unsigned n_node = nnode();
	// print out global coordinate for each node and value in the field variable
	// at each node
	for(unsigned i=0;i<n_node;i++)	{
		for(unsigned k=0; k<Dim+1; k++)	{
			outfile << nodal_position(i,k) << " "; 
		}
		outfile  << field_value_density(i) << std::endl;
	}
}
 
  unsigned DensityElement2D::self_test()
  {
	bool passed = true;
	if(FiniteElement::self_test()!=0)
	{
		passed = false;
	}
	if(passed)
	{
		return 0;
	}
	else
	{
		return 1;
	}
   }

