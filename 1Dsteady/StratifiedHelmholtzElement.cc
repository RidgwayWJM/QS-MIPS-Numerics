/*
//A class that implements a modified helmholtz equation element.
//Importantly, the Laplacian does not contain derivates with respect to the first variable, only
//the last one, which represent the spatial coordinate

// the modified helmholtz equation is
// c_t = D2 c_xx - beta*c + f(u,x,t),
// where beta is the decay coefficient, D2 the diffusion coeff, and f is a source function

*/

#include "StratifiedHelmholtzElement.h"

// function that evaluates residuals and, if flag==1, the jacobian as well if flag==0, only computes residuals
void StratifiedHelmholtzElement::fill_in_generic_residual_contribution_strat_helmholtz(Vector<double>& residuals, DenseMatrix<double>& jacobian, const unsigned& flag)
  {
	// number of nodes in the element
        const unsigned n_node = nnode();
        //allocate memeory for shape functions and derivatives
        Shape psi(n_node);
        DShape dpsi(n_node,2); // contains derivatives w.r.t (u,x1,...) (u is the first variable)

        // memory for local coords
        Vector<double> s(2);

        const unsigned n_intpt = integral_pt()->nweight();

        // integers to store local equation numbers and unknowns
        int local_eqn_number = 0;
        int local_dof_number = 0;

        // loop over integration points
        for(unsigned ipt=0; ipt<n_intpt; ipt++)
        {

                s[0] = integral_pt()->knot(ipt,0); // arguments are (i,j) where i is the index of the integration poitn and j is the jth coordinate
                s[1] = integral_pt()->knot(ipt,1);


                //get the integral weight
                double w = integral_pt()->weight(ipt);

                //call the derivatives of the shape and test
                double J = dshape_eulerian(s, psi, dpsi);
                double W = w*J;

                // Allocate memory for the solution, their derivatives, and global
                // poisition at the knot point
                double interpolated_x = 0.0, interpolated_u = 0.0, interpolated_c = 0.0,
		       interpolated_dcdx = 0.0, interpolated_dcdt=0.0;

                // Calculate the interpolated values by looping over shape functions
                for(unsigned j=0; j<n_node; j++)
                {
                        interpolated_u += nodal_position(j,0)*psi[j];
                        interpolated_x += nodal_position(j,1)*psi[j];
                        interpolated_c += field_value_strat_helmholtz(j)*psi[j];
                        interpolated_dcdx += field_value_strat_helmholtz(j)*dpsi(j,1);
			interpolated_dcdt += dcdt_strat_helmholtz(j)*psi[j];
                }
		
		//evaluate source function		
		double source = 0.0;
		Vector<double> interpolated_xVec(2);
		interpolated_xVec[0] = interpolated_u;
		interpolated_xVec[1] = interpolated_x;
		get_source_strat_helmholtz(ipt, interpolated_xVec, &source);
	
		//ASSEMBLE RESIDUALS AND JACOBIAN
                //loop over test functions (same as shape functions here)
                for (unsigned j=0; j<n_node; j++)
                {
                        // get the local equation number
                        local_eqn_number = nodal_local_eqn(j, field_index_strat_helmholtz());
                        // if not a Dirichlet boundary, then add contribution
                        if(local_eqn_number >=0)
                        {
                                // Steady-State Density Equation Residuals (Weak Form) //
                                // diffusion term
                                residuals[local_eqn_number] += diff_coeff()*interpolated_dcdx*dpsi(j,1)*W;
                                // decay term and source function
                                residuals[local_eqn_number] += (decay_coeff()*interpolated_c - source)*psi[j]*W;
				
				//time derivative term
				residuals[local_eqn_number] += interpolated_dcdt*psi[j]*W;
				
                                //calculate Jacobian if the flag is set
                                if(flag)
                                {
                                        // loop over shape functions
                                        for(unsigned i=0; i<n_node; i++)
                                        {
                                                //get local degree of freedom number
                                                local_dof_number = nodal_local_eqn(i, field_index_strat_helmholtz());
                                                // add contribution if not on a Dirichlet boundary
                                                if(local_dof_number>=0)
                                                {
                                                        // diffusion term
                                                        jacobian(local_eqn_number, local_dof_number) +=
                                                                diff_coeff()*dpsi(i,1)*dpsi(j,1)*W;
                                                        // kinetics term
                                                        jacobian(local_eqn_number, local_dof_number) +=
                                                                decay_coeff()*psi[i]*psi[j]*W;
							//time-derivative term
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
  void StratifiedHelmholtzElement::dinterpolated_c_ddata_strat_helmholtz(Vector<double>& s, 
		  Vector<double>& dc_ddata, Vector<unsigned>& global_eqn_numbers)	{
  
  	const unsigned n_node = nnode();	  
	const unsigned c_nodal_ind = field_index_strat_helmholtz();
	unsigned n_dof = 0;

	Shape psi(n_node);
	shape(s, psi);

	// find number of dofs associated with interpolated c
	for(unsigned j=0; j<n_node; j++)	{
		int global_eqn = this->node_pt(j)->eqn_number(c_nodal_ind);
		if(global_eqn>=0)	{ // if not a diirchlet boundary point, add to the count
			++n_dof;
		}
	}
	dc_ddata.resize(n_dof,0.0);
	global_eqn_numbers.resize(n_dof,0);

	//loop over nodes again and set eqn numbers and derivatives
	unsigned count = 0;
	for(unsigned j=0;j<n_node;j++)	{
		int global_eqn = this->node_pt(j)->eqn_number(c_nodal_ind);
		if(global_eqn>=0)	{ // if not a diirchlet boundary point, add to the count
			global_eqn_numbers[count] = global_eqn;
			dc_ddata[count] = psi[j];
			++count;
		}

	}
  }

  
	 
  void StratifiedHelmholtzElement::output(ostream &outfile)
  {
        unsigned n_node = nnode();
        // print out global coordinate for each node and value in the field variables
        // at each node
        for(unsigned i=0;i<n_node;i++)
        {
			Vector<double> s(2);
			local_coordinate_of_node(i,s);
			unsigned external_data_component_ind = int(s[1] + 1.5);
                	outfile << nodal_position(i,0) << " " <<  nodal_position(i,1) << " " 
				<< this->field_value_strat_helmholtz(i) << " " 
				<< external_data_pt(0)->value(external_data_component_ind) << std::endl;
        }
  }
 
  unsigned StratifiedHelmholtzElement::self_test()
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
  
