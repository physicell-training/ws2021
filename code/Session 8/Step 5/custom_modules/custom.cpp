/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed
	SeedRandom( parameters.ints("random_seed") );

	/*
	   Put any modifications to default cell definition here if you
	   want to have "inherited" by other cell types.

	   This is a good place to set default functions.
	*/

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL;
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;
	cell_defaults.functions.custom_cell_rule = NULL;
	cell_defaults.functions.contact_function = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	/*
	   This parses the cell definitions in the XML config file.
	*/

	initialize_cell_definitions_from_pugixml();

	/*
	   Put any modifications to individual cell definitions here.

	   This is a good place to set custom functions.
	*/

	cell_defaults.functions.update_phenotype = phenotype_function;
	cell_defaults.functions.custom_cell_rule = custom_function;
	cell_defaults.functions.contact_function = contact_function;

	Cell_Definition* pBacteria = find_cell_definition( "bacteria" );
	pBacteria->functions.update_phenotype = bacteria_phenotype;

	/*
	   This builds the map of cell definitions and summarizes the setup.
	*/

	build_cell_definitions_maps();
	display_cell_definitions( std::cout );

	return;
}

void setup_microenvironment( void )
{
	// set domain parameters

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.

	// initialize BioFVM

	initialize_microenvironment();

	return;
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Zmin = microenvironment.mesh.bounding_box[2];

	double Xmax = microenvironment.mesh.bounding_box[3];
	double Ymax = microenvironment.mesh.bounding_box[4];
	double Zmax = microenvironment.mesh.bounding_box[5];

	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0;
		Zmax = 0.0;
	}

	double Xrange = Xmax - Xmin;
	double Yrange = Ymax - Ymin;
	double Zrange = Zmax - Zmin;

	double center_x = 0.5*(Xmin+Xmax);
	double center_y = 0.5*(Ymin+Ymax);
	double center_z = 0.5*(Zmin+Zmax);

	// create some of each type of cell

	Cell* pC;

	// find cell definitions
	Cell_Definition* pBacteria =
		find_cell_definition( "bacteria" );

	Cell_Definition* pSupplier =
		find_cell_definition( "supplier" );



	for( int k=0; k<parameters.ints("number_of_bacteria"); k++ )
	{
		std::vector<double> position = {0,0,0};
		double r = NormalRandom(0,1) *
			parameters.doubles( "radius_bacteria_region" );
		double theta = 6.28318530718 * UniformRandom();

		position[0] = center_x + r*cos(theta);
		position[1] = center_y + r*sin(theta);
		position[2] = center_z;

		pC = create_cell( *pBacteria );
		pC->assign_position( position );
	}

	for( int k=0; k<parameters.ints("number_of_suppliers"); k++ )
	{
		std::vector<double> position = {0,0,0};
		position[0] = Xmin + UniformRandom()*Xrange;
		position[1] = Ymin + UniformRandom()*Yrange;
		position[2] = Zmin + UniformRandom()*Zrange;

		pC = create_cell( *pSupplier );
		pC->assign_position( position );

		pC->is_movable = false;
	}

	return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; }

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; }

void bacteria_phenotype( Cell* pCell, Phenotype& phenotype , double dt )
{
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; // don't bother doing this function again!
		return;
	}

	// find my cell definition
	// don't use static if you plan to use this for more than one cell type
	static Cell_Definition* pCD = find_cell_definition( pCell->type_name );

	// find the index of resource
	static int nR = microenvironment.find_density_index( "resource" );

	// find the index of quorum factor
	static int nQ = microenvironment.find_density_index( "quorum factor" );

	// index of necrotic death model
	static int nNecrosis = 1; // PhysiCell_constants::necrosis_death_model;

	// sample microenvironment at cell position to get resource
	double R = pCell->nearest_density_vector()[nR];
	double q = pCell->nearest_density_vector()[nQ];

	// check for necrotic death
	if( R < pCell->custom_data["R_necrosis"] )
	{ phenotype.death.rates[nNecrosis] = pCell->custom_data["necrosis_rate"]; }
	else
	{ phenotype.death.rates[nNecrosis] = 0.0; }

	// set birth rate

	// set proliferation
	// first, set to the cell line rate
	phenotype.cycle.data.transition_rate(0,1) =
		pCD->phenotype.cycle.data.transition_rate(0,1);

	// scale with R
	double scaling_factor = (R - pCell->custom_data["R_necrosis"])
		/ (pCell->custom_data["R_max_growth"] - pCell->custom_data["R_necrosis"]);
	if( scaling_factor > 1 )
	{ scaling_factor = 1.0; }
	if( scaling_factor < 0 )
	{ scaling_factor = 0.0; }

	// scale with quorum factor
	double Qlow = pCell->custom_data["quorum_low_prolif"];
	double Qhigh = pCell->custom_data["quorum_high_prolif"];
	double constant = 4.0 / pow( Qhigh - Qlow , 2.0 );

	// no proliferation if Q is low or high
	if( q < Qlow || q > Qhigh )
	{ scaling_factor = 0.0; }

	scaling_factor *= ( constant*(q-Qlow)*(Qhigh-q) );

	// multiply by scaling factor
	phenotype.cycle.data.transition_rate(0,1) *= scaling_factor;

	// get the cell line's motile speed
	phenotype.motility.migration_speed = pCD->phenotype.motility.migration_speed;

	// get a scaling factor
	scaling_factor = 1.0 - q / pCell->custom_data["quorum_motility_slowdown"];
	if( scaling_factor < 0.0 )
	{ scaling_factor = 0.0; }

	// scale migration speed
	phenotype.motility.migration_speed *= scaling_factor;

	// scale migration bias by the quorum factor
	scaling_factor = q / pCell->custom_data["quorum_motility_slowdown"];
	if( scaling_factor > 1.0 )
	{ scaling_factor = 1.0; }
	phenotype.motility.migration_bias = pCD->phenotype.motility.migration_bias;
	phenotype.motility.migration_bias *= scaling_factor;

	return;
}
