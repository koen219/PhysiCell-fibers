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
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 

    Cell_Definition* pCD = find_cell_definition( "ECM"); 
    pCD->functions.custom_cell_rule = fiber_custom_function; 
//    pCD->functions.update_migration_bias = rotating_migration_bias; 
    pCD->functions.contact_function = fiber_contact_function; 

    pCD->phenotype.mechanics.maximum_number_of_attachments = 2; 
    pCD->phenotype.mechanics.attachment_elastic_constant = 0.2; // 0.2; 
    pCD->phenotype.mechanics.attachment_rate = 1; // .1; 
    pCD->phenotype.mechanics.detachment_rate = 0; 

    pCD = find_cell_definition( "pusher"); 
    if( pCD )
    {
        pCD->functions.update_migration_bias = rotating_migration_bias;
        pCD->functions.update_phenotype = pusher_phenotype_function; 
    }

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
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
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
    
    std::vector<std::string> out = paint_by_number_cell_coloring(pCell); 

    if( pCell->type_name == "pusher" )
    {
        out[0] = "orange"; 
        out[2] = "orange";
        out[3] = "orange";
        return out; 
    }
    
    int n_attached= pCell->state.number_of_attached_cells(); 

    if( n_attached == 0 )
    { out[0] = "grey"; }

    if( n_attached == 1 )
    { out[0] = "red"; }

    if( n_attached == 2 )
    { out[0] = "blue"; }

    if( n_attached == 3 )
    { out[0] = "aquamarine"; }

    if( n_attached > 3 )
    { out[0] = "black"; }

    if( get_single_signal( pCell, "custom:crosslink") > 0.5 )
    {
        for( int n=0; n < pCell->state.number_of_attached_cells() ; n++ )
        {
            Cell* pTest = pCell->state.attached_cells[n]; 
            if( get_single_signal(pTest,"custom:crosslink" ) > 0.5 )
            {
                std::cout << "Yell something obscene!" << std::endl; 
                out[0] = "magenta"; 
            }

        }


    }


    out[2] = out[0]; 
    out[3] = out[0]; 

    return out; 
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

/*
#include <cmath>
struct Point {
    double x;
    double y;
    double z;
    double dot(Point P){
        return x*P.x + y*P.y + z*P.z;
    }
    double mag(){
        return sqrt(this.dot(this))
    }
}
Point cross(Point A, Point B){ // TODO
    return Point(0.0, 0.0, 0.0) ;
}

Point getAngularHarmonicForce_Monasse(Point A, Point B, Point C, double k, double theta0){
    Point BA = B - A;
    Point BC = B - C;
    Point CB = C - B;

    double BA.mag(); 
    double BC_mag = BC.mag();

//  If BA and BC are colinear, take an arbitrary orthogonal vector
    // https://math.stackexchange.com/q/3077100
    Point BA_cross_BC = np.cross( BA, BC );
    if ( not (BA_cross_BC > 0)){
        Point orthogonal_vector = {
            BA.y + BA.z,
            BA.z - BA.x,
            -BA.x - BA.y
        };
        Point p_a = cross( BA , orthogonal_vector );
        p_a = p_a / p_a.mag();
        Point p_c = cross( CB , orthogonal_vector );
        p_c = p_c / p_c.mag();
    }
    else{
        Point p_a = cross( BA , BA_cross_BC );
        p_a = p_a / p_a.mag();
        Point p_c = cross( CB , BA_cross_BC );
        p_c = p_c / p_c.mag();
    }

    // rounding is necessary to avoid numerical errors 
    // when vectors are colinear the argument can be = 1.0000000000000002
    // which is outside of the range of np.arccos, defined for [-1, 1]
    double costheta = BA.dot(BC) / BA_mag * BC_mag;
    double theta = arccos(round(costheta))

    double delta_theta = k*(theta - theta0)

    double F_a = delta_theta/BA_mag * p_a 
    double F_c = delta_theta/BC_mag * p_c 

    double F_b = -F_a - F_c

    return Point(F_a, F_b, F_c)

}
*/

/*
 * import numpy as np

def getAngularHarmonicForce_Monasse(A, B, C, k, theta0):
    """
    Monasse, Bernard, and Frédéric Boussinot. 
    "Determination of forces from a potential in molecular dynamics." 
    arXiv preprint arXiv:1401.1181 (2014).
    Returns the angular harmonic forces (F_a, F_b, F_c) in the bond
    A--B--C
    - A, B, and C are lists with coordinates of the points.
    - k is the spring constant.
    - theta0 the rest length of the bond in radians.
    The forces on A, B, and C are calculated as:
    
    F_a = k*(theta - theta0) * p_a # was: F_a = ( -k*(theta - theta0) / |BA| ) * p_a
    
    F_c = k*(theta - theta0) * p_c # was: F_c = ( -k*(theta - theta0) / |BC| ) * p_c 
    
    F_b = -F_a - F_c
    
    - p_a is a unit vector in the ABC plane orthogonal to (BA)
    - p_c is a unit vector in the ABC plane orthogonal to (CB)
    - theta is the angle ABC in radians
    """
    A = np.array(A)
    B = np.array(B)
    C = np.array(C)

    BA = B - A
    BC = B - C
    CB = C - B

    BA_mag = mag(BA)   # |BA|
    BC_mag = mag(BC)   # |BC|

    # If BA and BC are colinear, take an arbitrary orthogonal vector
    # https://math.stackexchange.com/q/3077100
    BA_cross_BC = np.cross( BA, BC )
    if mag(BA_cross_BC) == 0:
        orthogonal_vector = np.array( [BA[1] + BA[2], BA[2] - BA[0], - BA[0] - BA[1]] )
        p_a = np.cross( BA , orthogonal_vector )
        p_a = p_a / mag(p_a)
        p_c = np.cross( CB , orthogonal_vector )
        p_c = p_c / mag(p_c)
    else:
        p_a = np.cross( BA , BA_cross_BC )
        p_a = p_a / mag(p_a)
        p_c = np.cross( CB , BA_cross_BC )
        p_c = p_c / mag(p_c)

    # rounding is necessary to avoid numerical errors 
    # when vectors are colinear the argument can be = 1.0000000000000002
    # which is outside of the range of np.arccos, defined for [-1, 1]
    theta = np.arccos( round( BA.dot(BC) / ( BA_mag * BC_mag ), 10) )

    delta_theta = k*(theta - theta0)

    F_a = delta_theta/BA_mag * p_a 
    F_c = delta_theta/BC_mag * p_c 

    F_b = -F_a - F_c

    return(F_a, F_b, F_c)
 */

// std::vector<Point> getAngularHarmonicForce_Monasse(Point A, Point B, Point C, double k, double theta0); 

/* Force calculation based on:
	Monasse, Bernard, and Frédéric Boussinot. 
    "Determination of forces from a potential in molecular dynamics." 
    arXiv preprint arXiv:1401.1181 (2014).


    Returns the angular harmonic forces (F_a, F_b, F_c) in the bond
    A--B--C
    - A, B, and C are lists with coordinates of the points.
    - k is the spring constant.
    - theta0 the rest length of the bond in radians.
    The forces on A, B, and C are calculated as:
    
    F_a = k*(theta - theta0) * p_a
    
    F_c = k*(theta - theta0) * p_c
    
    F_b = -F_a - F_c
    
    - p_a is a unit vector in the ABC plane orthogonal to (BA)
    - p_c is a unit vector in the ABC plane orthogonal to (CB)
    - theta is the angle ABC in radians
*/
std::vector< std::vector<double> > compute_angular_force_contributions( Cell* pCell , Phenotype& phenotype , double dt )
{
    if( pCell->state.number_of_attached_cells() < 2 )
    { std::cout << "bad move genius" << std::endl; exit(-1); }

    double theta0 = get_single_signal( pCell , "custom:preferred_angle" ); 
    double k = get_single_signal( pCell , "custom:angular_spring_constant" ); 
    
    Cell* pLeft = pCell->state.attached_cells[0];
    Cell* pMiddle = pCell; 
    Cell* pRight = pCell->state.attached_cells[1]; 

    std::vector<double> A = pLeft->position; 
    std::vector<double> B = pMiddle->position;  
    std::vector<double> C = pRight->position;  

//    std::cout << A << " " << B << " " << C << std::endl; 
 //   std::cout << theta0 << " " << k << std::endl; 

    std::vector<double> BA = B-A; 
    std::vector<double> BC = B-C;
    std::vector<double> CB = C-B; 

    double BA_mag = norm(BA); 
    double BC_mag = norm(BC);

//    std::cout << BA << " " << BC << " " << CB <<  " " << BA_mag << " " << BC_mag << std::endl; 

    std::vector<double> BA_cross_BC = BioFVM::cross_product( BA, BC );

//    std::cout << BA_cross_BC << std::endl; 

    std::vector<double> p_a, p_b, p_c;

	// If BA and BC are colinear (cross product is 0 vector), 
	// take an arbitrary orthogonal vector (see: https://math.stackexchange.com/q/3077100)
    if ( not (norm(BA_cross_BC) > 0)){
        std::vector<double> orthogonal_vector = {
            BA[1] + BA[2],
            BA[2] - BA[0],
            -BA[0] - BA[1]
            };
        p_a = BioFVM::cross_product( BA , orthogonal_vector );
        p_a = (1/ norm(p_a)) * p_a ;
        p_c = BioFVM::cross_product( CB , orthogonal_vector );
        p_c = (1/ norm(p_c))* p_c  ;
    }
    else{
        p_a = BioFVM::cross_product( BA , BA_cross_BC );
        p_a = (1/norm(p_a)) * p_a  ;
        p_c = BioFVM::cross_product( CB , BA_cross_BC );
        p_c = (1 / norm(p_c) ) * p_c ;
    }

    // Capping the value is necessary to avoid numerical errors 
    // when vectors are colinear the argument can be = 1.0000000000000002
    // which is outside of the range of acos, defined for [-1, 1]
    double costheta = BioFVM::dot_product(BA,BC) / ( BA_mag * BC_mag );
    if( costheta > 1 )
    { costheta = 1; }
    if( costheta < -1 )
    { costheta = -1; }

    double theta = acos(costheta); 

    double delta_theta = k*(theta - theta0);
 //   std::cout << p_a << " " << p_c << std::endl; 
 //   std::cout << costheta << " " << theta << " " << delta_theta << std::endl; 


    std::vector<double> F_a = (delta_theta/BA_mag) * p_a  ;
    std::vector<double> F_c = (delta_theta/BC_mag) * p_c ;
    std::vector<double> F_b = (-1.0)*F_a - F_c;
    std::vector< std::vector<double> > out;
    out.push_back(F_a);
    out.push_back(F_b);
    out.push_back(F_c);

    return out;
}



void dynamic_spring_attachments( Cell* pCell , Phenotype& phenotype, double dt )
{
    // check for detachments 
    double detachment_probability = phenotype.mechanics.detachment_rate * dt; 
    for( int j=0; j < pCell->state.attached_cells.size(); j++ )
    {
        Cell* pTest = pCell->state.attached_cells[j]; 
        if( UniformRandom() <= detachment_probability )
        { detach_cells( pCell , pTest ); }
    }

    // check if I have max number of attachments 
    if( pCell->state.attached_cells.size() >= phenotype.mechanics.maximum_number_of_attachments )
    { return; }

    // check for new attachments; 
    double attachment_probability = phenotype.mechanics.attachment_rate * dt; 
    bool done = false; 
    int j = 0; 
    while( done == false && j < pCell->state.neighbors.size() )
    {
        Cell* pTest = pCell->state.neighbors[j]; 
        if( pTest->state.number_of_attached_cells() < pTest->phenotype.mechanics.maximum_number_of_attachments )
        {
            // double affinity = phenotype.mechanics.cell_adhesion_affinity[ pTest->type]; 
            std::string search_string = "adhesive affinity to " + pTest->type_name; 
            double affinity = get_single_behavior( pCell , search_string );

            double prob = attachment_probability * affinity; 
            if( UniformRandom() <= prob )
            {
                // attempt the attachment. testing for prior connection is already automated 
                attach_cells( pCell, pTest ); 
                if( pCell->state.attached_cells.size() >= phenotype.mechanics.maximum_number_of_attachments )
                { done = true; }
            }
        }
        j++; 
    }
    return; 
}

void dynamic_ECM_link_attachments( Cell* pCell , Phenotype& phenotype, double dt )
{
    // check for detachments 
    double detachment_probability = phenotype.mechanics.detachment_rate * dt; 
    for( int j=0; j < pCell->state.attached_cells.size(); j++ )
    {
        Cell* pTest = pCell->state.attached_cells[j]; 

        bool detachment_possible = true; 
        if( get_single_signal(pCell,"custom:crosslink") > 0 || 
            get_single_signal(pTest,"custom:crosslink") > 0 )
        { detachment_possible = false; }

        if( detachment_possible )
        {
            if( UniformRandom() <= detachment_probability )
            { detach_cells( pCell , pTest ); }
        }
    }

    // check if I have max number of attachments 
    if( pCell->state.attached_cells.size() >= phenotype.mechanics.maximum_number_of_attachments )
    { return; }

    // check for new attachments; 
    double attachment_probability = phenotype.mechanics.attachment_rate * dt; 
    bool done = false; 
    int j = 0; 
    while( done == false && j < pCell->state.neighbors.size() )
    {
        Cell* pTest = pCell->state.neighbors[j]; 
        if( pTest->state.number_of_attached_cells() < pTest->phenotype.mechanics.maximum_number_of_attachments )
        {
            // double affinity = phenotype.mechanics.cell_adhesion_affinity[ pTest->type]; 
            std::string search_string = "adhesive affinity to " + pTest->type_name; 
            double affinity = get_single_behavior( pCell , search_string );

            double prob = attachment_probability * affinity; 
            if( UniformRandom() <= prob )
            {
                // attempt the attachment. testing for prior connection is already automated 
                attach_cells( pCell, pTest ); 
                if( pCell->state.attached_cells.size() >= phenotype.mechanics.maximum_number_of_attachments )
                { done = true; }
            }
        }
        j++; 
    }
    return; 
}


void dynamic_ECM_crosslink_attachments( Cell* pCell , Phenotype& phenotype, double dt )
{
    if( pCell->state.number_of_attached_cells() < 1 )
    { return; }

    bool i_am_crosslink = (bool) get_single_signal( pCell , "custom:crosslink" ); 

    // check for crosslink detachments 
    double detachment_probability  = 0; 
    for( int j=0; j < pCell->state.attached_cells.size(); j++ )
    {
        Cell* pTest = pCell->state.attached_cells[j]; 
        bool test_is_crosslink = (bool) get_single_signal( pTest , "custom:crosslink" ); 
        if( i_am_crosslink )
        { detachment_probability = get_single_signal(pCell,"custom:crosslink_detachment_rate") * dt; }
        if( test_is_crosslink )
        { detachment_probability = get_single_signal(pTest,"custom:crosslink_detachment_rate") * dt; }

        if( UniformRandom() <= detachment_probability )
        {
            detach_cells( pCell , pTest ); 
            set_single_behavior( pCell , "custom:crosslink" , 0 ); 
            set_single_behavior( pTest , "custom:crosslink" , 0 ); 
        }
    }

/*
    // check if I have max number of attachments 
    if( pCell->state.attached_cells.size() >= phenotype.mechanics.maximum_number_of_attachments )
    { return; }
*/ 
    // allow me to form a cross link if I have one regular attachment 
    if( pCell->state.number_of_attached_cells() != 1 )
    { return; }

    // check for new attachments; 
    double attachment_probability = get_single_signal( pCell ,"custom:crosslink_attachment_rate" ) * dt; 
    bool done = false; 
    int j = 0; 
    while( done == false && j < pCell->state.neighbors.size() )
    {
        Cell* pTest = pCell->state.neighbors[j]; 
        bool crosslink_possible = true; 
        if( get_single_signal(pTest,"custom:crosslink") > 0.5 )
        { crosslink_possible = false; }
        if( pTest->state.number_of_attached_cells() != 2 ) 
        { crosslink_possible = false; }

        //    std::string search_string = "adhesive affinity to " + pTest->type_name; 
        //    double affinity = get_single_behavior( pCell , search_string );
        // double prob = attachment_probability * affinity; 

        if( crosslink_possible )
        {
            if( UniformRandom() <= attachment_probability )
            {
                // attempt the attachment. testing for prior connection is already automated 
                attach_cells( pCell, pTest ); 
                if( pCell->state.attached_cells.size() >= phenotype.mechanics.maximum_number_of_attachments )
                { done = true; set_single_behavior(pCell,"custom:crosslink" , 1); }
            }
        }
        j++; 
    }
    return; 
}

void fiber_custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
    dynamic_ECM_link_attachments(pCell, phenotype, dt ); 
    dynamic_ECM_crosslink_attachments(pCell, phenotype, dt );

    if( pCell->state.number_of_attached_cells() >= 2 )
    {
        std::vector< std::vector<double> > forces = compute_angular_force_contributions(pCell,phenotype,dt); 
        Cell* pLeft = pCell->state.attached_cells[0];
        Cell* pMiddle = pCell;
        Cell* pRight = pCell->state.attached_cells[1];

        pMiddle->velocity += forces[1]; 
        #pragma omp critical 
        {
            pLeft->velocity += forces[0]; 
            pRight->velocity += forces[2]; 
        }
    }

    if( get_single_signal( pCell, "time") > 500 )
    {
        set_single_behavior( pCell , "quorum factor secretion" , 10); 
        set_single_behavior( pCell , "migration speed" , 1); 
        set_single_behavior( pCell , "chemotactic sensitivity to quorum factor" , -1 ); 
        return; 
    }

    if( PhysiCell_globals.current_time > 200 )
    {
        set_single_behavior( pCell , "migration speed" , 0); 
        return; 
    }

    // just test code from here on down. 


    if( pCell->state.attached_cells.size() == 0 )
    { 
        set_single_behavior( pCell , "chemotactic sensitivity to quorum factor" , 1 );  
        set_single_behavior( pCell , "migration speed" , 2 ); 

        set_single_behavior( pCell , "quorum factor secretion" , 5 ); // 1

        return; 
    }

    if( pCell->state.attached_cells.size() == 1 )
    {
        set_single_behavior( pCell , "chemotactic sensitivity to quorum factor" , -1 ); 

        set_single_behavior( pCell , "quorum factor secretion" , 100 ); 
        set_single_behavior( pCell , "migration bias" , 1 ); 

        phenotype.mechanics.attachment_rate = 0; 

        Cell* pOther = pCell->state.attached_cells[0]; 
        int n_other_attached = pOther->state.number_of_attached_cells(); 

        // if I'm in a chain of at least 3
        if( n_other_attached == 2 ) // && PhysiCell_globals.current_time > 10 )
        {
            phenotype.motility.migration_speed = 4; 
        }
        else
        {
            if( pCell < pOther )
            {
                phenotype.mechanics.attachment_rate = 0; 
                set_single_behavior( pCell , "chemotactic sensitivity to quorum factor" , 1 ); 
                set_single_behavior( pCell , "migration speed" , 2 ); 

            }
            else
            {
                phenotype.mechanics.attachment_rate = .1; 
                set_single_behavior( pCell , "chemotactic sensitivity to quorum factor" , 1 ); 
                set_single_behavior( pCell , "migration speed" , 2 ); 
            }
        }


        return; 
    }

    if( pCell->state.attached_cells.size() == 2 )
    {
        set_single_behavior( pCell , "quorum factor secretion" , 0 ); 
        set_single_behavior( pCell , "migration speed" , 0 ); 
    }

    return; 
}

void fiber_contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{
    // spring link 

    standard_elastic_contact_function(pMe,phenoMe,pOther,phenoOther,dt); 

    // angular spring calculations 
    
    return; 
} 

void rotating_migration_bias( Cell* pCell, Phenotype& phenotype , double dt )
{
    double t =  PhysiCell_globals.current_time; 
    if( t > 200 && t < 400 )
    {
        set_single_behavior( pCell , "migration speed" , 0 ); 
        set_single_behavior( pCell , "quorum factor secretion" , 10 ); 
    }
    else
    {
        double sp = get_single_base_behavior( pCell, "migration speed"); 
        set_single_behavior( pCell , "migration speed" , sp ); 
        set_single_behavior( pCell , "quorum factor secretion" , 0 );         
    }

    double x = pCell->position[0]; 
    double y = pCell->position[1]; 
    double radial = 0.0;
    phenotype.motility.migration_bias_direction = { -y + radial*x , x + radial*y , 0} ; 
    normalize ( &(phenotype.motility.migration_bias_direction) );  

    return; 
}

void pusher_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ 
    if( get_single_signal( pCell , "time") > 400 )
    { set_single_behavior( pCell ,"apoptosis" , 9e99 ); }
    
    return; 
}
