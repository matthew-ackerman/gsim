#include "interface.h"
#include "map_file.h"
#include "state.h"
#include "sample_name.h"
#include "correl_data.h"

#include <cstring>
#include <sstream>
#include <tuple>
#include <map>
#include <fstream>
#include <ctime>
#include <omp.h>

#include "Eigen/Core"

#define UINT uint32_t
#define WORD 32

#define BLOCK 256

#define t(a)	a.transpose()
#define Av(a,b) Eigen::VectorXf(a.array() * b.array() )
#define Am(a,b) Eigen::MatrixXf( a.array() * b.array() )
#define center(a,b) Eigen::MatrixXf( a.colwise() - b )

int main (int argc, char **argv){

	std::cerr << __FILE__ << std::endl;


	std::string names_file="", input_file="";
	int indX=-1, indY=-1, model=0;
	std::string namex="", namey="";
	double sa=1, sd=1, rho=0;
	double ma=0, md=0;

	Environment env;
	env.set_name("gwas");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("gwas on blups. Please direct questions to matthew.s.ackerman@gmail.com");

	env.optional_arg('s',"state",  state_name,      "please .", ".");
	env.optional_arg('b',"blups",  blups_name,      "please .", ".");

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);

	//Eigen::initParallel();
	Eigen::setNbThreads(4);
	std::cerr << "using " << Eigen::nbThreads( ) << " threads.\n";

	Flat_file <Relatedness> rel_file;
	Relatedness rel;

	if (rel_name=="")
	{
		rel_file.open(READ);
	}
	else
	{
		rel_file.open(rel_name.c_str(), READ);
	}
	rel=state_file.read_header();
	rel_file.read(Pstates);
	rel_file.close();

	Flat_file <Phenotype> pheno_file;
	Phenotype p;

	if (pheno_name=="")
	{
		pheno_file.open(READ);
	}
	else
	{
		pheno_file.open(pheno_name.c_str(), READ);
	}
	p=pheno_file.read_header();
	pheno_file.read(p);
	pheno_file.close();

	std::cerr << "Sample: " << rel.sample_size() << ", " << " genome: "<<  rel.genome_size() << std::endl;
	
	int N(Pstates.sample_size());
	int LEN(Pstates.genome_size()*32);
	
	S = ( va*MtM-( (va*ma2)/(va+ma2) )*2*Mt1*t(Mt1) ) + ( vd*HtH-( (vd*md2)/(vd+md2) )*2*Ht1*t(Ht1) ) + 4*( (vad+sad*ma*md)/(sqrt(va*vd)+ma*md )*MtH-(sad*ma*md)/(sqrt(va*vd)+ma*md)*( Ht1*t(Mt1) + Mt1*t(Ht1) ) )
	p S
}
