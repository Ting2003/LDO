#include "main.h"

const char * usage="Usage: %s [-eorbILifl] benchmark\n\
    -i input file\n\
    -f output file\n\
    -l log file (default to screen)\n"
;

const char * usage2="Usage: %s -i input -a input_ldo -f output\n";

int main(int argc, char * argv[]){
//#if 0
	int c;
	//char * logfile="/dev/null";
	char * logfile=NULL;
	char * input=NULL, * output=NULL;
	char * input_ldo = NULL;
	bool input_flag = false, output_flag = false;
	bool input_flag_ldo = false;
//#if 0
	while( ( c = getopt(argc, argv, "i:a:f:e:o:r:b:l:LI")) != -1 ){
		switch(c){
		case 'l':
			logfile = optarg;
			break;
		case 'i':
			input = optarg;
			input_flag = true;
			break;
		case 'a':
			input_ldo = optarg;
			input_flag_ldo = true;
			break;
		case 'f':
			output = optarg;
			output_flag = true;
			break;
		case '?':
		default:
			fprintf(stderr, usage2, argv[0]);
			exit(EXIT_FAILURE);
		}
	}
	//if( argc == optind ) report_exit(usage2);
	if( !input_flag || !input_flag_ldo || ! output_flag ){
		fprintf(stderr, usage2, argv[0]);
		exit(EXIT_FAILURE);
	}

	open_logfile(logfile);
//#endif
	//output = "./out";
	if( freopen(output, "w", stdout) == NULL )
		report_exit("Ouptut file error\n");

	Tran tran;
	// start to parfile
	Chip chip;
	//vector<Circuit *> cktlist;
	Parser parser(&chip);
	// parser.parse_ldo(input_ldo);
	parser.parse(input, tran);

	// do the job
	size_t i=0;
//#pragma omp parallel for private(i)	
	/*for(i=0;i<cktlist.size();i++){
		Circuit * ckt = cktlist[i];
		// functions for transient solve
		ckt->solve(tran);

		//ckt->solve_DC();	
		delete ckt;
	}*/
	//tran.print_tr_nodes();

	close_logfile();

//#endif
	return 0;
}
