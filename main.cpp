#include "main.h"

const char * usage="Usage: %s [-eorbILifl] benchmark\n\
    -i input file\n\
    -f output file\n\
    -l log file (default to screen)\n"
;

const char * usage2="Usage: %s -i input -f output\n";

int main(int argc, char * argv[]){
//#if 0
	int c;
	//char * logfile="/dev/null";
	char * logfile=NULL;
	char * input=NULL, * output=NULL;
	bool input_flag = false, output_flag = false;
//#if 0
	while( ( c = getopt(argc, argv, "i:f:e:o:r:b:l:LI")) != -1 ){
		switch(c){
		case 'l':
			logfile = optarg;
			break;
		case 'i':
			input = optarg;
			input_flag = true;
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
	if( !input_flag || ! output_flag ){
		fprintf(stderr, usage2, argv[0]);
		exit(EXIT_FAILURE);
	}

	open_logfile(logfile);
//#endif
	//output = "./out";
	if( freopen(output, "w", stdout) == NULL )
		report_exit("Ouptut file error\n");

	Tran tran;
	vector<CKT_TOP> cktlist; 
	// start to parfile
	// vector<Circuit *> cktlist;
	Parser parser(&cktlist);	
	parser.parse(input, tran);
	return 0;
	// do the job
	size_t i=0;
//#pragma omp parallel for private(i)	
	for(i=0;i<cktlist.size();i++){
		Circuit * ckt1 = cktlist[i].ckt1;
		Circuit * ckt2 = cktlist[i].ckt2;
		ckt1->print_netlist();
		clog<<endl;
		ckt2->print_netlist();
		// functions for transient solve
		//ckt->solve(tran);

		//ckt->solve_DC();	
		//delete cktlist[i];
	}
	cktlist.clear();
	//tran.print_tr_nodes();

	close_logfile();

//#endif
	return 0;
}
