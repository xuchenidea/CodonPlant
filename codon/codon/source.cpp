#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>   
#include <sys/timeb.h>
#include <sys/types.h>
#include "main.h"

#define CODON_ROW 30
#define CODON_COLOMN 6

int SEQ_NUM = 20;					//5-1024
int DUP_LENGTH = 10;				//7-11
int ITERATION_NUM = 10000;			
float MUTION_RATIO = 5.0f;			//1-100

char codon[CODON_ROW][CODON_COLOMN][4] = {
	/*
	{ "TTT", "TTC" },
	{ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" },
	{ "ATT", "ATC", "ATA" },
	{ "ATG" },
	{ "GTT", "GTC", "GTA", "GTG" },
	{ "AGT", "AGC", "TCT", "TCC", "TCA", "TCG" },
	{ "CCT", "CCC", "CCA", "CCG" },
	{ "ACT", "ACC", "ACA", "ACG" },
	{ "GCT", "GCC", "GCA", "GCG" },
	{ "TAT", "TAC" },
	{ "TAA", "TAG","TGA" },
	{ "CAT", "CAC" },
	{ "CAA", "CAG" },
	{ "AAT", "AAC" },
	{ "AAA", "AAG" },
	{ "GAT", "GAC" },
	{ "GAA", "GAG" },
	{ "TGT", "TGC" },
	{"TGG"},
	{ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" },
	{ "GGT", "GGC", "GGA", "GGG" },
	*/
	'\0'
};

float opt_codon[CODON_ROW][CODON_COLOMN][2] = {
/*	{ { 4844 ,1.04 },{ 4508 ,0.96 } },
	{ { 1532, 0.47 },{ 4816, 1.48 },{ 5533, 1.70 },{ 3422, 1.05 },{ 1551, 0.48 },{ 2646 ,0.81 } },
	{ { 5782, 1.46 },{ 4300, 1.08 },{ 1830, 0.46 } },
	//	{	{ 5293, 1.00 }},
	{ { 6607, 1.64 },{ 3272, 0.81 },{ 1759, 0.44 },{ 4482, 1.11 } },
	{ { 2353 ,0.86 },{ 2652, 0.97 },{ 4583, 1.67 },{ 2649, 0.96 },{ 3389, 1.23 },{ 861, 0.31 } },

	{ { 4482, 1.58 },{ 1883, 0.66 },{ 4175, 1.47 },{ 822, 0.29 }, }, //ccg	
	{ { 4577, 1.56 },{ 3297, 1.12 },{ 3142, 1.07 },{ 733, 0.25 } },
	{ { 8246, 1.85 },{ 3640, 0.82 },{ 4894, 1.10 },{ 1057, 0.24 } }, //gcg

	{ { 3532, 1.02 },{ 3394, 0.98 } },
	{ { 291, 1.04 },{ 169, 0.60 },{ 383, 1.36 } },
	{ { 2507, 1.08 },{ 2117, 0.92 } },

	{ { 3642, 0.97 },{ 3846, 1.03 } },
	{ { 4553, 1.03 },{ 311, 0.97 } },

	{ { 5487, 0.66 },{ 11221, 1.34 } },
	{ { 7946, 1.33 },{ 3972, 0.67 } },  //gat gac
	{ { 6161, 0.88 },{ 7856, 1.12 } },

	{ { 1664, 0.95 },{ 1840, 1.05 } },
	//	{	{ 2457, 1.00 }},

	{ { 1966, 1.05 },{ 1288, 0.69 },{ 943, 0.50 },{ 719, 0.38 },{ 2881, 1.54 },{ 3455, 1.84 } },
	{ { 5735, 1.38 },{ 3062, 0.74 },{ 5163, 1.24 },{ 2689, 0.65 } },//ggg*/
};

// special codon should not exceed 99
char condon_special[100][7] = { "AATAAA", "ATAAT", "AATTAA", "AACCAA", "ATTA ", "ATTTA", "ATAAAA", "ATGAAA", "AAGCAT", "ATATAA", "AATCAA",
"ATACTA", "ATACAT", "AAAATA", "ATTAAA", "AATTAA", "AATACA", "CATAAA", "AATAAT", "AATCAA", "AATGAA", "ATGGAA", "AATTAA",
"TATAAA", "ATGTAA", "TGTGAA", "AATGCT", "GATATG", "ATGCAA", "AATGTG", "AAAGAT", "ATTAA", "AATAAA", "AATAAT", "CATTG", "ATATA" }; 

char    buffer[256 * 1024];
int     codon_row = 0;
int     codon_length[CODON_ROW] = { 0 };
float   cai[1024] = { 0 };     
float   f_seq[1024] = { 0 };
char*   seq[1024] = { NULL };
int     current_seq_len = 0;  //string length that don't include the '\0' 
float	opt_codon_max[CODON_ROW][CODON_COLOMN] = { 0.0f };
int		codon_max_index[CODON_ROW] = { 0 };
int		codon_special_num = 36; 

float   cal_cai(char*);
int     generic(char*);
void    save_file(char* str, FILE* out);
void    get_codon(FILE* in);
void    get_opt_codon(FILE* in);
void    get_special_codon(FILE* in);
void	get_codon_max();
void	generic_mutation3(char* c, int index, int length);
int		check_dup(char * c, int length, int* x);
int		kmp(const char T[], const char P[], int next[], int len);
char	* origstr = NULL;
int main(int   argc, char   *argv[]) {

	FILE   *in, *out, *codon_in, *codon_opt_in, *codon_special_in;
	char* str = NULL;
/*
	if (argc<9) {
		printf("please input file for analysis as the following rule.\n");
		printf("./codon [seqence file name] [codon file name] [codon_optimazation_file] [codon_special] [seq num] [dup length] [iteration num] [mution ration].\n");
		printf("example: codon.exe seq.txt codon.txt opt_codon.txt special_codon.txt 20 10 100 5.0f\n");
		exit(0);
	}

	SEQ_NUM = atoi(argv[5]);
	DUP_LENGTH = atoi(argv[6]);
	ITERATION_NUM = atoi(argv[7]);
	MUTION_RATIO = atof(argv[8]);
	*/
///*
	argv[1] = "seq.txt";
	argv[2] = "codon.txt";
	argv[3] = "opt_codon.txt";
	argv[4] = "special_codon.txt";	
	SEQ_NUM = 20;
	DUP_LENGTH = 10;
	ITERATION_NUM = 10000;
	MUTION_RATIO = 5.0f;
	//*/

	if ((codon_in = fopen(argv[2], "r")) == NULL) {
		printf("Can't open codon file.");
		exit(0);
	}else {
		get_codon(codon_in);
		fclose(codon_in);
	}

	if ((codon_opt_in = fopen(argv[3], "r")) == NULL) {
		printf("Can't open codon file.");
		exit(0);
	}else {
		get_opt_codon(codon_opt_in);
		fclose(codon_opt_in);
	}

	if ((codon_special_in = fopen(argv[4], "r")) == NULL) {
		printf("Can't open codon file.");
 		exit(0);
	}else {
	//	get_special_codon(codon_special_in);
		fclose(codon_special_in);
	}

	if ((in = fopen(argv[1], "r")) == NULL) {
		printf("Can't open this file.");
		exit(0);
	}
	if ((out = fopen("output.txt", "w")) == NULL) {
		printf("Can't write to file.");
		exit(0);
	}
	getchar();

	get_codon_max();
	while (fgets(buffer, sizeof(buffer), in) != NULL) {
		if (strlen(buffer)<30)
			continue;

		str = (char*)malloc(sizeof(char)*(strlen(buffer) + 1));
		strncpy(str, buffer, strlen(buffer));
		memset(str + strlen(buffer), '\0', 1);
		origstr = str;
		current_seq_len = strlen(str);

		cal_cai(str);
		generic(str);

		free(str);
	}
	fclose(out);
}

void save_file(char* str, FILE* out) {
	fputs(str, out);
	fputs("\n", out);
}


void	get_codon_max() {
	for (int i = 0; i < codon_row; i++) {
		float max = 0;
		for (int j = 0; j < codon_length[i]; j++) {
			if (max < opt_codon[i][j][0]) {
				max = opt_codon[i][j][0];
			}
		}

		for (int j = 0; j < codon_length[i]; j++) {
			opt_codon_max[i][j] = opt_codon[i][j][0] / max;
		}
	}
}

void get_codon(FILE* in) {
	char buf[1024] = { "\0" };
	int i = 0, j = 0;
	char * token = NULL;


	while (fgets(buf, sizeof(buf), in) != NULL) {
		token = strtok(buf, ":");
		if (token != NULL) {
			strncpy(&codon[i][j][0], token, 3);
			memset(&codon[i][j][3], '\0', 1);
			j++;
		}
		else
			break;
		while ((token = strtok(NULL, ":")) != NULL) {
			if (*token == '\n'|| *token == '\r')
				break;
			strncpy(&codon[i][j][0], token, 3);
			memset(&codon[i][j][3], '\0', 1);
			j++;
		}

		codon_length[i] = j;

		i++; j = 0;
		memset(buf, '\0', 1024);
	}
	codon_row = i;
	printf("Loading codon as following:\n");
	for (int m = 0; m<i; m++) {
		for (j = 0; j<CODON_COLOMN; j++) {
			if (*codon[m][j] != NULL)
				printf("[%d %d] %s ", m, j, codon[m][j]);
		}
		printf("[num] %d\n", codon_length[m]);
		printf("\n");
	}
}

void get_special_codon(FILE* in) {
	char buf[1024] = { "\0" };
	int i = 0, j = 0;
	char * token = NULL;

	int a = 0;
	while (fgets(buf, sizeof(buf), in) != NULL) {
		int index = 0;
		a = 0;

		for (int i = 0; i < 1024 ; i++) {
			if (buf[i] == ':' || buf[i] == '\0') {
				strncpy( condon_special[a], &buf[index], i-index   );
				condon_special[a][i - index] = '\0';
				a++; 
				index = i+1;
			}
			if (buf[i] == '\0')
				break;
		}		
	}
	codon_special_num = a;
	printf("Loading codon special as following:\n");
	for (int j = 0; j < a; j++) {
		printf("[%d] %s ", j,condon_special[j]);
	}
	printf("\n");
	
}

void get_opt_codon(FILE* in) {
	char buf[1024] = { "\0" };
	int i = 0, j = 0;
	char * token = NULL;

	char *delim = "\t{}, ";
	while (fgets(buf, sizeof(buf), in) != NULL) {
		token = strtok(buf, delim);
		if (token != NULL) {
			opt_codon[i][j][0] = atof(token);
			if( (token = strtok(NULL, delim)) != NULL)
				opt_codon[i][j][1] = atof(token);
			j++;
		}
		else
			break;

		while ((token = strtok(NULL, delim)) != NULL) {
			opt_codon[i][j][0] = atof(token);
			if ((token = strtok(NULL, delim)) != NULL)
				opt_codon[i][j][1] = atof(token);
			j++;
		}	

		i++; j = 0;		
	}
	
	printf("Loading opt codon as following:\n");
	for (int m = 0; m<i; m++) {
		for (j = 0; j<CODON_COLOMN; j++) {
			if (opt_codon[m][j][0] == 0)
				break;
			printf("[%d,%d], %.2f %.2f ", m, j, opt_codon[m][j][0], opt_codon[m][j][1]);
		}		
		printf("\n");
	}
}



float cal_cai(char* str) {
	/*
	float codon_num[CODON_ROW][CODON_COLOMN] = { 0.0f };
	float codon_max[CODON_ROW] = { 0.0f };

	int i = 0, j = 0, k = 0;
	int len = current_seq_len;

	int sum = 0;
	for (i = 0; i<len; i += 3) {
	char tmp[3] = { "\0" };
	strncpy(tmp, &str[i], 3);
	for (j = 0; j<codon_row; j++)
	for (k = 0; k<codon_length[j]; k++) {
	if (codon[j][k] == NULL)
	break;
	if (strncmp(tmp, codon[j][k], 3) == 0) {
	codon_num[j][k]++;
	sum++;
	}
	}

	}

	for (int j = 0; j<codon_row; j++) {
	for (int k = 0; k<CODON_COLOMN; k++) {
	if (codon_num[j][k]<0.1f)
	break;

	if (codon_max[j] < codon_num[j][k]) {
	codon_max[j] = codon_num[j][k];
	codon_max_index[j] = k;
	}
	printf("[%d:%d] = %0.1f ", j, k, codon_num[j][k]);
	}
	}
	printf("sum = %d\n", sum);

	float cai = 0.0f;
	for (int j = 0; j<codon_row; j++) {
	float a = 0.0f;
	for (int k = 0; k<codon_length[j]; k++) {
	if (codon_num[j][k]<0.1f)
	break;
	float w = codon_num[j][k] / codon_max[j];
	//float w = codon_num[j][k] / codon_num[j][0];

	w = log(w);
	w = codon_num[j][k] * w;
	a += w;
	}
	cai += a;
	}
	cai = cai / sum;
	cai = exp(cai);

	*/

	float cai = 0.0f;

	int i = 0, j = 0, k = 0;
	int len = current_seq_len;

	double sum = 1.0f;

	for (i = 0; i<len; i += 3) {
		char tmp[4] = { "\0" };
		strncpy(tmp, &str[i], 3);
		bool found = false;
		for (j = 0; j < codon_row; j++) {
			for (k = 0; k < codon_length[j]; k++) {
				if (codon[j][k] == NULL)
					continue;
				if (strncmp(tmp, codon[j][k], 3) == 0) {
					sum = sum* opt_codon_max[j][k];
					found = true;
					break;
				}
			}
			if (found)
				break;
		}
	}

	sum = pow(sum, (double)1 / (current_seq_len / 3));
	//	printf("[cai = %0.3f]\n", sum);
	return sum;
}

int get_random(int max_int) {

	int i;
	double j;

	i = rand();
	j = ((double)i / (double)RAND_MAX);
	i = (int)(j * (double)max_int);
	if (i == max_int) {
		i--;
	}

	return i;
}


int get_one_seq(char* str, char* orig) {
	int i = 0;
	int len = current_seq_len;
	strncpy(str, orig, len);

	for (i = 0; i<len; i += 3) {
		char tmp[4] = { "\0" };
		strncpy(tmp, &orig[i], 3);

		for (int j = 0; j<codon_row; j++)
			for (int k = 0; k<codon_length[j]; k++) {
				if (codon[j][k] == NULL)
					break;
				if (strncmp(tmp, codon[j][k], 3) == 0) {
					int a = get_random(codon_length[j]);
					if (a == 6)
						printf("%d \n", a);
					strncpy(&str[i], codon[j][a], 3);
				}
			}
	}
	str[len] = '\0';
	log("%s\n", str);
	return 0;
}

float f(char* str, int length) {
	float b = 1.0f;
	b = cal_cai(str);
	int num_special = 0;
	for (int i = 0; i < codon_special_num; i++) {
		char* tmp = NULL;
		tmp = strstr(str, condon_special[i]);
		while (tmp) {
			num_special++;

			tmp = tmp + strlen(condon_special[i]);
			tmp = strstr(tmp, condon_special[i]);
		}
	}

	if (num_special > 3)
		num_special = 3;

	int a = 0;
	int index = check_dup(str, current_seq_len, &a);
	float x = 0.0f;
	if (index != -1) {

		x = 1.0f;
	}

	return b - 0.1*num_special - 0.05*x; //make sure that f result is positive

}

int get_seq(char * str, int num) {
	log("org seq is: %s\n", str);

	unsigned int seedVal;
	struct timeb timeBuf;
	ftime(&timeBuf);
	seedVal = (((unsigned int)timeBuf.time & 0xFF) +
		(unsigned int)timeBuf.millitm) ^
		(unsigned int)timeBuf.millitm;
	srand(seedVal);

	for (int i = 0; i<num; i++) {
		seq[i] = (char*)malloc(sizeof(char)*(current_seq_len + 1));
		if (seq[i] == NULL) {
			printf("malloc error, exit program.\n");
			return -1;
		}
		f_seq[i] = 0.0f;
		while (f_seq[i] < 0.01f) {
			log("[%d] ",i);
			get_one_seq(seq[i], str);
			f_seq[i] = f(seq[i], current_seq_len);
		}
	}
	return 0;
}


int generic_select(int * x, int* y) {
	int sum = 0;
	float q[1024 +1] = { 0.f };
	int i;
	for (i = 1; i <= SEQ_NUM; i++) {
		q[i] = q[i - 1] + (f_seq[i - 1] * 10);
		sum += (f_seq[i - 1] * 10);
	}

	int a, b;
	double j1, j2;
	a = rand();
	b = rand();
	j1 = ((double)a / (double)RAND_MAX);
	j2 = ((double)b / (double)RAND_MAX);

	a = j1* sum;
	b = j2* sum;

	for (i = 1; i <= SEQ_NUM; i++) {
		if (q[i] > a) {
			break;
		}
	}

	*x = i - 1;

	for (i = 1; i <= SEQ_NUM; i++) {
		if (q[i] > b) {
			break;
		}
	}
	*y = i - 1;

	return 0;
}

int generic_crossover(int x, int y, char* c, int len) {
	char *a = seq[x];
	char *b = seq[y];

	strncpy(c, a, len / 2);
	strncpy(c + len / 2, b + len / 2, len - len / 2);

	c[len] = '\0';
	return 0;
}

int generic_mutation2(char *c) {
	for (int i = 0; i<codon_special_num; i++) {
		char* tmp = NULL;
		tmp = strstr(c, condon_special[i]);
		while (tmp)
		{
			int index = (tmp - c) / sizeof(char);
			int end = index + strlen(condon_special[i]) - 1;
			index = index / 3 * 3;
			end = end / 3 * 3;
			int a = get_random((end - index) / 3 + 1);

			char b[4] = { "\0" };
			strncpy(b, &c[index + a * 3], 3);

			for (int j = 0; j < codon_row; j++) {
				for (int k = 0; k < codon_length[j]; k++) {
					if (codon[j][k] == NULL)
						break;
					if (strncmp(b, codon[j][k], 3) == 0) {
						//get a random codon at j row;
						int d = get_random(codon_length[j]);
						if (d == k) {
							d++;
							if (d >= codon_length[j]) {
								d = d - codon_length[j];
							}
						}
						strncpy(c + index + a * 3, codon[j][d], 3);
					}
				}
			}

			tmp = tmp + strlen(condon_special[i]);
			tmp = strstr(tmp, condon_special[i]);
		}
	}
	return 0;

}

int generic_mutation(char * c) {
	int a = get_random(current_seq_len);
	//int b = get_random(current_seq_len);

	a = a / 3 * 3;
	char tmp[4] = { "\0" };
	strncpy(tmp, &c[a], 3);

	for (int j = 0; j<codon_row; j++)
		for (int k = 0; k<codon_length[j]; k++) {
			if (codon[j][k] == NULL)
				break;
			if (strncmp(tmp, codon[j][k], 3) == 0) {
				int index = get_random(codon_length[j]);
				strncpy(c + a, codon[j][index], 3);
			}
		}
	return 0;
}

int replace_seq(char* a, int b, float new_value) {

	strncpy(seq[b], a, current_seq_len);
	f_seq[b] = new_value;
	return 0;
}

int output() {
	int maxindex = 0;
	float max = 0.0f;

	for (int i = 0; i < SEQ_NUM; i++) {
		if (f_seq[i] >max) {
			max = f_seq[i];
			maxindex = i;
		}
	}

	int num_special = 0;
	for (int i = 0; i < codon_special_num; i++) {
		char* tmp = NULL;
		tmp = strstr(seq[maxindex], condon_special[i]);
		while (tmp) {
			num_special++;

			tmp = tmp + strlen(condon_special[i]);
			tmp = strstr(tmp, condon_special[i]);
		}
	}
	printf("[max value] %0.4f, [special num] %d\n", max, num_special);
	printf("[seq string]: %s\n", seq[maxindex]);
	return maxindex;
}

void check_codon(int index) {

	for (int i = 0; i < current_seq_len; i += 3) {
		char tmp[4] = { "\0" };
		strncpy(tmp, &seq[index][i], 3);

		bool bfound = false;
		for (int j = 0; j < codon_row; j++) {
			for (int k = 0; k < codon_length[j]; k++) {
				if (codon[j][k] == NULL)
					break;
				if (strncmp(tmp, codon[j][k], 3) == 0) {
					char* a = origstr + i;
					int x = 0;
					for (x = 0; x < codon_length[j]; x++) {
						if (strncmp(a, codon[j][x], 3) == 0) {
							bfound = true;
							break;
						}
					}
					if (x == codon_length[j]) {
						printf("Can't find string at %d.\n ", i);
						return;
					}
					if (bfound)
						break;
				}
			}
			if (bfound)
				break;
		}
		//	if (!bfound) {
		//		printf("can't find %s at %d.\n", tmp,i);
		//	}
	}
	int max = 0;
	int x = 0;
	int res = check_dup(seq[index], current_seq_len, &x);
	if (res != -1) {
		printf("have duplicate codon\n");
	}else	{
		printf("[dup] 0");
	}
}

int  generic(char * str) {
	int mutation_num = 20;
	int ret = get_seq(str, SEQ_NUM);
	if (ret == -1)
		return -1;
	output();
	unsigned int seedVal;
	struct timeb timeBuf;

	ftime(&timeBuf);
	seedVal = (((unsigned int)timeBuf.time & 0xFF) +
		(unsigned int)timeBuf.millitm) ^
		(unsigned int)timeBuf.millitm;
	srand(seedVal);
	log("seed = %d\n", seedVal);

	char* c = (char*)malloc(sizeof(char)*(current_seq_len + 1));
	if (c == NULL) {
		printf("malloc error, exit program.\n");
		return NULL;
	}

	for (int i = 0; i<ITERATION_NUM; i++) {
		printf("[%d] ", i);		
		int x = 0, y = 0;
		generic_select(&x, &y);    
		log("the x = %d ,y = %d\n", x, y);
		generic_crossover(x, y, c, current_seq_len);
		log("crossover, str = %s\n", c);
		mutation_num = (int)(MUTION_RATIO / 100 * current_seq_len / 3);
		for (int m = 0; m < mutation_num; m++) {
			generic_mutation(c);
		}
		log("1st mutation, str = %s\n", c);
		generic_mutation2(c);
		log("2rd mutation, str = %s\n", c);
		int a = 0;
		int index = check_dup(c, current_seq_len, &a);
		if (index != -1) {

			/*
			for (int m = 0; m < DUP_LENGTH; m++) {
			printf("%c ", *(c + index + m));
			}
			printf("| ");
			for (int m = 0; m < DUP_LENGTH; m++) {
			printf("%c ", *(c + a + m));
			}
			printf("\n");
			*/

			generic_mutation3(c, index, DUP_LENGTH);
			generic_mutation3(c, a, DUP_LENGTH);

		}
		log("3rd mutation, str = %s\n", c);

		float f_new = f(c, current_seq_len);
		//printf("This new seq f = %0.3f\n", f_new);
		float n = f_seq[x] < f_seq[y] ? f_seq[x] : f_seq[y];
		int min = f_seq[x] < f_seq[y] ? x : y;
		log("after mutation, new sequece is %.4f, the ori %d = %.4f, the orig %d is %.4f\n", f_new, x, f_seq[x], y, f_seq[y]);
		if (f_new > n) {
			replace_seq(c, min, f_new);
		}
//		int aa = output();
	}
	printf("\n");

	int index = output();
	check_codon(index);

	if (c != NULL) {
		free(c);
		c = NULL;
	}

	for (int i = 0; i < SEQ_NUM; i++) {
		if (seq[i] != NULL) {
			free(seq[i]);
			seq[i] = NULL;
		}
	}
	getchar();
	return 0;
}
void generic_mutation3(char* c, int index, int length) {
	index = index / 3 * 3;
	int a = get_random(length / 3 + 1);

	char b[4] = { "\0" };
	strncpy(b, &c[index + a * 3], 3);

	for (int j = 0; j < codon_row; j++) {
		for (int k = 0; k < codon_length[j]; k++) {
			if (codon[j][k] == NULL)
				break;
			if (strncmp(b, codon[j][k], 3) == 0) {
				//get a random codon at j row;
				int d = get_random(codon_length[j]);
				if (d == k) {
					d++;
					if (d >= codon_length[j]) {
						d = d - codon_length[j];
					}
				}
				strncpy(c + index + a * 3, codon[j][d], 3);
			}
		}
	}
}


int check_dup(char * c, int length, int* x) {

	int* next = (int*)malloc(sizeof(int) * DUP_LENGTH);
	for (int i = 0; i < length + 1 - 2 * DUP_LENGTH; i++) {
		char * index = c;
		int ret = kmp(index + i + DUP_LENGTH, index + i, next, current_seq_len - (i + DUP_LENGTH));
		if (ret != -1) {
			free(next);
			*x = i;
			return ret + i + DUP_LENGTH;
		}
	}

	char *str = (char*)malloc(sizeof(char)*(current_seq_len + 1));

	for (int i = 0; i < length; i++) {
		str[i] = c[length - 1 - i];
	}
	str[length] = '\0';
	for (int i = 0; i < length - DUP_LENGTH; i++) {
		char * index = str;

		int ret = kmp(index, c + i, next, length - i - 1);
		if (ret != -1) {
			free(next);
			free(str);
			*x = i;
			return length - ret - DUP_LENGTH;
		}
	}
	/*
	for (int i = 0; i < length; i++) {
	if (c[length - 1 - i] == 'A')
	str[i] = 'T';
	else if (c[length - 1 - i] == 'T')
	str[i] = 'A';
	else if (c[length - 1 - i] == 'C')
	str[i] = 'G';
	else if (c[length - 1 - i] == 'G')
	str[i] = 'C';
	}
	str[length] = '\0';
	for (int i = 0; i < length - DUP_LENGTH; i++) {
	char * index = c;
	int ret = kmp(index + i + DUP_LENGTH, index + i, next, current_seq_len - (i + DUP_LENGTH));
	if (ret != -1) {
	free(next);
	free(str);
	*x = i;
	return ret;
	}
	}

	for (int i = 0; i < length; i++) {
	if (c[i] == 'A')
	str[i] = 'T';
	else if (c[i] == 'T')
	str[i] = 'A';
	else if (c[i] == 'C')
	str[i] = 'G';
	else if (c[i] == 'G')
	str[i] = 'C';
	}
	str[length] = '\0';
	for (int i = 0; i < length - DUP_LENGTH; i++) {
	char * index = c;
	int ret = kmp(index + i + DUP_LENGTH, index + i, next, current_seq_len - (i + DUP_LENGTH));
	if (ret != -1) {
	free(next);
	free(str);
	*x = i;
	return ret;
	}
	}


	free(str);
	*/
	free(next);
	return -1;
}
void makeNext(const char P[], int next[])
{
	int q, k;
	int m = DUP_LENGTH;
	next[0] = 0;
	for (q = 1, k = 0; q < m; ++q)
	{
		while (k > 0 && P[q] != P[k])
			k = next[k - 1];
		if (P[q] == P[k])
		{
			k++;
		}
		next[q] = k;
	}
}

int kmp(const char T[], const char P[], int next[], int len)
{

	int i, q;
	int m = DUP_LENGTH;
	if (len <= DUP_LENGTH)
		return -1;

	makeNext(P, next);
	for (i = 0, q = 0; i < len; ++i) {
		while (q > 0 && P[q] != T[i])
			q = next[q - 1];
		if (P[q] == T[i]) {
			q++;
		}
		if (q == m) {
			return i - m + 1;
		}
	}
	return -1;
}
