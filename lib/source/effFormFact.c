/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on arXiv:1802.03174
effFormFact.c

Early versions of this code were developed by M. Peiró
and E. Gerstmayr. 

Using arXiv:1203.3542 and then the unit conventions of
arXiv:1308.6288.

nuclear models (lib/nuclear_responses/<model>.txt, sources in the README.md there):
  - Fitz2012: arXiv:1203.3542
  - Anand2013: DMFormFactor v6 of arXiv:1308.6288 with its default one-body density matrices
  - McElvain2022: DMFormFactor v6 with the one-body density matrices of arXiv:2208.07945 (Elastic)
  - Khaleq2023: arXiv:2311.15764
The coefficients of the nuclear responses are read from RAPIDD_for_DM/lib/nuclear_responses
(read_nuclear_responses).
  
##########################################*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "effFormFact.h"
#include "phys_consts.h"

#define SQR(X) ((X)*(X))
#define ABS(X) ((X) > 0 ? (X) : (-(X)))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define VERBOSE 0
//static const double B[199] = {7.309204019430945, 7.315148063632845, 7.463325019607083, 7.6147938871555, 7.755089818125042, 7.883294539176039, 8.000702993350325, 8.108833807657158, 8.209031480136986, 8.3024131489692, 8.389893277928199, 8.472220967542217, 8.550013802433732, 8.623785234903465, 8.693965947425994, 8.760920351728842, 8.824959341635076, 8.886350205384998, 8.945324389668524, 9.002083633107096, 9.056804854096278, 9.109644079567662, 9.16073962909474, 9.210214715909135, 9.258179587520988, 9.304733299880532, 9.349965197591908, 9.393956156603048, 9.436779633620716, 9.478502557213162, 9.51918608842395, 9.558886273191588, 9.597654604555899, 9.635538509243588, 9.67258177054525, 9.708824897263094, 9.744305446800508, 9.779058309088583, 9.813115956930426, 9.846508667436574, 9.879264718482132, 9.911410563505587, 9.942970987464498, 9.973969246344618, 10.004427192269866, 10.034365385968668, 10.063803198106669, 10.092758900789109, 10.121249750361008, 10.149292062484639, 10.176901280347233, 10.20409203674345, 10.230878210684573, 10.257272979106373, 10.283288864179013, 10.308937776662704, 10.334231055701403, 10.359179505401958, 10.383793428507028, 10.408082657436173, 10.432056582939442, 10.455724180581774, 10.479094035253476, 10.502174363881696, 10.524973036499986, 10.547497595817205, 10.569755275412948, 10.591753016674263, 10.613497484577346, 10.634995082408002, 10.656251965505897, 10.677274054109738, 10.698067045373433, 10.718636424617069, 10.738987475870694, 10.759125291764038, 10.779054782810428, 10.79878068612926, 10.818307573647477, 10.837639859817196, 10.856781808883582, 10.87573754173417, 10.894511042358516, 10.913106163944565, 10.931526634636146, 10.94977606297412, 10.967857943041912, 10.98577565933463, 11.003532491369489, 11.021131618054, 11.038576121827079, 11.055868992587243, 11.07301313142095, 11.090011354143211, 11.106866394661843, 11.123580908175787, 11.140157474217359, 11.156598599547458, 11.172906720912314, 11.189084207669627, 11.205133364291555, 11.221056432751421, 11.236855594800623, 11.252532974141772, 11.268090638503775, 11.283530601624058, 11.298854825143007, 11.31406522041522, 11.329163650241938, 11.344151930528827, 11.359031831872871, 11.373805081082079, 11.388473362631377, 11.403038320057911, 11.417501557298769, 11.431864639973963, 11.446129096617378, 11.460296419858205, 11.474368067555226, 11.488345463886233, 11.502230000394697, 11.516023036995687, 11.529725902942959, 11.543339897758996, 11.556866292129701, 11.570306328765385, 11.583661223229507, 11.59693216473668, 11.610120316921268, 11.623226818577903, 11.636252784375113, 11.649199305543284, 11.662067450538002, 11.67485826567989, 11.687572775771905, 11.700211984695045, 11.71277687598339, 11.725268413379336, 11.737687541369812, 11.750035185704311, 11.762312253895407, 11.774519635702557, 11.78665820359976, 11.798728813227795, 11.810732303831621, 11.822669498683512, 11.8345412054925, 11.846348216800648, 11.858091310366689, 11.86977124953746, 11.881388783607653, 11.892944648168307, 11.904439565444413, 11.915874244622149, 11.927249382166009, 11.938565662126289, 11.949823756437222, 11.961024325206159, 11.972168016994056, 11.983255469087634, 11.994287307763496, 12.005264148544454, 12.016186596448382, 12.027055246229848, 12.037870682614756, 12.048633480528263, 12.059344205316199, 12.070003412960189, 12.080611650286736, 12.091169455170428, 12.101677356731486, 12.112135875527835, 12.122545523741907, 12.13290680536229, 12.143220216360467, 12.153486244862735, 12.163705371317539, 12.173878068658274, 12.184004802461796, 12.1940860311027, 12.204122205903564, 12.214113771281239, 12.224061164889333, 12.233964817757013, 12.243825154424204, 12.253642593073353, 12.263417545657816, 12.27315041802699, 12.282841610048289, 12.292491515726056, 12.30210052331753, 12.311669015445874, 12.321197369210488, 12.33068595629454, 12.340135143069901, 12.349545290699542, 12.358916755237406, 12.368249887725913, 12.377545034291105};
static const double B2[199] = {53.424463397665484, 53.511391192871365, 55.70122034829307, 57.98508594386077, 60.1414180871867, 62.14633279140276, 64.01124838780486, 65.75318572020367, 67.38819784188003, 68.93006409617665, 70.39030921502477, 71.77852812286199, 73.10273602180732, 74.36967177773901, 75.58504389500276, 76.75372540933662, 77.87990738151218, 78.967219972746, 80.01882843639856, 81.03750973745466, 82.0257141651819, 82.98561525640216, 83.91915055206684, 84.82805511314919, 85.7138892747903, 86.57806178190565, 87.42184919617989, 88.24641227218031, 89.05280985351874, 89.84201072709645, 90.61490378604407, 91.37230678381056, 92.11497390835305, 92.84360236311615, 93.55883810788428, 94.26128088571572, 94.95148864054606, 95.62998141255446, 96.29724478416256, 96.95373293790357, 97.59987137784583, 98.23605935837014, 98.86267205756074, 99.48006252702822, 100.0885634454287, 100.68848869912614, 101.28013481022202, 101.86378222945778, 102.43969650918277, 103.0081293696137, 103.56931966993314, 104.1234942943311, 104.67086896186036, 105.21164896790573, 105.7460298641481, 106.27419808310337, 106.79633151262334, 107.31260002513994, 107.82316596590574, 108.32818460402363, 108.82780454965014, 109.3221681404024, 109.81141179968499, 110.29566636937392, 110.77505741905175, 111.24970553376973, 111.71972658211985, 112.18523196622834, 112.64632885512965, 113.10312040284239, 113.55570595234829, 114.004181226565, 114.44863850730505, 114.88916680312778, 115.32585200690761, 115.75877704387659, 116.18802201082838, 116.61366430711834, 117.03577875803836, 117.4544377310985, 117.86971124570546, 118.28166707668599, 118.69037085207165, 119.09588614552487, 119.49827456375947, 119.89759582928102, 120.29390785874756, 120.68726683722922, 121.07772728862403, 121.46534214246957, 121.85016279737096, 122.23223918125206, 122.6116198086208, 122.98835183502534, 123.36248110886856, 123.73405222073288, 124.10310855036089, 124.4696923114243, 124.83384459420756, 125.19560540632183, 125.55501371155978, 125.91210746699205, 126.26692365840206, 126.61949833414788, 126.96986663753641, 127.31806283778657, 127.6641203596574, 128.0080718118093, 128.34994901396323, 128.68978302292092, 129.02760415750114, 129.3634420224485, 129.6973255313644, 130.02928292870916, 130.3593418109198, 130.687529146687, 131.01387129643095, 131.3383940310148, 131.66112254973103, 131.9820814975954, 132.3012949819798, 132.61878658861536, 132.93457939699383, 133.24869599519465, 133.56115849416372, 133.8719885414683, 134.18120733455092, 134.48883563350418, 134.794893773388, 135.0994016761086, 135.40237886187757, 135.70384446027015, 136.00381722089793, 136.30231552371404, 136.5993573889646, 136.89496048680155, 137.18914214657124, 137.48191936579116, 137.7733088188281, 138.06332686528935, 138.35198955813806, 138.63931265154508, 138.92531160848552, 139.21000160809177, 139.4933975527718, 139.77551407510143, 140.05636554449987, 140.3359660736959, 140.61432952499396, 140.89146951634606, 141.16739942723774, 141.44213240439515, 141.71568136731835, 141.98805901364946, 142.25927782437944, 142.5293500689009, 142.79828780991141, 143.06610290817346, 143.33280702713537, 143.5984116374187, 143.86292802117651, 144.1263672763268, 144.38874032066573, 144.6500578958649, 144.91033057135584, 145.1695687481066, 145.4277826622934, 145.6849823888706, 145.9411778450436, 146.19637879364635, 146.45059484642755, 146.70383546724844, 146.95610997519495, 147.20742754760658, 147.45779722302555, 147.7072279040677, 147.95572836021915, 148.2033072305589, 148.44997302641212, 148.695734133934, 148.9405988166285, 149.184575217802, 149.42767136295555, 149.66989516211638, 149.91125441211088, 150.15175679878146, 150.391409899148, 150.63022118351608, 150.86819801753364, 151.10534766419707, 151.34167728580942, 151.57719394588997, 151.81190461103944, 152.0458161527594, 152.27893534922882, 152.51126888703922, 152.7428233628879, 152.97360528523205, 153.20362107590438};

/*########################################################################################
Simpy the HELM FORM Factor for tests and simplicity
########################################################################################*/
double fsi(int A, double Er){
    double s = 0.9;                                         // nuclear skin thickness in fm
    double a = 0.52;                                        // surface thickness of the nucleus in fm
    double c=1.23*pow((double)A*1.0,(1.0/3.0))-0.6;             // parameter of the Fermi distribution
    double rn=sqrt(c*c-5.0*s*s+7.0/3.0*M_PI*M_PI*a*a);      // effective radius
    double q=(6.92e-3)*sqrt((double)A*1.0)*sqrt(Er);            // transfered momentum
    
    return 3.0*((sin(q*rn)-q*rn*cos(q*rn))/pow(q*rn,3.0))*exp(-pow(q*s,2.0)/2.0);
}

/*########################################################################################
Generic polynomial function used to parametrize each nuclear response in arXiv:1308.6288v1
########################################################################################*/
double f(double Er, int A, int Z, double a0, double a1, double a2, double a3, double a4,
         double a5, double a6, double a7, double a8, double a9, double a10)
{
    double mn = nucleus_mass(A, Z); // GeV
    double q2 = 2.*mn*Er*1.e-6;//GeV
    double b2 = B2[A-1];

    const double y   = 0.25*q2*b2;
    const double y2  = y*y;
    const double y3  = y2*y;
    const double y4  = y2*y2;
    const double y5  = y2*y3;
    const double y6  = y3*y3;
    const double y7  = y3*y4;
    const double y8  = y4*y4;
    const double y9  = y4*y5;
    const double y10 = y5*y5;

    double FormFact = exp(-2.*y) * ( a0 +  a1*y + a2*y2 + a3*y3 + a4*y4 + a5*y5 + a6*y6 + a7*y7 + a8*y8 + a9*y9 + a10*y10);

    return FormFact;
}

/*########################################################################################
Nuclear responses
The fit coefficients of the nuclear responses of each nuclear model are tabulated in
RAPIDD_for_DM/lib/nuclear_responses/<model>.txt (format and sources in the README.md there)
and are read by read_nuclear_responses.
Labels of the responses ("Symbol1" and "Symbol2"):
- M = M
- Sp = Sigma'
- Spp = Sigma''
- D = Delta
- Ppp = Phi''
- Pp = Phi', Pt (they enter O_12 and O_13 as Pp + Pt)
Labels of the nucleons ("N1" and "N2"): p, n for the models in the proton-neutron basis, +, - for
the models in the isospin basis, with +- = (p +- n)/2 (the basis of each model is declared in
nuclear_tables and returned by nuclear_model_basis).
########################################################################################*/

#define N_COEFF 11       /* a_0 ... a_10 of f */
#define LABEL_LEN 8      /* labels of at most 7 characters */
#define LINE_LEN 4096

enum nuclear_basis {BASIS_PN, BASIS_ISO};

/* one row of <model>.txt: squared response F_{O1,O2}^{N1,N2} = norm * f(E_R; a_0 ... a_10) */
struct nuclear_response {
	int A;
	int Z;
	char o1[LABEL_LEN], o2[LABEL_LEN], n1[LABEL_LEN], n2[LABEL_LEN];
	double norm;    /* pn: the prefactor column; iso: 16 pi / (2 j_N + 1) from the j_N column */
	double a[N_COEFF];
};

struct nuclear_table {
	const char *name;           /* name of the model (NR_framework of FF); the file is <name>.txt */
	enum nuclear_basis basis;
	struct nuclear_response *rows;
	int n_rows;
};

/* The nuclear models: name and basis (proton-neutron or isospin) of the nucleon labels and of the column c */
static struct nuclear_table nuclear_tables[] = {
	{"Fitz2012", BASIS_PN, NULL, 0},       /* arXiv:1203.3542 */
	{"Anand2013", BASIS_ISO, NULL, 0},     /* DMFormFactor v6, default OBDMs */
	{"McElvain2022", BASIS_ISO, NULL, 0},  /* DMFormFactor v6, OBDMs of the Elastic repository */
	{"Khaleq2023", BASIS_PN, NULL, 0},     /* arXiv:2311.15764 */
};
#define N_NUCLEAR_TABLES ((int)(sizeof(nuclear_tables) / sizeof(nuclear_tables[0])))

static int nuclear_loaded = 0;

static int valid_nucleon(enum nuclear_basis basis, const char *label){
	if (basis == BASIS_PN) return strcmp(label, "p") == 0 || strcmp(label, "n") == 0;
	return strcmp(label, "+") == 0 || strcmp(label, "-") == 0;
}

static int read_label(const char **s, char *label){
	/*Reads the next whitespace-separated label of *s (at most LABEL_LEN - 1 characters).*/
	int n = 0;
	char fmt[16];
	snprintf(fmt, sizeof(fmt), "%%%ds%%n", LABEL_LEN - 1);
	if (sscanf(*s, fmt, label, &n) != 1 || ((*s)[n] != '\0' && !isspace((unsigned char)(*s)[n]))) return 0;
	*s += n;
	return 1;
}

static int read_numbers(const char **s, double *x, int count){
	/*Reads count numbers from *s.*/
	for (int k = 0; k < count; k++) {
		int n = 0;
		if (sscanf(*s, "%lf%n", &x[k], &n) != 1) return 0;
		*s += n;
	}
	return 1;
}

static int only_spaces(const char *s){
	while (*s != '\0') {
		if (!isspace((unsigned char)*s)) return 0;
		s++;
	}
	return 1;
}

static FILE * open_table(const char *dir, const char *name, char *path, size_t path_len){
	snprintf(path, path_len, "%s/%s.txt", dir, name);
	FILE *fp = fopen(path, "r");
	if (fp == NULL) fprintf(stderr, "read_nuclear_responses: cannot open %s\n", path);
	return fp;
}

static const char * next_row(FILE *fp, char *line, int *line_number){
	/*Returns the next data line of fp (comments # and empty lines are skipped), NULL at the end
	of the file or for a line longer than LINE_LEN (*line_number < 0).*/
	while (fgets(line, LINE_LEN, fp) != NULL) {
		(*line_number)++;
		if (strchr(line, '\n') == NULL && !feof(fp)) {
			*line_number = -*line_number;
			return NULL;
		}
		const char *s = line;
		while (isspace((unsigned char)*s)) s++;
		if (*s != '\0' && *s != '#') return s;
	}
	return NULL;
}

static void * grow(void *rows, int n_rows, int *capacity, size_t size){
	/*Returns rows with room for one more element (NULL if out of memory).*/
	if (n_rows < *capacity) return rows;
	*capacity = (*capacity == 0) ? 64 : 2 * *capacity;
	return realloc(rows, (size_t)*capacity * size);
}

static int same_response(const struct nuclear_response *r, int A, int Z, const char *Symbol1,
						 const char *Symbol2, const char *N1, const char *N2){
	return r->A == A && r->Z == Z && strcmp(Symbol1, r->o1) == 0 && strcmp(Symbol2, r->o2) == 0
		&& strcmp(N1, r->n1) == 0 && strcmp(N2, r->n2) == 0;
}

static int read_nuclear_table(const char *dir, const struct nuclear_table *t,
							  struct nuclear_response **rows_out, int *n_out){
	/*Reads <dir>/<t->name>.txt. Each data line is
	  A Z O1 O2 N1 N2 c a_0 ... a_10
	with N1, N2 = p, n and c the prefactor of the response (pn models) or N1, N2 = +, - and
	c = j_N the spin of the nucleus, the prefactor being 16 pi / (2 j_N + 1) (iso models).
	Each row (A, Z, O1, O2, N1, N2) must appear once.
	Output: 0 on success, 2 if the file is missing, 3 if a line is malformed or repeats a row,
	4 if out of memory.*/
	char path[1024], line[LINE_LEN];
	FILE *fp = open_table(dir, t->name, path, sizeof(path));
	if (fp == NULL) return 2;
	struct nuclear_response *rows = NULL;
	int n_rows = 0, capacity = 0, line_number = 0, err = 0;
	const char *s;
	while ((s = next_row(fp, line, &line_number)) != NULL) {
		struct nuclear_response r;
		double c;
		int n = 0;
		if (sscanf(s, "%d %d%n", &r.A, &r.Z, &n) != 2 || r.A < 1 || r.A > 199) { err = 3; break; }
		s += n;
		if (!read_label(&s, r.o1) || !read_label(&s, r.o2) || !read_label(&s, r.n1) || !read_label(&s, r.n2)
			|| !valid_nucleon(t->basis, r.n1) || !valid_nucleon(t->basis, r.n2)
			|| !read_numbers(&s, &c, 1) || !read_numbers(&s, r.a, N_COEFF) || !only_spaces(s)) { err = 3; break; }
		r.norm = (t->basis == BASIS_PN) ? c : 16*M_PI/(2*c + 1);
		for (int i = 0; i < n_rows && err == 0; i++) {
			if (same_response(&rows[i], r.A, r.Z, r.o1, r.o2, r.n1, r.n2)) err = 5;
		}
		if (err != 0) break;
		struct nuclear_response *tmp = grow(rows, n_rows, &capacity, sizeof(*rows));
		if (tmp == NULL) { err = 4; break; }
		rows = tmp;
		rows[n_rows++] = r;
	}
	if (err == 0 && line_number < 0) err = 3;
	if (err == 3) fprintf(stderr, "read_nuclear_responses: %s, line %d is malformed\n", path, abs(line_number));
	if (err == 5) {
		fprintf(stderr, "read_nuclear_responses: %s, line %d repeats an earlier row\n", path, line_number);
		err = 3;
	}
	fclose(fp);
	if (err != 0) {
		free(rows);
		return err;
	}
	*rows_out = rows;
	*n_out = n_rows;
	return 0;
}

int read_nuclear_responses(const char *dir){
	/*Reads the nuclear responses of all the nuclear models from dir (e.g. .../lib/nuclear_responses):
	<model>.txt for each model of nuclear_tables (Fitz2012, Anand2013, McElvain2022, Khaleq2023). The tables are replaced only if all the files are read (otherwise the
	previous ones, if any, are kept).
	Output: 0 on success, 2 if a file is missing, 3 if a line is malformed, 4 if out of memory.*/
	struct nuclear_response *rows[N_NUCLEAR_TABLES];
	int n_rows[N_NUCLEAR_TABLES];
	int err = 0, n_read = 0;
	for (; n_read < N_NUCLEAR_TABLES; n_read++) {
		err = read_nuclear_table(dir, &nuclear_tables[n_read], &rows[n_read], &n_rows[n_read]);
		if (err != 0) break;
	}
	if (err != 0) {
		for (int t = 0; t < n_read; t++) free(rows[t]);
		return err;
	}
	for (int t = 0; t < N_NUCLEAR_TABLES; t++) {
		free(nuclear_tables[t].rows);
		nuclear_tables[t].rows = rows[t];
		nuclear_tables[t].n_rows = n_rows[t];
	}
	nuclear_loaded = 1;
	return 0;
}

int nuclear_responses_loaded(void){
	/*Returns 1 if the nuclear responses have been read (read_nuclear_responses), 0 otherwise.*/
	return nuclear_loaded;
}

static void require_nuclear_responses(void){
	if (!nuclear_loaded) {
		fprintf(stderr, "Error: the nuclear responses have not been read (read_nuclear_responses).\n");
		exit(1);
	}
}

static const struct nuclear_table * nuclear_table(const char *name){
	for (int t = 0; t < N_NUCLEAR_TABLES; t++) {
		if (strcmp(name, nuclear_tables[t].name) == 0) return &nuclear_tables[t];
	}
	return NULL;
}

static double FF_table(const struct nuclear_table *t, int A, int Z, const char *Symbol1, const char *Symbol2,
					   const char *N1, const char *N2, double Er){
	/*Response of the row of t with the same A, Z, Symbol1, Symbol2, N1 and N2 (0 if none).*/
	for (int i = 0; i < t->n_rows; i++) {
		const struct nuclear_response *r = &t->rows[i];
		if (same_response(r, A, Z, Symbol1, Symbol2, N1, N2)) {
			return r->norm * f(Er, A, Z, r->a[0], r->a[1], r->a[2], r->a[3], r->a[4], r->a[5],
							   r->a[6], r->a[7], r->a[8], r->a[9], r->a[10]);
		}
	}
	return 0.0;
}

/*########################################################################################
Form factor function for all isotopes: response of the nucleus (A, Z).
NR_framework (nuclear model):
- Fitz2012: fits of arXiv:1203.3542 (p, n basis; outdated, they diverge at E_R > 100 keV)
- Anand2013: DMFormFactor v6 (arXiv:1308.6288) with its default one-body density matrices (+, - basis)
- McElvain2022: DMFormFactor v6 with the one-body density matrices of the Elastic repository
  (arXiv:2208.07945) (+, - basis)
- Khaleq2023: fits of arXiv:2311.15764 (SN100PN interaction, Xe only; p, n basis)
########################################################################################*/
int nuclear_model_basis(const char* NR_framework){
	/*Basis of the nuclear model: 0 proton-neutron (labels p, n), 1 isospin (labels +, -),
	-1 if NR_framework is not a nuclear model.*/
	const struct nuclear_table *t = nuclear_table(NR_framework);
	if (t == NULL) return -1;
	return (t->basis == BASIS_ISO) ? 1 : 0;
}

double FF(char* NR_framework, int A, int Z, char* Symbol1, char * Symbol2, char*N1, char * N2, double Er){
	const struct nuclear_table *t = nuclear_table(NR_framework);
	if (t == NULL) {
		printf("Error: NR_framework %s not recognized. Please use Fitz2012, Anand2013, McElvain2022 or Khaleq2023.\n", NR_framework);
		exit(1);
	}
	require_nuclear_responses();
	return FF_table(t, A, Z, Symbol1, Symbol2, N1, N2, Er);
}

/*########################################################################################
General form factors
########################################################################################*/

/*The following form factors are those with zeroth order dependence on the velocity.
Other dependences should be included separately like FormFact_v2 function, the coefficient of
(v_T^perp)^2 = v^2 - vmin^2 (also for inelastic scattering, arXiv:1409.0536): difrate_isotope_dEr
integrates it with halo(vmin, 2) - vmin^2 halo(vmin, 0).
The expressions for the prefactors have been extrapolated from eq. 37 and 38 of arXiv:1308.6288v1.*/

double FormFact_v0(char* NR_framework, int A, int Z, int i, int j, char* N1, char* N2, double Er, double mchi, double jchi) {
    
    double mnucleus = nucleus_mass(A, Z);
	double q2 = 2. * mnucleus * Er * 1.e-6;
	double q4 = q2 * q2;
    double mN = mnucleon;
	double mN2 = mN * mN;
	double mN4 = mN2 * mN2;
    double Cj = jchi * (jchi + 1.) / 12.;

    if (i == 1 && j == i) {
        return FF(NR_framework, A, Z, "M", "M", N1, N2, Er);
    }
    if (i == 3 && j == i) {
        return q4 / (4. * mN4) * FF(NR_framework, A, Z, "Ppp", "Ppp", N1, N2, Er);
	}
    if (i == 4 && j == i) {
        return Cj * (FF(NR_framework, A, Z, "Spp", "Spp", N1, N2, Er) + FF(NR_framework, A, Z, "Sp", "Sp", N1, N2, Er));
    }
    if (i == 5 && j == i) {
        return 4. * Cj * q4 / mN4 * FF(NR_framework, A, Z, "D", "D", N1, N2, Er);
    }
    if (i == 6 && j == i) {
        return Cj * q4 / mN4 * FF(NR_framework, A, Z, "Spp", "Spp", N1, N2, Er);
    }
    if (i == 8 && j == i) {
        return 4 * Cj * q2 / mN2 * FF(NR_framework, A, Z, "D", "D", N1, N2, Er);
    }
    if (i == 9 && j == i) {
        return Cj * q2 / mN2 * FF(NR_framework, A, Z, "Sp", "Sp", N1, N2, Er);
    }
    if (i == 10 && j == i) {
        return  q2 / (4. * mN2) * FF(NR_framework, A, Z, "Spp", "Spp", N1, N2, Er);
    }
    if (i == 11 && j == i) {
        return 4 * Cj * q2 / mN2 * FF(NR_framework, A, Z, "M", "M", N1, N2, Er);
    }
	if (i == 12 && j == i) {
		return Cj * (q2 / mN2) * (FF(NR_framework, A, Z, "Ppp", "Ppp", N1, N2, Er) + FF(NR_framework, A, Z, "Pp", "Pp", N1, N2, Er) + FF(NR_framework, A, Z, "Pt", "Pt", N1, N2, Er));
	}
	if (i == 13 && j == i) {
		return Cj * q4 / mN4 * (FF(NR_framework, A, Z, "Pp", "Pp", N1, N2, Er) + FF(NR_framework, A, Z, "Pt", "Pt", N1, N2, Er));
	}
    if (i == 15 && j == i) {
        return Cj * (q2/mN2)*(q4/mN4) * FF(NR_framework, A, Z, "Ppp", "Ppp", N1, N2, Er);
    }
    if (i == 1 && j == 3) {
        return q2 / mN2 * FF(NR_framework, A, Z, "M", "Ppp", N1, N2, Er);
    }
    if (i == 4 && j == 5) {
        return 4. * Cj * q2 / mN2 * FF(NR_framework, A, Z, "Sp", "D", N1, N2, Er);
    }
    if ((i == 4 && j == 6) || (i == 6 && j == 4)) {
        return Cj * q2 / mN2 * FF(NR_framework, A, Z, "Spp", "Spp", N1, N2, Er);
    }
    if (i == 9 && j == 8) {
        return -4. * Cj * q2 / mN2 * FF(NR_framework, A, Z, "Sp", "D", N1, N2, Er);
    }
	if (i == 11 && j == 12) {
		return 4. * Cj * q2 / mN2 * FF(NR_framework, A, Z, "M", "Ppp", N1, N2, Er);
	}
	if (i == 11 && j == 15) {
		return -4. * Cj * q4 / mN4 * FF(NR_framework, A, Z, "M", "Ppp", N1, N2, Er);
	}
	if (i == 12 && j == 15) {
		return -2 * Cj * q4 / mN4 * FF(NR_framework, A, Z, "Ppp", "Ppp", N1, N2, Er);
	}
    return 0.0;
}

double FormFact_v2(char* NR_framework, int A, int Z, int i, int j, char* N1, char* N2, double Er, double mchi, double jchi) {

    double mnucleus = nucleus_mass(A, Z);
    double mN = mnucleon;
	double mN2 = mN * mN;
	double mN4 = mN2 * mN2;
    double q2 = 2. * mnucleus * Er * 1.e-6;
    double q4 = q2 * q2;
    double Cj = jchi * (jchi + 1.) / 12.;

    if (i == 3 && j == i) {
        return 1./ 8. * q2 / mN2 * FF(NR_framework, A, Z, "Sp", "Sp", N1, N2, Er);
    }
    if (i == 5 && j == i) {
        return 4 * Cj * q2 / mN2 * FF(NR_framework, A, Z, "M", "M", N1, N2, Er);
    }
    if (i == 7 && j == i) {
        return 1. / 8. * FF(NR_framework, A, Z, "Sp", "Sp", N1, N2, Er);
    }
    if (i == 8 && j == i) {
        return 4 * Cj * FF(NR_framework, A, Z, "M", "M", N1, N2, Er);
    }
	if (i == 12 && j == i) {
		return Cj * (FF(NR_framework, A, Z, "Spp", "Spp", N1, N2, Er) + 0.5 * FF(NR_framework, A, Z, "Sp", "Sp", N1, N2, Er));
	}
    if (i == 13 && j == i) {
        return Cj * q2 / mN2 * FF(NR_framework, A, Z, "Spp", "Spp", N1, N2, Er);
    }
    if (i == 14 && j == i) {
        return Cj * q2 / (2. * mN2) * FF(NR_framework, A, Z, "Sp", "Sp", N1, N2, Er);
    }
    if (i == 15 && j == i) {
        return Cj * q4 / mN4 * 0.5 * FF(NR_framework, A, Z, "Sp", "Sp", N1, N2, Er);
    }
	if (i == 12 && j == 15) {
		return -Cj * q2 / mN2 * FF(NR_framework, A, Z, "Sp", "Sp", N1, N2, Er);
	}
    return 0.0;
}
