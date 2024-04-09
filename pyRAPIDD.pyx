cdef extern from "coeffs_eft.h":
    void read_coeffs(const char* fname)
    double Cp(int i)
    double Cn(int i)
    double give_mass()
    void set_any_coeffs(double C, int i)
    void set_any_Ncoeff(double input, int i, char* nucleon)
    void set_med_mass(double med_mass)
    double give_med_mass()

cdef extern from "halo.h":
    double halo (double vmin, int i)
    void read_halo(char*path)

cdef extern from "effFormFact.h":
    double FormFact_v0(int A,int Z,int i, int j, char * N1, char * N2, double Er, double mchi, double jchi)

def py_FormFact_v0(int A, int Z, int i, int j, bytes N1, bytes N2, float Er, float mchi, float jchi) -> float:
    return FormFact_v0 ( A, Z, i, j, N1, N2, Er, mchi, jchi)

cdef extern from "MenenFF.h":
    double FormFactMen_v0(int A,int Z, int i, int j, char * N1, char * N2, double Er, double mchi, double jchi)

def py_FormFactMenen(int A, int Z, int i, int j, bytes N1, bytes N2, float Er, float mchi, float jchi) -> float:
    return FormFactMen_v0( A, Z, i, j, N1, N2, Er, mchi, jchi)


def py_read_coeffs(bytes fname) -> None:
    read_coeffs(fname)

def py_give_mass() -> double:
    return give_mass()

def py_Cp(int i) -> double:
    return Cp(i)

def py_Cn(int i) -> double:
    return Cn(i)

def py_set_any_coeffs(float C, int i):
    set_any_coeffs(C, i)

def py_set_any_Ncoeff(float input, int i, bytes nucleon):
    set_any_Ncoeff(input, i, nucleon)

def py_set_med_mass(float med_mass):
    set_med_mass(med_mass)
    
def py_med_mass():
    return give_med_mass()
    
    
def py_read_halo(bytes path) -> void:
    read_halo(path) 
    return


cdef extern from "efficiency_curve.h":
    void read_efficiency(char* path)
    double efficiency ( double E_r)

def py_read_efficiency(bytes path) -> void:
    read_efficiency(path)
    return

def py_efficiency(float Er) -> double:
    return efficiency(Er)


cdef extern from "binning_general.h":
    double counts_bin_python( char* halo_path, double rhochi, double mass, double exposure, double E1, double E2, char*model, char*Target, char*basis)
    double counts_bin_python_noread(double rhochi, double mass, double exposure, double E1, double E2, char*model, char*Target)

def py_counts_bin(bytes halo_path, float rhochi, float mass, float exposure, float E1, float E2, bytes model, bytes Target, bytes basis) -> float:
    return counts_bin_python( halo_path, rhochi,  mass,  exposure, E1,  E2, model, Target, basis)

def py_counts_bin_noread(float rhochi, float mass, float exposure, float E1, float E2, bytes model, bytes Target) -> float:
    return counts_bin_python_noread(rhochi,  mass,  exposure, E1,  E2, model, Target)


cdef extern from "specific_exps.h":
    double bin_Xenon1T( double rhochi, double mass, double E1, double E2, char* model, char* basis)
    double bin_DS50( double rhochi, double mass, double E1, double E2, char* model, char* basis)
    double bin_LZ( double rhochi, double mass, double E1, double E2, char* model, char* basis)
    double bin_DS20k( double rhochi, double mass, double E1, double E2, char* model, char* basis)
    double res_fn_lux_nr(int A, int Z, double ERnr)
    double counts_effres_bin_Xenon1T( double rhochi, double mass, double E1, double E2, char*model, char*basis)
    double counts_effres_bin_LZ( double rhochi, double mass, double E1, double E2, char*model, char*basis)
    void read_DS50_LEFF(char* path)
    double DS50_LEFF ( double E_r)
    double res_fn_DS20_nr(double ERnr)
    double counts_effres_bin_DS50( double rhochi, double mass, double E1, double E2, char*model, char*basis)
    double counts_effres_bin_DS20k( double rhochi, double mass, double E1, double E2, char*model, char*basis)

def py_bin_Xe1t( float rhochi, float mass, float E1, float E2, bytes model, bytes basis) -> float:
    return bin_Xenon1T(rhochi, mass, E1, E2, model, basis)

def py_bin_DS50( float rhochi, float mass, float E1, float E2, bytes model, bytes basis ) -> float:
    return bin_DS50(rhochi, mass, E1, E2, model, basis)

def py_bin_LZ( float rhochi, float mass, float E1, float E2, bytes model, bytes basis) -> float:
    return bin_LZ(rhochi, mass, E1, E2, model, basis)

def py_bin_DS20k( float rhochi, float mass, float E1, float E2, bytes model, bytes basis) -> float:
    return bin_DS20k(rhochi, mass, E1, E2, model, basis)

def py_read_DS50_LEFF( bytes path) -> void :
    return read_DS50_LEFF(path)

def py_DS50_LEFF( float E_r) -> float:
    return DS50_LEFF(E_r)

def py_res_fn_DS20_nr(float ERnr) -> float:
    return res_fn_DS20_nr(ERnr)

def py_bin_DS20k_res(float rhochi, double mass, double E1, double E2, bytes model, bytes basis) -> float:
    return counts_effres_bin_DS20k(rhochi, mass, E1, E2, model, basis)

def py_bin_DS50_res(float rhochi, double mass, double E1, double E2, bytes model, bytes basis) -> float:
    return counts_effres_bin_DS50(rhochi, mass, E1, E2, model, basis)
    
cdef extern from "difRateGen.h":
    double difrate_dER_python(double rhochi, double mass, double logenergy, char*model, char*Target, char*ISO_switch)

cdef extern from "difRateGen.h":
    double difrate_dER_python_2(double rhochi, double mass, double energy, char*model, char*Target, char*ISO_switch)


def py_difrate_dER(float rhochi, float mass, float logenergy, bytes model, bytes Target, bytes ISO_switch) -> float:
    return difrate_dER_python(rhochi, mass, logenergy, model, Target, ISO_switch)

def py_difrate_dER_notlog(float rhochi, float mass, float energy, bytes model, bytes Target, bytes ISO_switch) -> float:
    return difrate_dER_python_2(rhochi, mass, energy, model, Target, ISO_switch)


def check_reslux(int A, int Z, float ERnr)->float:
    return res_fn_lux_nr(A, Z, ERnr)

def py_bin_LZ_res(float rhochi, double mass, double E1, double E2, bytes model, bytes basis)-> float:
    return counts_effres_bin_LZ(rhochi, mass, E1, E2, model, basis)

def py_bin_Xenon1T_res(float rhochi, double mass, double E1, double E2, bytes model, bytes basis)-> float:
    return counts_effres_bin_Xenon1T(rhochi, mass, E1, E2, model, basis)
