/**
 * Implementation of a class which handles the calculation of the slowing-down frequency nu_s.
 * The A^p component of the collision operator is given by -p*nu_s.
*/

/**
 * The calculations of the electron-ion contribution are based on Eq (2.31) from
 * L Hesslow et al., Generalized collision operator for fast electrons
 * interacting with partially ionized impurities, J Plasma Phys 84 (2018).
 * The relativistic thermal ee contribution is based on the expressions given in
 * Pike & Rose, Dynamical friction in a relativistic plasma, Phys Rev E 89 (2014).
 * The non-linear contribution corresponds to the isotropic component of the
 * non-relativistic operator following Rosenbluth, Macdonald & Judd, Phys Rev (1957),
 * and is described in doc/notes/theory.pdf Appendix B.
 */
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/FVMException.hpp"
#include <cmath>

using namespace DREAM;

/**
 * Mean excitation energy atomic data for ions from 
 * Sauer, Sabin, Oddershede J Chem Phys 148, 174307 (2018)
 */

// List of mean excitation energies in units of eV
const real_t SlowingDownFrequency::MEAN_EXCITATION_ENERGY_DATA[MAX_Z][MAX_Z] = {
/* H  */ { 14.99, NAN,   NAN,   NAN,   NAN,   NAN,   NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* He */ { 42.68, 59.88, NAN,   NAN,   NAN,   NAN,   NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Li */ { 33.1,  108.3, 134.5, NAN,   NAN,   NAN,   NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Be */ { 42.2,  76.9,  205.0, 240.2, NAN,   NAN,   NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* B  */ { 52.6,  82.3,  136.9, 330.4, 374.6, NAN,   NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* C  */ { 65.9,  92.6,  134.8, 214.2, 486.2, 539.5, NAN,   NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* N  */ { 81.6,  107.4, 142.4, 200.2, 308.7, 672.0, 734.3, NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* O  */ { 97.9,  125.2, 157.2, 202.2, 278.6, 420.7, 887.8, 959.0,  NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* F  */ { 116.5, 144.0, 176.4, 215.6, 272.3, 370.2, 550.0, 1133.5, 1213.7, NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Ne */ { 137.2, 165.2, 196.9, 235.2, 282.8, 352.6, 475.0, 696.8,  1409.2, 1498.4, NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Na */ { 125.7, 189.2, 220.4, 256.8, 301.9, 358.7, 443.5, 593.3,  861.2,  1715.6, 1813.9, NAN,    NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Mg */ { 128.0, 173.7, 246.8, 282.5, 324.3, 376.7, 443.8, 544.8,  724.8,  1043.2, 2051.5, 2158.8, NAN,    NAN,    NAN,    NAN,    NAN,    NAN},
/* Al */ { 132.2, 172.7, 225.8, 310.8, 351.0, 398.8, 459.2, 537.4,  656.4,  869.6,  1242.7, 2417.2, 2533.5, NAN,    NAN,    NAN,    NAN,    NAN},
/* Si */ { 140.8, 177.2, 221.2, 283.1, 381.4, 426.5, 480.6, 549.7,  640.1,  778.6,  1027.9, 1459.8, 2813.0, 2938.3, NAN,    NAN,    NAN,    NAN},
/* P  */ { 151.6, 185.3, 225.2, 274.3, 345.9, 458.5, 508.8, 569.7,  648.2,  751.7,  911.2,  1199.7, 1694.6, 3238.8, 3373.1, NAN,    NAN,    NAN},
/* S  */ { 162.4, 195.7, 232.8, 277.3, 332.4, 414.4, 542.1, 598.0,  666.2,  754.6,  872.2,  1054.5, 1384.9, 1947.0, 3694.5, 3837.8, NAN,    NAN},
/* Cl */ { 174.9, 206.8, 242.9, 284.1, 333.8, 395.5, 488.6, 632.1,  694.0,  769.9,  869.1,  1001.8, 1208.2, 1583.7, 2217.2, 4180.2, 4332.5, NAN},
/* Ar */ { 188.7, 219.5, 254.0, 293.7, 339.4, 394.9, 463.9, 568.6,  728.8,  797.0,  881.1,  991.6,  1140.3, 1372.6, 1796.0, 2505.0, 4695.9, 4857.2 }};

const real_t SlowingDownFrequency::MEAN_EXCITATION_ENERGY_FUNCTION_D[MAX_NE] = {0, 0.00, 0.24, 0.34, 0.41, 0.45, 0.48, 0.50, 0.51, 0.52, 0.55, 0.57, 0.58, 0.59};
const real_t SlowingDownFrequency::MEAN_EXCITATION_ENERGY_FUNCTION_S_0[MAX_NE] = {0, 0.30, 1.51, 2.32, 3.13, 3.90, 4.67, 5.44, 6.21, 6.97, 8.10, 9.08, 10.03, 10.94};
// according to Berger et al., J of the ICRU os19 22 (1984), all neutral ions with Z >= 19 have I~10*Z eV (8.8 to 11.1 eV) 
const real_t SlowingDownFrequency::HIGH_Z_EXCITATION_ENERGY_PER_Z = 10.0; 
const real_t SlowingDownFrequency::HYDROGEN_MEAN_EXCITATION_ENERGY = 14.99; // Mean excitation energy for neutral H

// Tabulated values of the integral 'I' of 'bremsIntegrand' from 0 to 'X' used in the bremsstrahlung formula. 
const real_t SlowingDownFrequency::BREMS_INTEGRAL_X[BREMS_INTEGRAL_N] = {0.00000000000000,  0.00149639561318158,    0.0119711649054527,	0.0404026815559027,	0.0957693192436213,	0.187049451647698,	0.323221452447222,	0.513263695321283,	0.766154553948971,	1.09087240200937,	1.49639561318158,	1.99170256114469,	2.58577161957778,	3.28758116215994,	4.10610956257026,	5.05033519448785,	6.12923643159177,	7.35179164756112,	8.72697921607499,	10.2637775108125,	11.9711649054527,	13.8581197736746,	15.9336204891575,	18.2066454255803,	20.6861729566222,	23.3811814559622,	26.3006492972795,	29.4535548542531,	32.8488765005621,	36.4955926098856,	40.4026815559028,	44.5791217122926,	49.0338914527341,	53.7759691509066,	58.8143331804889,	64.1579619151604,	69.8158337285999,	75.7969269944867,	82.1102200864998,	88.7646913783184,	95.7693192436213,	103.133082056088,	110.864958189397,	118.973926017228,	127.468963913260,	136.359050251172,	145.653163404643,	155.360281747351,	165.489383652978,	176.049447495200,	187.049451647698,	198.498374484150,	210.405194378236,	222.778889703635,	235.628438834025,	248.962820143086,	262.791012004497,	277.121992791937,	291.964740879085,	307.328234639620,	323.221452447222,	339.653372675569,	356.632973698340,	374.169233889215,	392.271131621873,	410.947645269992,	430.207753207252,	450.060433807332,	470.514665443912,	491.579426490669,	513.263695321283,	535.576450309434,	558.526669828799,	582.123332253060,	606.375415955894,	631.291899310981,	656.881760691999,	683.153978472628,	710.117531026547,	737.781396727435,	766.154553948971,	795.245981064834,	825.064656448703,	855.619558474258,	886.919665515177,	918.973955945140,	951.791408137825,	985.381000466912,	1019.75171130608,	1054.91251902901,	1090.87240200937,	1127.64033862086,	1165.22530723714,	1203.63628623190,	1242.88225397881,	1282.97218885156,	1323.91506922382,	1365.71987346928,	1408.39557996160,	1451.95116707448,	1496.39561318158,	1541.73789665660,	1587.98699587320,	1635.15188920507,	1683.24155502589,	1732.26497170933,	1782.23111762908,	1833.14897115881,	1885.02751067220,	1937.87571454293,	1991.70256114469,	2046.51702885114,	2102.32809603598,	2159.14474107287,	2216.97594233550,	2275.83067819754,	2335.71792703268,	2396.64666721459,	2458.62587711696,	2521.66453511346,	2585.77161957778,	2650.95610888358,	2717.22698140455,	2784.59321551437,	2853.06378958672,	2922.64768199528,	2993.35387111372,	3065.19133531573,	3138.16905297498,	3212.29600246516,	3287.58116215994,	3364.03351043300,	3441.66202565802,	3520.47568620868,	3600.48347045866,	3681.69435678164,	3764.11732355129,	3847.76134914130,	3932.63541192535,	4018.74849027711,	4106.10956257026,	4194.72760717849,	4284.61160247547,	4375.77052683488,	4468.21335863039,	4561.94907623570,	4656.98665802448,	4753.33508237040,	4851.00332764715,	4950.00037222840,	5050.33519448785,	5152.01677279915,	5255.05408553599,	5359.45611107206,	5465.23182778102,	5572.39021403657,	5680.94024821237,	5790.89090868212,	5902.25117381948,	6015.03002199813,	6129.23643159176,	6244.87938097405,	6361.96784851867,	6480.51081259930,	6600.51725158963,	6721.99614386332,	6844.95646779406,	6969.40720175554,	7095.35732412141,	7222.81581326538,	7351.79164756112,	7482.29380538230,	7614.33126510260,	7747.91300509571,	7883.04800373530,	8019.74523939505,	8158.01369044864,	8297.86233526975,	8439.30015223206,	8582.33611970925,	8726.97921607499,	8873.23841970297,	9021.12270896687,	9170.64106224036,	9321.80245789713,	9474.61587431084,	9629.09028985519,	9785.23468290385,	9943.05803183050,	10102.5693150088,	10263.7775108125,	10426.6915976152,	10591.3205537906,	10757.6733577124,	10925.7589877542,	11095.5864222898,	11267.1646396928,	11440.5026183369,	11615.6093365958,	11792.4937728432,	11971.1649054527,	12151.6317127980,	12333.9031732528,	12517.9882651908,	12703.8959669856,	12891.6352570110,	13081.2151136406,	13272.6445152481,	13465.9324402071,	13661.0878668914,	13858.1197736746,	14057.0371389305,	14257.8489410326,	14460.5641583547,	14665.1917692705,	14871.7407521535,	15080.2200853776,	15290.6387473163,	15503.0057163435,	15717.3299708326,	15933.6204891575,	16151.8862496918,	16372.1362308091,	16594.3794108832,	16818.6247682878,	17044.8812813965,	17273.1579285829,	17503.4636882209,	17735.8075386840,	17970.1984583459,	18206.6454255803,	18445.1574187609,	18685.7434162614,	18928.4123964555,	19173.1733377168,	19420.0352184189,	19669.0070169357,	19920.0977116407,	20173.3162809077,	20428.6717031103,	20686.1729566222,	20945.8290198171,	21207.6488710686,	21471.6414887505,	21737.8158512364,	22006.1809369000,	22276.7457241150,	22549.5191912550,	22824.5103166938,	23101.7280788050,	23381.1814559622,	23662.8794265393,	23946.8309689098,	24233.0450614474,	24521.5306825258,	24812.2968105188,	25105.3524237999,	25400.7065007428,	25698.3680197213,	25998.3459591090,	26300.6492972795,	26605.2870126066,	26912.2680834640,	27221.6014882253,	27533.2962052642,	27847.3612129543,	28163.8054896694,	28482.6380137832,	28803.8677636693,	29127.5037177014,	29453.5548542531,	29782.0301516982,	30112.9385884103,	30446.2891427632,	30782.0907931304,	31120.3525178857,	31461.0832954028,	31804.2921040553,	32149.9879222169,	32498.1797282613,	32848.8765005621,	33202.0872174931,	33557.8208574279,	33916.0863987403,	34276.8928198038,	34640.2490989921,	35006.1642146790,	35374.6471452381,	35745.7068690432,	36119.3523644678,	36495.5926098856,	36874.4365836704,	37255.8932641958,	37639.9716298355,	38026.6806589632,	38416.0293299525,	38808.0266211772,	39202.6815110109,	39600.0029778272,	40000.0000000000};
const real_t SlowingDownFrequency::BREMS_INTEGRAL_I[BREMS_INTEGRAL_N] = {0.00000000000000,	0.00149583618521427,	0.0119355270535575,	0.0400017530787400,	0.0935690203127896,	0.178961336450119,	0.300288707338392,	0.459144861209056,	0.654784791245859,	0.884674105514979,	1.14517796911550,	1.43218481044150,	1.74156328864825,	2.06944254526409,	2.41235219009931,	2.76726847914547,	3.13160560969440,	3.50317918552194,	3.88015856850110,	4.26101750896821,	4.64448780986870,	5.02951806276985,	5.41523800649615,	5.80092830684630,	6.18599522056955,	6.56994949650123,	6.95238887051496,	7.33298356692374,	7.71146429402562,	8.08761229878136,	8.46125111724258,	8.83223972023391,	9.20046680729849,	9.56584604658333,	9.92831209518956,	10.2878172646886,	10.6443287211107,	10.9978261287288,	11.3482996632334,	11.6957483331185,	12.0401785588655,	12.3816029682817,	12.7200393735085,	13.0555099010744,	13.3880402511668,	13.7176590662469,	14.0443973923792,	14.3682882193367,	14.6893660877630,	15.0076667535216,	15.3232269008955,	15.6360838975852,	15.9462755855226,	16.2538401024166,	16.5588157297040,	16.8612407632109,	17.1611534033701,	17.4585916622920,	17.7535932853712,	18.0461956854391,	18.3364358877491,	18.6243504843210,	18.9099755963692,	19.1933468437179,	19.4744993202513,	19.7534675745762,	20.0302855951830,	20.3049867994873,	20.5776040262143,	20.8481695306579,	21.1167149824115,	21.3832714652141,	21.6478694786062,	21.9105389411274,	22.1713091948218,	22.4302090108496,	22.6872665960256,	22.9425096001327,	23.1959651238746,	23.4476597273512,	23.6976194389554,	23.9458697646016,	24.1924356972104,	24.4373417263820,	24.6806118482007,	24.9222695751203,	25.1623379458874,	25.4008395354656,	25.6377964649281,	25.8732304112923,	26.1071626172719,	26.3396139009281,	26.5706046652019,	26.8001549073153,	27.0282842280277,	27.2550118407401,	27.4803565804384,	27.7043369124699,	27.9269709411487,	28.1482764181859,	28.3682707509433,	28.5869710105082,	28.8043939395889,	29.0205559602313,	29.2354731813566,	29.4491614061222,	29.6616361391056,	29.8729125933158,	30.0830056970316,	30.2919301004717,	30.4997001822976,	30.7063300559535,	30.9118335758450,	31.1162243433613,	31.3195157127418,	31.5217207967931,	31.7228524724581,	31.9229233862403,	32.1219459594881,	32.3199323935413,	32.5168946747438,	32.7128445793250,	32.9077936781544,	33.1017533413720,	33.2947347428976,	33.4867488648225,	33.6778065016864,	33.8679182646436,	34.0570945855194,	34.2453457207624,	34.4326817552930,	34.6191126062525,	34.8046480266558,	34.9892976089485,	35.1730707884731,	35.3559768468460,	35.5380249152476,	35.7192239776280,	35.8995828738314,	36.0791103026407,	36.2578148247442,	36.4357048656281,	36.6127887183953,	36.7890745465137,	36.9645703864950,	37.1392841505071,	37.3132236289210,	37.4863964927943,	37.6588102962936,	37.8304724790568,	38.0013903684971,	38.1715711820515,	38.3410220293728,	38.5097499144699,	38.6777617377947,	38.8450642982801,	39.0116642953274,	39.1775683307471,	39.3427829106535,	39.5073144473124,	39.6711692609480,	39.8343535815039,	39.9968735503649,	40.1587352220377,	40.3199445657921,	40.4805074672646,	40.6404297300239,	40.7997170771017,	40.9583751524865,	41.1164095225851,	41.2738256776493,	41.4306290331707,	41.5868249312440,	41.7424186418996,	41.8974153644058,	42.0518202285430,	42.2056382958476,	42.3588745608305,	42.5115339521672,	42.6636213338617,	42.8151415063861,	42.9660992077934,	43.1164991148087,	43.2663458438936,	43.4156439522906,	43.5643979390425,	43.7126122459922,	43.8602912587591,	44.0074393076973,	44.1540606688305,	44.3001595647692,	44.4457401656088,	44.5908065898072,	44.7353629050461,	44.8794131290727,	45.0229612305256,	45.1660111297424,	45.3085666995512,	45.4506317660465,	45.5922101093476,	45.7333054643448,	45.8739215214264,	46.0140619271946,	46.1537302851648,	46.2929301564528,	46.4316650604461,	46.5699384754644,	46.7077538394046,	46.8451145503756,	46.9820239673181,	47.1184854106144,	47.2545021626856,	47.3900774685769,	47.5252145365314,	47.6599165385547,	47.7941866109661,	47.9280278549413,	48.0614433370433,	48.1944360897452,	48.3270091119405,	48.4591653694462,	48.5909077954952,	48.7222392912197,	48.8531627261262,	48.9836809385605,	49.1137967361658,	49.2435128963310,	49.3728321666303,	49.5017572652580,	49.6302908814508,	49.7584356759062,	49.8861942811910,	50.0135693021443,	50.1405633162712,	50.2671788741320,	50.3934184997219,	50.5192846908456,	50.6447799194853,	50.7699066321602,	50.8946672502835,	51.0190641705092,	51.1430997650753,	51.2667763821398,	51.3900963461120,	51.5130619579772,	51.6356754956158,	51.7579392141179,	51.8798553460913,	52.0014261019656,	52.1226536702891,	52.2435402180228,	52.3640878908291,	52.4842988133543,	52.6041750895076,	52.7237188027350,	52.8429320162888,	52.9618167734931,	53.0803750980036,	53.1986089940647,	53.3165204467605,	53.4341114222652,	53.5513838680833,	53.6683397132943,	53.7849808687845,	53.9013092274819,	54.0173266645851,	54.1330350377863,	54.2484361874939,	54.3635319370502,	54.4783240929464,	54.5928144450323,	54.7070047667249,	54.8208968152128,	54.9344923316572,	55.0477930413898,	55.1608006541084,	55.2735168640680,	55.3859433502714,	55.4980817766525,	55.6099337922633,	55.7215010314511,	55.8327851140382,	55.9437876454946,	56.0545102171142,	56.1649544061804,	56.2751217761365,	56.3850138767479,	56.4946322442645,	56.6039784015837,	56.7130538584037,	56.8218601113806,	56.9303986442812,	57.0386709281328,	57.1466784213728,	57.2544225699913,	57.3619048076811,	57.4691265559732,	57.5760892243828,	57.6827942105403,	57.7892429003336};


/**
 * Constructor
 */
SlowingDownFrequency::SlowingDownFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                CoulombLogarithm *lnLee,CoulombLogarithm *lnLei,
                enum OptionConstants::momentumgrid_type mgtype,  struct collqty_settings *cqset)
                : CollisionFrequency(g,u,ih,lnLee,lnLei,mgtype,cqset){
    hasIonTerm = false;
    gsl_ad_w = gsl_integration_workspace_alloc(1000);
    bremsSpline = gsl_spline_alloc(gsl_interp_steffen, BREMS_INTEGRAL_N);
    gsl_spline_init(bremsSpline, BREMS_INTEGRAL_X, BREMS_INTEGRAL_I, BREMS_INTEGRAL_N);
    gsl_acc = gsl_interp_accel_alloc();
}


/**
 * Destructor.
 */
SlowingDownFrequency::~SlowingDownFrequency(){
    gsl_integration_workspace_free(gsl_ad_w);
    gsl_spline_free(bremsSpline);
}


/**
 * Evaluates the matched Bethe formula according to Eq (2.31) in the Hesslow paper.
 * Modification: Moved the -beta^2 contribution inside the interpolation term in order
 * to preserve positivity of the contribution.
 */
real_t SlowingDownFrequency::evaluateScreenedTermAtP(len_t iz, len_t Z0, real_t p, OptionConstants::collqty_collfreq_mode collfreq_mode){
    len_t Z = ionHandler->GetZ(iz); 
    len_t ind = ionHandler->GetIndex(iz,Z0);
    if (Z==Z0)
        return 0;
    real_t p2 = p*p;
    real_t gamma = sqrt(1+p2);
    real_t beta2 = p2/(1+p2);
    real_t h = (p2/sqrt(1+gamma))/atomicParameter[ind];
    real_t NBound = Z - Z0;

    if (collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)
        return NBound*log(1+pow(h*exp(-beta2),kInterpolate))/kInterpolate ;
    else 
        return NBound*log(exp(1)+h*exp(-beta2));

// previous expression as given in the paper:
//    return NBound*(log(1+pow(h,kInterpolate))/kInterpolate-beta2) ;
}


/**
 * Returns the mean excitation energy of ion species with index iz and charge state Z0.
 * When data is missing, it uses the approximate analytical 2-parameter formula (8) from 
 * Sauer, Sabin, Oddershede J Chem Phys 148, 174307 (2018)
 */
real_t SlowingDownFrequency::GetAtomicParameter(len_t iz, len_t Z0){
    len_t Z = ionHandler->GetZ(iz);

    real_t I;
    real_t D_N;
    real_t S_N0;

    if (Z == Z0){
        return NAN;
    }
    if (Z <= MAX_Z){ /* use tabulated data */
        I = MEAN_EXCITATION_ENERGY_DATA[Z-1][Z0];
    }else{ /* use the formula instead */
        len_t Ne = Z-Z0;
        if (Ne <= MAX_NE){
            D_N = MEAN_EXCITATION_ENERGY_FUNCTION_D[Ne-1]; 
            S_N0 = MEAN_EXCITATION_ENERGY_FUNCTION_S_0[Ne-1];
        }else{
            D_N = MEAN_EXCITATION_ENERGY_FUNCTION_D[MAX_NE-1]; 
            S_N0 = Ne - sqrt(Ne*HIGH_Z_EXCITATION_ENERGY_PER_Z / HYDROGEN_MEAN_EXCITATION_ENERGY); // S_N0: for a neutral atom with Z=N
        }
        real_t A_N = (1-D_N) * (1-D_N);
        real_t B_N = 2*(1-D_N) * (Ne*D_N - S_N0);
        real_t C_N = (Ne*D_N - S_N0) * (Ne*D_N - S_N0);

        I = HYDROGEN_MEAN_EXCITATION_ENERGY * (A_N*Z*Z + B_N*Z + C_N);
    }
    return I / Constants::mc2inEV;
}


/**
 * Helper function to calculate the partial contribution to evaluateAtP from free electrons
 */
real_t SlowingDownFrequency::evaluateElectronTermAtP(len_t ir, real_t p,OptionConstants::collqty_collfreq_mode collfreq_mode){
    if (collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
        if(p==0)
            return 0;
        real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
        real_t gamma = sqrt(1+p*p);
        real_t gammaMinusOne = p*p/(gamma+1); // = gamma-1
        real_t Theta = T_cold[ir] / Constants::mc2inEV;
        
        real_t M = gamma*gamma* evaluatePsi1(ir,p) - Theta * evaluatePsi0(ir,p);
        M +=  (Theta*gamma - 1) * p * exp( -gammaMinusOne/Theta );
        M /= gamma*gamma*evaluateExp1OverThetaK(Theta,2.0);
        return  M;
    } else 
        return 1;
    
}


/**
 * Helper function for integral term in bremsstrahlung formula 
 */
real_t bremsIntegrand(real_t x, void*){
    return log(1+x)/x;
}


/**
 * Evaluates the bremsstrahlung stopping power formula. Using the non-screened 
 * formula given as (4BN) in H W Koch and J W Motz, Rev Mod Phys 31, 920 (1959).
 */
real_t SlowingDownFrequency::evaluateBremsstrahlungTermAtP(len_t iz, len_t /*Z0*/, real_t p, OptionConstants::eqterm_bremsstrahlung_mode brems_mode, OptionConstants::collqty_collfreq_type /*collfreq_type*/){
    if(brems_mode != OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER)
        return 0;
    else if(p==0)
        return 0;

    real_t preFactor = constPreFactor * Constants::alpha / (4*M_PI);
    len_t Z = ionHandler->GetZ(iz); 
    real_t gamma = sqrt(1+p*p);
    real_t gammaMinus1OverP = p/(gamma+1);
    preFactor *= Z*Z * gammaMinus1OverP;

    // The formula from ecritpaper Eq (18)
    // return preFactor * 4*M_PI*( 0.35+0.2*log(gamma) );

    real_t integralTerm;
    real_t X = 2*p*(gamma+p);
    real_t X_max = BREMS_INTEGRAL_X[BREMS_INTEGRAL_N-1];
    if(X<=X_max)
        integralTerm = gsl_spline_eval(bremsSpline,X,gsl_acc);
    else {
        real_t error;
        gsl_function GSL_Func;
        
        GSL_Func.function = &(bremsIntegrand);
        GSL_Func.params = nullptr; 
        real_t epsabs=0, epsrel=3e-3;
        gsl_integration_qag(&GSL_Func,X_max,X,epsabs,epsrel,gsl_ad_w->limit,QAG_KEY,gsl_ad_w,&integralTerm,&error);
        integralTerm += BREMS_INTEGRAL_I[BREMS_INTEGRAL_N-1]; // add value of integral from 0 to X_max
    }
    real_t logTerm = log(gamma+p);
    real_t Term1 = (4.0/3.0) * (3*gamma*gamma+1)/(gamma*p) * logTerm;
    real_t Term2 = -(8*gamma+6*p)/(3*gamma*p*p)*logTerm*logTerm - 4/3;
    real_t Term3 = 2.0/(gamma*p) * integralTerm;

    return preFactor*(Term1+Term2+Term3);
}


/**
 * Helper function to calculate a partial contribution to evaluateAtP
 */
real_t SlowingDownFrequency::evaluateDDTElectronTermAtP(len_t ir, real_t p,OptionConstants::collqty_collfreq_mode collfreq_mode){
    if ( (collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) && p){
        real_t T_cold = unknowns->GetUnknownData(id_Tcold)[ir];
        real_t gamma = sqrt(1+p*p);
        real_t gammaMinusOne = p*p/(gamma+1); // = gamma-1
        real_t Theta = T_cold / Constants::mc2inEV;
        real_t DDTheta = 1/Constants::mc2inEV;

        real_t Psi0 = evaluatePsi0(ir,p);
        real_t Psi1 = evaluatePsi1(ir,p);
        real_t Psi2 = evaluatePsi2(ir,p);
        real_t DDTPsi0 = DDTheta / (Theta*Theta) * (Psi1-Psi0);
        real_t DDTPsi1 = DDTheta / (Theta*Theta) * (Psi2-Psi1);

        real_t Denominator = gamma*gamma*evaluateExp1OverThetaK(Theta,2.0);
        real_t DDTDenominator = DDTheta/(Theta*Theta) * (gamma*gamma*evaluateExp1OverThetaK(Theta,1.0) - (1-2*Theta) * Denominator);

        real_t Numerator = gamma*gamma* Psi1 - Theta * Psi0;
        Numerator +=  (Theta*gamma - 1) * p * exp( -gammaMinusOne/Theta );
        
        real_t DDTNumerator = gamma*gamma* DDTPsi1 - (DDTheta * Psi0 + Theta * DDTPsi0 );
        DDTNumerator +=  (gamma + gammaMinusOne/(Theta*Theta) *(Theta*gamma - 1) ) * DDTheta * p * exp( -gammaMinusOne/Theta ) ;

        return  DDTNumerator  / Denominator - Numerator*DDTDenominator /(Denominator*Denominator);
    } else 
        return 0;
    
}


/**
 * Evaluates the purely momentum dependent prefactor for nu_s
 */
real_t SlowingDownFrequency::evaluatePreFactorAtP(real_t p,OptionConstants::collqty_collfreq_mode collfreq_mode){
    if(p==0) 
        return 0; 
    else if (collfreq_mode != OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_ULTRA_RELATIVISTIC)
        return constPreFactor * (1+p*p)/(p*p*p);
    else
        return constPreFactor / p;
}


/**
 * evaluates lim_{p\to 0} p^3nu_s, for use in the evaluation of the 
 * particle flux through the p=0 boundary. It corresponds to the 
 * evaluateAtP calculation evaluated at 0, if we set preFactor to 
 * constPreFactor instead of evaluatePreFactorAtP.
 */
real_t SlowingDownFrequency::GetP3NuSAtZero(len_t ir){
    real_t *ncold = unknowns->GetUnknownData(id_ncold);    
    real_t ntarget = ncold[ir];
    if (isNonScreened)
        ntarget += ionHandler->GetBoundElectronDensity(ir);

    real_t preFactor = constPreFactor;
    real_t lnLee0 = lnLambdaEE->evaluateAtP(ir,0);
    real_t p3nuS0 = lnLee0 * evaluateElectronTermAtP(ir,0,collQtySettings->collfreq_mode) * ntarget;

    // The partially screened term below will vanish identically for the 
    // formulas given in the Hesslow paper; we keep it for completeness in
    // case the model is changed in the future.
    if(isPartiallyScreened)
        for(len_t iz = 0; iz<nZ; iz++)
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                len_t ind = ionIndex[iz][Z0];
                p3nuS0 += evaluateScreenedTermAtP(iz,Z0,0,collQtySettings->collfreq_mode) * ionDensities[ir][ind];
            }
    p3nuS0 *= preFactor;
    return p3nuS0;
}


/**
 * Evaluates partial derivatives of lim_{p\to 0} p^3nu_s.
 */
real_t* SlowingDownFrequency::GetPartialP3NuSAtZero(len_t derivId){
    real_t preFactor = constPreFactor;
    len_t nMultiples = 1;
    if(derivId == id_ni)
        nMultiples = nzs;
    real_t *dP3nuS = new real_t[nr*nMultiples];
    for(len_t i = 0; i<nr*nMultiples; i++)
        dP3nuS[i] = 0;

    // Set partial n_cold 
    if(derivId == id_ncold)
        for(len_t ir=0; ir<nr; ir++){
            real_t lnLee0 = lnLambdaEE->evaluateAtP(ir,0);
            dP3nuS[ir] = preFactor * lnLee0 * evaluateElectronTermAtP(ir,0,collQtySettings->collfreq_mode);
        }
    else if(derivId == id_ni)    
        for(len_t ir=0; ir<nr; ir++){
            real_t electronTerm = preFactor*evaluateElectronTermAtP(ir,0,collQtySettings->collfreq_mode);
            real_t ntarget = unknowns->GetUnknownData(id_ncold)[ir];
            if (isNonScreened)
                ntarget += ionHandler->GetBoundElectronDensity(ir);
            for(len_t iz=0; iz<nZ; iz++)
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                    len_t indZ = ionIndex[iz][Z0];
                    dP3nuS[indZ*nr + ir] += ntarget*electronTerm*lnLambdaEE->evaluatePartialAtP(ir,0,id_ni,indZ);
                }
            if(isNonScreened)
                for(len_t iz=0; iz<nZ; iz++)
                    for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                        len_t indZ = ionIndex[iz][Z0];
                        dP3nuS[indZ*nr + ir] += (Zs[iz] - Z0) * electronTerm * lnLambdaEE->evaluateAtP(ir,0);
                    }
            else if(isPartiallyScreened)
                for(len_t iz=0; iz<nZ; iz++)
                    for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                        len_t indZ = ionIndex[iz][Z0];
                        dP3nuS[indZ*nr + ir] += preFactor * evaluateScreenedTermAtP(iz,Z0,0,collQtySettings->collfreq_mode);
                    }
        }
    else if(derivId == id_Tcold) 
        for(len_t ir=0; ir<nr; ir++){
            real_t ntarget = unknowns->GetUnknownData(id_ncold)[ir];
            if (isNonScreened)
                ntarget += ionHandler->GetBoundElectronDensity(ir);
           
            real_t lnLee0 = lnLambdaEE->evaluateAtP(ir,0);
            real_t dLnLee0 = lnLambdaEE->evaluatePartialAtP(ir,0,id_Tcold,0);
            dP3nuS[ir] = preFactor * ntarget * (lnLee0 * evaluateDDTElectronTermAtP(ir,0,collQtySettings->collfreq_mode)
                                                + dLnLee0 * evaluateElectronTermAtP(ir,0,collQtySettings->collfreq_mode));
        }
    return dP3nuS;
}


/**
 * Calculates a Rosenbluth potential matrix defined such that when it is muliplied
 * by the f_hot distribution vector, yields the slowing down frequency.
 */
void SlowingDownFrequency::calculateIsotropicNonlinearOperatorMatrix(){
    if( !(isPXiGrid && (mg->GetNp2() == 1)) )
        throw NotImplementedException("Nonlinear collisions only implemented for hot tails (np2=1) and p-xi grid");
    
    const real_t *p_f = mg->GetP1_f();
    const real_t *p = mg->GetP1();

    // See doc/notes/theory.pdf appendix B for details on discretization of integrals;
    // uses a trapezoidal rule
    real_t p2, p2f;
    real_t weightsIm1, weightsI;
    for (len_t i = 1; i<np1+1; i++){
        p2f = p_f[i]*p_f[i];
        p2  = p[0]*p[0];
        nonlinearMat[i][0] = 4*M_PI/p_f[i] * constPreFactor*( (p[1]-p[0])/2 + p[0]/3 )*p2/p2f;
        for (len_t ip = 1; ip < i-1; ip++){
            p2 = p[ip]*p[ip];
            nonlinearMat[i][ip] = 4*M_PI/p_f[i] * constPreFactor* trapzWeights[ip]*p2/p2f;
        } 
        p2 = p[i-1]*p[i-1];
        weightsIm1 = (p[i-1]-p[i-2])/2 + (p_f[i]-p[i-1])/(p[i]-p[i-1])*( (2*p[i]-p_f[i]-p[i-1])/2 );
        nonlinearMat[i][i-1] = 4*M_PI/p_f[i] * constPreFactor * weightsIm1*p2/p2f ;
        p2 = p[i]*p[i];
        weightsI = (p_f[i]-p[i-1])*(p_f[i]-p[i-1])/(p[i]-p[i-1]);
        nonlinearMat[i][i]   = 4*M_PI/p_f[i] * constPreFactor * (1.0/2)* weightsI *p2/p2f;
    }
}
