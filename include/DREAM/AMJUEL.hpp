#ifndef _DREAM_AMJUEL_HPP
#define _DREAM_AMJUEL_HPP

#include "DREAM/Constants.hpp"

namespace DREAM {
	class AMJUEL {
		private:
		static constexpr real_t recLyOpaque[81]={-2.959696621207E+01,-2.370057688281E-01,2.485234780243E-01,-9.938245216461E-02,1.980881578608E-02,-2.122222479009E-03,1.218203616198E-04,-3.464708585003E-06,3.763195232065E-08,
											-2.261509350573E+00,-3.916834765592E-01,4.175284638738E-01,-1.818480115491E-01,3.961587768727E-02,-4.770502751429E-03,3.183363649173E-04,-1.099361683574E-05,1.532076900830E-07,
											-4.674937331875E-01,-5.569001933269E-03,1.424684098594E-02,-7.573368055249E-03,2.345252431484E-03,-3.903412960562E-04,3.495619343912E-05,-1.549589929805E-06,2.649470795327E-08,
											2.507869795516E-01,2.497608379269E-02,-3.407068292654E-02,1.848151065349E-02,-4.874151697039E-03,6.859699832287E-04,-5.207997134901E-05,2.003342466410E-06,-3.056485646618E-08,
											2.069706780864E-02,-6.227904899439E-03,7.567752788769E-03,-3.519012590703E-03,7.900468003264E-04,-9.404610318098E-05,5.971364301743E-06,-1.910514915873E-07,2.423960414424E-09,
											-2.504106136665E-02,5.231970346733E-03,-4.961906824405E-03,1.530766046402E-03,-1.980846746944E-04,7.679462128899E-06,5.127742944183E-07,-5.042762997588E-08,1.117985320964E-09,
											4.740060719354E-03,-1.583429487117E-03,1.485117205295E-03,-4.566472727887E-04,6.072295908223E-05,-2.912042904116E-06,-6.966012349475E-08,1.046266371637E-08,-2.448127591728E-10,
											-3.716199599046E-04,1.790096302797E-04,-1.687816759412E-04,5.315093984335E-05,-7.528931894689E-06,4.568686753692E-07,-4.197741084618E-09,-6.390754230254E-10,1.836925499320E-11,
											1.078074419507E-05,-6.889626133438E-06,6.523611088216E-06,-2.087698447004E-06,3.070917297414E-07,-2.078805542480E-08,4.554850332229E-10,1.165730429588E-11,-4.891755053806E-13};
		
		static constexpr real_t ionizLyOpaque[81]={-2.842625123610E+01,3.816926440645E-02,-2.093090374769E-02,1.173278205549E-02,-3.195307082806E-03,4.891090254225E-04,-3.972967626237E-05,1.591035727928E-06,-2.476205052687E-08,
											1.212167851020E+01,-8.864661558973E-02,2.578404580790E-02,-1.125027025423E-02,3.132061119992E-03,-5.478419611750E-04,4.881508340119E-05,-2.055128753437E-06,3.275823973472E-08,
											-6.815821411657E+00,1.136458676986E-01,-3.209328253150E-02,7.428129132200E-03,-1.359196323415E-03,2.455552537690E-04,-2.356276122712E-05,1.024674585628E-06,-1.635599873170E-08,
											2.625844925126E+00,-7.176672848153E-02,2.221528694064E-02,-3.813499462891E-03,3.488231091499E-04,-5.136595197719E-05,5.236557713321E-06,-2.300633304095E-07,3.500516071428E-09,
											-6.666700835468E-01,2.362874407172E-02,-7.539637780701E-03,1.139599698461E-03,-5.231484976569E-05,3.851106631395E-06,-4.712865883416E-07,2.015344675266E-08,-2.228472576727E-10,
											1.063576010855E-01,-4.303004854934E-03,1.340877575599E-03,-1.830544200085E-04,3.418266263749E-06,2.695420674804E-07,-1.533579453796E-09,2.312045651251E-10,-2.651886075183E-11,
											-1.019791186281E-02,4.388098086177E-04,-1.289968813905E-04,1.595087521180E-05,-5.244682697928E-08,-4.069155577999E-08,6.001489992879E-10,-5.139558636314E-11,3.525452470540E-12,
											5.357498344762E-04,-2.353712235761E-05,6.392439908462E-06,-7.334538727835E-07,7.314879444968E-09,-1.321258886395E-09,3.249620872276E-10,-1.169248684183E-11,1.979747513777E-14,
											-1.183601163067E-05,5.183054638931E-07,-1.286865471654E-07,1.468328421917E-08,-8.136369296580E-10,2.036456697308E-10,-2.321096956876E-11,9.053581287584E-13,-1.015243451999E-14};
		
		static constexpr real_t ionizLossLyOpaque[81]={-2.431395592098E+01,-2.395941007384E-01,2.591565194903E-01,-1.087825788137E-01,2.304158174814E-02,-2.708007727528E-03,1.777247119352E-04,-6.080187141189E-06,8.438023009891E-08,
											1.113429718187E+01,1.849722545279E-01,-2.423728103974E-01,9.801393644159E-02,-1.946235514182E-02,2.079309382971E-03,-1.231134698172E-04,3.827948966411E-06,-4.895455818239E-08,
											-6.654446687338E+00,2.195366491981E-02,6.179888393756E-02,-2.979923286006E-02,5.913142010781E-03,-5.549143094821E-04,2.684044174965E-05,-6.633266274667E-07,6.895128780616E-09,
											2.747075059275E+00,-6.901300857989E-02,1.006521130909E-02,1.702193123225E-03,-6.298136767386E-04,4.093195127037E-05,-4.598542858756E-09,-5.542757391998E-08,9.032937294372E-10,
											-7.372137934626E-01,2.943069908340E-02,-7.829953302392E-03,6.752425702279E-04,2.856788083851E-05,-7.732295592583E-07,-3.124180764582E-07,1.188861911457E-08,1.649196327692E-11,
											1.227074461193E-01,-6.057176837025E-03,1.694081601397E-03,-1.621646872661E-04,-2.526799196487E-06,2.543500809367E-07,1.826904076362E-08,4.273639168838E-10,-5.875558984747E-11,
											-1.216484938923E-02,6.701189895777E-04,-1.793079853259E-04,1.506237912301E-05,4.384573257407E-07,-1.137888186554E-08,-4.163110915594E-09,3.345977887826E-11,5.495071058936E-12,
											6.566484532457E-04,-3.779995638455E-05,8.950721046575E-06,-4.568263582843E-07,-5.697388918232E-08,-1.072747527249E-09,6.106046869680E-10,-2.063138694762E-11,-9.709542597057E-16,
											-1.483802107723E-05,8.441118380444E-07,-1.530047985372E-07,-6.560997301392E-09,3.141908302488E-09,3.111434993423E-12,-2.241669215814E-11,9.998773490453E-13,-9.539526341658E-15};

		static constexpr real_t recRadLyOpaque[81]={-2.626461971500E+01,-1.141522828006E-01,1.076393602102E-01,-3.951176865624E-02,7.434448907130E-03,-7.651735363726E-04,4.256550578560E-05,-1.183957210565E-06,1.267430439071E-08,
											-8.898849653304E-02,-3.942292858703E-01,4.197250110873E-01,-1.819873891755E-01,3.954265794034E-02,-4.746673604459E-03,3.160257175597E-04,-1.090132407225E-05,1.518762884136E-07,
											-4.795279065913E-01,-1.335489318623E-01,1.698162748512E-01,-7.777618558640E-02,1.791770402443E-02,-2.265419802670E-03,1.597834560300E-04,-5.867107114568E-06,8.703048627852E-08,
											-6.457641001473E-02,4.613800494941E-02,-6.378456442460E-02,3.368756666920E-02,-8.588159967647E-03,1.161221464744E-03,-8.476482397035E-05,3.148378986370E-06,-4.661383263216E-08,
											6.800392305050E-02,2.637589098726E-02,-2.927766641228E-02,1.215210455092E-02,-2.536449841814E-03,2.949052179931E-04,-1.945865547958E-05,6.794359168371E-07,-9.714166652574E-09,
											-7.780596827160E-03,-1.696354616404E-03,4.101919602398E-03,-2.853313172787E-03,8.331295381374E-04,-1.219354343653E-04,9.388539733773E-06,-3.625759741848E-07,5.530677587059E-09,
											-9.252420142124E-04,-2.153234057341E-03,1.724087525894E-03,-3.889280715498E-04,1.445906496031E-05,5.400239087927E-06,-7.463824586159E-07,3.657036916983E-08,-6.329837178101E-10,
											2.115742192807E-04,3.847511359996E-04,-3.531396348434E-04,1.117507761644E-04,-1.626382515941E-05,1.138105849294E-06,-3.277550512167E-08,-2.762994890889E-11,1.306756395427E-11,
											-9.909336050813E-06,-1.829769520002E-05,1.745121101762E-05,-5.935692251130E-06,9.750922276399E-07,-8.514634639125E-08,4.007910080852E-09,-9.447800370353E-11,8.488148236190E-13};
		public:
		AMJUEL(){}
		~AMJUEL(){}
        real_t getRecLyOpaque(len_t Z0, real_t n, real_t T){
        	if(Z0==1){
				real_t lnRecAMJUEL=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++)
						lnRecAMJUEL+=recLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
				return exp(lnRecAMJUEL)/1e6;
			}else{
				return 0;
			}        			
        }
        
        real_t getRecLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
        	if(Z0==1){
				real_t lnRecAMJUEL=0;
				real_t derivFactor=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++){
						lnRecAMJUEL+=recLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
						if(j>0)
							derivFactor+=recLyOpaque[i*9+j]*pow(log(T),i)*j*pow(log(n/1e14),j-1)/n;
					}
				return derivFactor*exp(lnRecAMJUEL)/1e6;
			}else{
				return 0;
			}        			
        }
        
        real_t getRecLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
        	if(Z0==1){
				real_t lnRecAMJUEL=0;
				real_t derivFactor=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++){
						lnRecAMJUEL+=recLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
						if(i>0)
							derivFactor+=recLyOpaque[i*9+j]*i*pow(log(T),i-1)/T*pow(log(n/1e14),j);
					}
				return derivFactor*exp(lnRecAMJUEL)/1e6;
			}else{
				return 0;
			}        			
        }
        
        real_t getRecRadLyOpaque(len_t Z0, real_t n, real_t T){
        	if(Z0==1){
				real_t lnRecRadAMJUEL=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++)
						lnRecRadAMJUEL+=recRadLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
				return Constants::ec*exp(lnRecRadAMJUEL)/1e6;
			}else{
				return 0;
			}        			
        }
        
        real_t getRecRadLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
        	if(Z0==1){
				real_t lnRecRadAMJUEL=0;
				real_t derivFactor=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++){
						lnRecRadAMJUEL+=recRadLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
						if(j>0)
							derivFactor+=recRadLyOpaque[i*9+j]*pow(log(T),i)*j*pow(log(n/1e14),j-1)/n;
					}
				return Constants::ec*derivFactor*exp(lnRecRadAMJUEL)/1e6;
			}else{
				return 0;
			}        			
        }
        
        real_t getRecRadLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
        	if(Z0==1){
				real_t lnRecRadAMJUEL=0;
				real_t derivFactor=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++){
						lnRecRadAMJUEL+=recRadLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
						if(i>0)
							derivFactor+=recRadLyOpaque[i*9+j]*i*pow(log(T),i-1)/T*pow(log(n/1e14),j);
					}
				return Constants::ec*derivFactor*exp(lnRecRadAMJUEL)/1e6;
			}else{
				return 0;
			}        			
        }
        
        real_t getIonizLyOpaque(len_t Z0, real_t n, real_t T){
        	if(Z0==0){
				real_t lnIonizAMJUEL=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++)
						lnIonizAMJUEL+=ionizLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
				return exp(lnIonizAMJUEL)/1e6;
			}else{
				return 0;
			}			
        }
        
        real_t getIonizLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
        	if(Z0==0){
				real_t lnIonizAMJUEL=0;
				real_t derivFactor=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++){
						lnIonizAMJUEL+=ionizLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
						if(j>0)
							derivFactor+=ionizLyOpaque[i*9+j]*pow(log(T),i)*j*pow(log(n/1e14),j-1)/n;
					}
				return derivFactor*exp(lnIonizAMJUEL)/1e6;
			}else{
				return 0;
			}			
        }
        
        real_t getIonizLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
        	if(Z0==0){
				real_t lnIonizAMJUEL=0;
				real_t derivFactor=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++){
						lnIonizAMJUEL+=ionizLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
						if(i>0)
							derivFactor+=ionizLyOpaque[i*9+j]*i*pow(log(T),i-1)/T*pow(log(n/1e14),j);
					}
				return derivFactor*exp(lnIonizAMJUEL)/1e6;
			}else{
				return 0;
			}			
        }
        
        real_t getIonizLossLyOpaque(len_t Z0, real_t n, real_t T){
        	if(Z0==0){
				real_t lnIonizLossAMJUEL=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++)
						lnIonizLossAMJUEL+=ionizLossLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
				return Constants::ec*exp(lnIonizLossAMJUEL)/1e6;
			}else{
				return 0;
			}        			
        }
        
        real_t getIonizLossLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
        	if(Z0==0){
				real_t lnIonizLossAMJUEL=0;
				real_t derivFactor=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++){
						lnIonizLossAMJUEL+=ionizLossLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
						if(j>0)
							derivFactor+=ionizLossLyOpaque[i*9+j]*pow(log(T),i)*j*pow(log(n/1e14),j-1)/n;
					}
				return Constants::ec*derivFactor*exp(lnIonizLossAMJUEL)/1e6;
			}else{
				return 0;
			}        			
        }
        
        real_t getIonizLossLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
        	if(Z0==0){
				real_t lnIonizLossAMJUEL=0;
				real_t derivFactor=0;
				for(len_t i=0;i<9;i++)
					for(len_t j=0;j<9;j++){
						lnIonizLossAMJUEL+=ionizLossLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
						if(i>0)
							derivFactor+=ionizLossLyOpaque[i*9+j]*i*pow(log(T),i-1)/T*pow(log(n/1e14),j);
					}
				return Constants::ec*derivFactor*exp(lnIonizLossAMJUEL)/1e6;
			}else{
				return 0;
			}        			
        }
    };
}
#endif/*_DREAM_EQUATION_AMJUEL_HPP*/
        
