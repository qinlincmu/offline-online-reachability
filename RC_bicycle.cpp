#include "../flowstar/Continuous.h"
#include <tuple>
#include <unordered_map>


using namespace flowstar;
using namespace std;
#include <iostream>
#include <fstream>
#include <boost/numeric/interval.hpp>
#include <gmp.h>
#include <mpfr.h>
#define M_PI 3.1415
#define LR 1.364
#define LF 1.364

struct full_rounding:
  boost::numeric::interval_lib::rounded_arith_opp<double>
{
private:
  typedef int mpfr_func(mpfr_t, const __mpfr_struct*, mp_rnd_t);
  double invoke_mpfr(double x, mpfr_func f, mp_rnd_t r) {
    mpfr_t xx;
    mpfr_init_set_d(xx, x, GMP_RNDN);
    f(xx, xx, r);
    double res = mpfr_get_d(xx, r);
    mpfr_clear(xx);
    return res;
  }
public:
# define GENR_FUNC(name) \
  double name##_down(double x) { return invoke_mpfr(x, mpfr_##name, GMP_RNDD); } \
  double name##_up  (double x) { return invoke_mpfr(x, mpfr_##name, GMP_RNDU); }
  GENR_FUNC(exp)
  GENR_FUNC(log)
  GENR_FUNC(sin)
  GENR_FUNC(cos)
  GENR_FUNC(tan)
  GENR_FUNC(asin)
  GENR_FUNC(acos)
  GENR_FUNC(atan)
  GENR_FUNC(sinh)
  GENR_FUNC(cosh)
  GENR_FUNC(tanh)
  GENR_FUNC(asinh)
  GENR_FUNC(acosh)
  GENR_FUNC(atanh)
};

namespace dummy {
  using namespace boost;
  using namespace numeric;
  using namespace interval_lib;
  typedef save_state<full_rounding> R;
  typedef checking_strict<double> P;
  typedef interval<double, policies<R, P> > IT;
};

typedef dummy::IT IT;


struct Point
{
	float x;
	float y;
};

tuple<double, double> EstimateBeta(float delta_low, float delta_up){
	IT delta(delta_low, delta_up);
	float beta_low;
	float beta_up;
	IT result = atanh((LR/(LF+LR))*tan(delta));
	beta_low = result.lower();
	beta_up = result.upper();
	return make_tuple(beta_low, beta_up);
}

int overlap(Point l1, Point r1, Point l2, Point r2){
    if(l1.x >= r2.x or l2.x >= r1.x)
        return 0;
    if(l1.y <= r2.y or l2.y <= r1.y)
        return 0;
    return 1;
}

														   

tuple<vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>> GetPipesPerCycle(double x_0_min, double x_0_max, double y_0_min,
                                                               double y_0_max, double psi_0_min,  double psi_0_max, 
                                                               double delta_0_min, double delta_0_max, double v_0_min, double v_0_max, double horizon, double stepsize)
{
	unsigned int numVars = 7;

	std::vector<double> x_tube;
	std::vector<double> y_tube;
	std::vector<double> yaw_tube;
	std::vector<double> vx_tube;
	std::vector<double> vy_tube;
	std::vector<double> yawdot_tube;
    
	Variables _stateVars;
	stateVars = _stateVars;
	int x_id = stateVars.declareVar("x");
	int y_id = stateVars.declareVar("y");
	int psi_id = stateVars.declareVar("psi");
	int delta_id = stateVars.declareVar("delta");
	int v_id = stateVars.declareVar("v");
	int beta_id = stateVars.declareVar("beta");
	int t_id = stateVars.declareVar("t");





	// define the dynamics
	Expression_AST<Real> ode_expression_x("v * cos(psi+beta)");
	Expression_AST<Real> ode_expression_y("v * sin(psi+beta)");
	Expression_AST<Real> ode_expression_psi("v * sin(beta)/"+std::to_string(LR));
	Expression_AST<Real> ode_expression_delta("0");
	Expression_AST<Real> ode_expression_beta("0");
	Expression_AST<Real> ode_expression_v("0");
	Expression_AST<Real> ode_expression_t("1");



	vector<Expression_AST<Real> > ode_rhs(numVars);
	ode_rhs[x_id] = ode_expression_x;
	ode_rhs[y_id] = ode_expression_y;
	ode_rhs[psi_id] = ode_expression_psi;
	ode_rhs[delta_id] = ode_expression_delta;
	ode_rhs[v_id] = ode_expression_v;
	ode_rhs[beta_id] = ode_expression_beta;
	ode_rhs[t_id] = ode_expression_t;



	Deterministic_Continuous_Dynamics dynamics(ode_rhs);



	// set the reachability parameters
	Computational_Setting setting;

	// set the stepsize and the order
	setting.setFixedStepsize(stepsize, 5);

	// set the time horizon
	setting.setTime(horizon);

	// set the cutoff threshold
	setting.setCutoffThreshold(1e-8);

	// set the queue size for the symbolic remainder, it is 0 if symbolic remainder is not used
	setting.setQueueSize(100);

	// print out the computation steps
	setting.printOff();

	// set up the remainder estimation
	Interval I(-0.01, 0.01);
	vector<Interval> remainder_estimation(numVars, I);
	setting.setRemainderEstimation(remainder_estimation);

	// call this function when all of the parameters are defined
	setting.prepare();

    //tuple<vector<double>, vector<double>> beta;
	float beta_0_min;
	float beta_0_max;
	tie(beta_0_min, beta_0_max) = EstimateBeta(delta_0_min, delta_0_max);


	// define the initial set which is a box



	Interval init_x(x_0_min, x_0_max), init_y(y_0_min, y_0_max), init_psi(psi_0_min, psi_0_max),
			init_delta(delta_0_min, delta_0_max), init_v(v_0_min, v_0_max), init_beta(beta_0_min, beta_0_max), init_t;

	vector<Interval> initial_box(numVars);
	initial_box[x_id] = init_x;
	initial_box[y_id] = init_y;
	initial_box[psi_id] = init_psi;
	initial_box[delta_id] = init_delta;
	initial_box[v_id] = init_v;
	initial_box[beta_id] = init_beta;
	initial_box[t_id] = init_t;

	Flowpipe initialSet(initial_box);


	// empty unsafe set
	vector<Constraint> unsafeSet;

	/*
	 * The structure of the class Result_of_Reachability is defined as below:
	 * nonlinear_flowpipes: the list of computed flowpipes
	 * tmv_flowpipes: translation of the flowpipes, they will be used for further analysis
	 * fp_end_of_time: the flowpipe at the time T
	 */
	Result_of_Reachability result;

	// run the reachability computation
	clock_t begin, end;
	//begin = clock();

	dynamics.reach(result, setting, initialSet, unsafeSet);
//	dynamics.reach(result, setting, result.fp_end_of_time, unsafeSet);
//	dynamics.reach(result, setting, result.fp_end_of_time, unsafeSet);

	//end = clock();
	//printf("Inner time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	// flowpipes should be translated to single Taylor model vectors before plotting
	result.transformToTaylorModels(setting);

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result.nonlinear_flowpipes.begin();
	std::vector<unsigned int> varIDs({0, 1, 2});
	

	IT beta(beta_0_min, beta_0_max);
	IT v(v_0_min, v_0_max);
	float beta_low;
	float beta_up;
	IT result_vx = v*cos(beta);
	IT result_vy = v*sin(beta);
	IT result_yawdot = v*sin(beta)/LR;

	for(; tmvIter != result.tmv_flowpipes.end() ; ++tmvIter, ++fpIter)
	{
		std::vector<Interval> box;
		tmvIter->intEval(box, fpIter->domain, varIDs);

		Interval X = box[0], Y = box[1], PSI = box[2];
		x_tube.push_back(X.inf());
		x_tube.push_back(X.sup());

		y_tube.push_back(Y.inf());
		y_tube.push_back(Y.sup());

		yaw_tube.push_back(PSI.inf());
		yaw_tube.push_back(PSI.sup());

		vx_tube.push_back(result_vx.lower());
		vx_tube.push_back(result_vx.upper());

		vy_tube.push_back(result_vy.lower());
		vy_tube.push_back(result_vy.upper());

		yawdot_tube.push_back(result_yawdot.lower());
		yawdot_tube.push_back(result_yawdot.upper());

	}

	return make_tuple(x_tube, y_tube, yaw_tube, vx_tube, vy_tube, yawdot_tube);
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>> GetPipesMultiCycle(double x_0_min, double x_0_max, double y_0_min,
                                                               double y_0_max, double psi_0_min,  double psi_0_max, 
                                                               vector<double> delta, vector<double> v, double horizon, double stepsize)
{
    std::vector<double> x_tube;
    std::vector<double> y_tube;
    std::vector<double> psi_tube;
	std::vector<double> vx_tube;
	std::vector<double> vy_tube;
	std::vector<double> psidot_tube;

	std::vector<double> x_tube_i; //tube in one control cycle
    std::vector<double> y_tube_i;
    std::vector<double> psi_tube_i;
	std::vector<double> vx_tube_i;
	std::vector<double> vy_tube_i;
	std::vector<double> psidot_tube_i;

	for(int i=0; i<delta.size(); i++){
        tie(x_tube_i, y_tube_i, psi_tube_i, vx_tube_i, vy_tube_i, psidot_tube_i) = GetPipesPerCycle(x_0_min, x_0_max, y_0_min,
                                                               y_0_max, psi_0_min, psi_0_max, 
                                                               delta.at(i), delta.at(i), v.at(i), v.at(i), horizon, stepsize);
		int size_per_cycle = x_tube_i.size(); //number of boxes in one cycle

		//extract the last box 
		x_0_min = x_tube_i.at(size_per_cycle-2);
		x_0_max = x_tube_i.at(size_per_cycle-1);
		y_0_min = y_tube_i.at(size_per_cycle-2);
		y_0_max = y_tube_i.at(size_per_cycle-1);
		psi_0_min = psi_tube_i.at(size_per_cycle-2);
		psi_0_max = psi_tube_i.at(size_per_cycle-1);
		float vx_min = vx_tube_i.at(size_per_cycle-2);
		float vx_max = vx_tube_i.at(size_per_cycle-1);
		float vy_min = vy_tube_i.at(size_per_cycle-2);
		float vy_max = vy_tube_i.at(size_per_cycle-1);
		float psidot_min = psidot_tube_i.at(size_per_cycle-2);
		float psidot_max = psidot_tube_i.at(size_per_cycle-1);
		
		x_tube.push_back(x_0_min);
		x_tube.push_back(x_0_max);
		y_tube.push_back(y_0_min);
		y_tube.push_back(y_0_max);
		psi_tube.push_back(psi_0_min);
		psi_tube.push_back(psi_0_max);
		vx_tube.push_back(vx_min);
		vx_tube.push_back(vx_max);
		vy_tube.push_back(vy_min);
		vy_tube.push_back(vy_max);
		psidot_tube.push_back(psidot_min);
		psidot_tube.push_back(psidot_max);


		x_tube_i.clear();
		y_tube_i.clear();
		psi_tube_i.clear();
		vx_tube_i.clear();
		vy_tube_i.clear();
		psidot_tube_i.clear();
	}
	return make_tuple(x_tube, y_tube, psi_tube, vx_tube, vy_tube, psidot_tube);
}	

std::vector<float> rangeByStep(float start, float end, float step){
	/*
	Function: generate interval with even step 
	*/
    std::vector<float> range;
	float endi = start;
    while(endi <= end){
		range.push_back(endi);
		endi += step;
	}
	return range;
}

void offline_tube_generate(std::vector<float> delta){
	/*
    Function: generate tubes offline
    Input: 
        delta: control sampling space   
        Output: store lookup table: yaw_min, yaw_max, steering angle, x1_low, x1_high...x3_low, x3_high, same for y, yaw, vx, vy, yawdot
	*/
	std::vector<float> psi; //yaw range
    psi = rangeByStep(-1*M_PI/2, 1*M_PI/2, 0.1);
	cout<<"yaw-size"<<psi.size()<<endl;
	float v = 10;
	float x_uncertainty = 0.05;
	float y_uncertainty = 0.05;
	float psi_uncertainty = 0.05;
    float x0 = 0;
	float y0 = 0;
	int cycle = 3; //number of control cycles
	float step_size = 0.005; //stepsize in one control cycle
	float horizon_per_cycle = 0.025; //horizon in one control cycle, control sampling time
	std::vector<double> x_tube;
	std::vector<double> y_tube;
	std::vector<double> psi_tube;
	std::vector<double> vx_tube;
	std::vector<double> vy_tube;
	std::vector<double> psidot_tube;
	ofstream lookup_tablefile;
	std::vector<double> delta_queue; //sequence of steering angle, length=cycle
	std::vector<double> v_queue; ////sequence of velocity command, length=cycle
    lookup_tablefile.open ("lookup.txt");
    cout<<"delta-size"<<delta.size()<<endl;
    for(unsigned i=0; i<psi.size(); i++){
		for(unsigned j=0; j<delta.size(); j++){
            for(unsigned ci=0; ci<cycle; ci++){
                delta_queue.push_back(delta.at(j));
				v_queue.push_back(v);
			}
            tie(x_tube, y_tube, psi_tube, vx_tube, vy_tube, psidot_tube) = GetPipesMultiCycle(x0-x_uncertainty, x0+x_uncertainty, y0-y_uncertainty, y0+y_uncertainty,
			                                psi.at(i)-psi_uncertainty, psi.at(i)+psi_uncertainty, delta_queue, v_queue, horizon_per_cycle, step_size);
		    lookup_tablefile << std::to_string(psi.at(i)-0.05)+" "+std::to_string(psi.at(i)+0.05)+" "+std::to_string(delta.at(j))+" ";
			for(unsigned k=0; k<x_tube.size(); k++){
			    lookup_tablefile << std::to_string(x_tube.at(k))+" ";
			}
			for(unsigned k=0; k<y_tube.size(); k++){
                lookup_tablefile << std::to_string(y_tube.at(k))+" ";
			}
			for(unsigned k=0; k<psi_tube.size(); k++){
                lookup_tablefile << std::to_string(psi_tube.at(k))+" ";
			}
			for(unsigned k=0; k<vx_tube.size(); k++){
                lookup_tablefile << std::to_string(vx_tube.at(k))+" ";
			}
			for(unsigned k=0; k<vy_tube.size(); k++){
                lookup_tablefile << std::to_string(vy_tube.at(k))+" ";
			}
			for(unsigned k=0; k<psidot_tube.size(); k++){
				if(k<psidot_tube.size()-1)
                    lookup_tablefile << std::to_string(psidot_tube.at(k))+" ";
				else
                    lookup_tablefile << std::to_string(psidot_tube.at(k));
			}
            lookup_tablefile << "\n";
			delta_queue.clear();
			v_queue.clear();
		}
	}
	lookup_tablefile.close();
}


unordered_map<string, vector<float>> load_table(){
	/*
	Function: reading lookup table file and store it as an unordered map
	          key of the map (yaw_low_bound, yaw_up_bound, delta)
			  value vector of tube: x1_low, x1_up, x2_low,x2_up ... y1_low, y1_up ...
	*/
    char buf[1000];
    FILE *fp;
    int len;
    const char sep[2] = " ";
    char *token;
    std::vector<string> temp;
	std::vector<float> data;
    if((fp = fopen("lookup.txt", "r")) == NULL){
      perror("fail to read");
      exit(1);
    }
	unordered_map<string, vector<float>> map;
	std::string key;
    while (fgets(buf, 1000, fp) != NULL) {

        len = strlen(buf);
        buf[len-1] = '\0';
        token = strtok(buf, sep);
        while (token != NULL){
          temp.push_back(token);
          token = strtok(NULL, sep);}
        key = temp[0]+" "+temp[1]+" "+temp[2];
		int tubesize = temp.size()-3;
		for (int i=3; i<tubesize; i++){
			data.push_back(std::stof(temp[i]));
		}
        map[key] = data;
        temp.clear();
        data.clear();
		}
	return map;
}


vector<float> split_string(string input, char c){
	std::istringstream ss(input);
	std::string token;
	vector<float> split_vector;
	while(std::getline(ss, token, c)){
		split_vector.push_back(stof(token));
	}
	return split_vector;
}

tuple<float, float> get_yaw_bound(float yaw, unordered_map<string, vector<float>> map){
	/*
	Function: given yaw value, find its low and up bound in the lookuptable
	*/
	float yaw_min=-1000;
	float yaw_max=1000;
	vector<float> key_vector;
	for(auto e : map) {
	    string key = e.first;
        key_vector = split_string(key, ' ');
		if ((yaw >= key_vector[0]) && (yaw <= key_vector[1])){
			yaw_min = key_vector[0];
			yaw_max = key_vector[1];
			break;
		}
	}
	return make_tuple(yaw_min, yaw_max);
}

int check_safety(vector<float> x, vector<float> y, float xp, float yp, Point l2, Point r2){
	/*
    Function: given xy tube and xy position, check its overlap with obstacle
	       x,y: xy tube
		   xp,yp: xy position
           l2: left upper point of obstacle represented by rectangle
           r2: right lower point of obstacle represented by rectangle  
    Output: flag, safe-1; unsafe-0
	*/
	int flag = 1;
	for(int j=0; j<x.size(); j++){
	    Point l1;
		Point r1;
		l1.x = x[j*2]+xp;
		l1.y = y[j*2+1]+yp;
		r1.x = x[j*2+1]+xp;
		r1.y = y[j*2]+yp;
		if(overlap(l1, r1, l2, r2) == 1){
		    flag = 0;
			return flag;
		}
		else
		    continue;
		}
	return flag;
}

void online_check_for_single_grid(unordered_map<string, vector<float>> map, float xp, float yp, float yaw, std::vector<float> delta, Point l2, Point r2){
	/*
    Function: for a single point, load the offline tubes and check their overlap with obstacle,
    if multiple tubes found, we don't need to take over; if only one tube is safe, must take over 
    Input: map: lookup table
	       x,y,yaw of the point
           delta: control sampling space
           l2: left upper point of obstacle represented by rectangle
           r2: right lower point of obstacle represented by rectangle  
    Output: None
	*/
    float yaw_min;
	float yaw_max;
	vector<float> pos;
	int key_n=0;
	tie(yaw_min, yaw_max) = get_yaw_bound(yaw, map);

	std::vector<float> data;
	int flag;
	string keyi;
	std::vector<float> xd;
	std::vector<float> yd;
    for(int i=0; i<delta.size(); i++){
		xd.clear();
		yd.clear();
		flag = 1;
		keyi = to_string(yaw_min)+" "+to_string(yaw_max)+" "+to_string(delta[i]);
        if(map.find(keyi) != map.end()){
            std::cout << "Element Found "<< std::endl;
			data = map.at(keyi);
	        for(int di=0; di<data.size()/2; di++){
	        xd.push_back(data.at(di));
		    yd.push_back(data.at(data.size()/2+di));
	        }
		    flag = check_safety(xd, yd, xp, yp, l2, r2);
        }
        else
        {
            std::cout << "Element Not Found " << delta[i] << std::endl;
		 	continue;
        }
		if(flag == 1)
		    pos.push_back(delta[i]);
	}
	if(pos.size() <= 1)
	    std::cout<<"Must take over to execute fallback"<<"\n";
    return;
}

void online_check_for_interval(unordered_map<string, vector<float>> map, float x_min, float x_max, float y_min, float y_max, float yaw_min, float yaw_max, std::vector<float> delta, Point l2, Point r2){
	/*
	 Function: for x,y,yaw intervals, load the offline tubes and check each cell's overlap with obstacle,
        if multiple tubes found, we don't need to take over; if only one tube is safe, must take over 
        Input: map: lookup table
		       x,y,yaw intervals
               delta: control sampling space
               l2: left upper point of obstacle represented by rectangle
               r2: right lower point of obstacle represented by rectangle  
        Output: None, can be extended based on how you want to use it
	*/
	std::vector<float> x_r;
	std::vector<float> y_r;
	std::vector<float> yaw_r;
	std::vector<float> pos;
	std::vector<float> xd; //x coordinates of a tube 
	std::vector<float> yd; //y coordinates of a tube
	std::vector<float> data; //tube data
	float yaw_left;
	float yaw_right;
	float xp;
	float yp;
	string key;
	int flag;


    x_r = rangeByStep(x_min, x_max, 0.1);
    y_r = rangeByStep(y_min, y_max, 0.1);
    yaw_r = rangeByStep(yaw_min, yaw_max, 0.1);
    cout<<yaw_r.size()<<" "<<delta.size()<<" "<<x_r.size()<<" "<<y_r.size()<<endl;
    for(int i=0; i<yaw_r.size(); i++){
		tie(yaw_left, yaw_right) = get_yaw_bound(yaw_r.at(i), map);
		for(int j=0; j<delta.size(); j++){
			key = to_string(yaw_left)+" "+to_string(yaw_right)+" "+to_string(delta[j]);
			xd.clear();
			yd.clear();
            if(map.find(key) != map.end()){
                data = map.at(key);
	            for(int di=0; di<data.size()/2; di++){
	            xd.push_back(data.at(di));
		        yd.push_back(data.at(data.size()/2+di));
	            }
			for(int l=0; l<x_r.size(); l++){
				for(int m=0; m<y_r.size(); m++){
					xp = x_r.at(l);
					yp = y_r.at(m);
				    flag = check_safety(xd, yd, xp, yp, l2, r2);
				}
			}
			}
			else{
		 	    continue;
            }
        }
	}
	return;
}
int main(){
	clock_t begin, end;

	std::vector<float> delta;

	delta = rangeByStep(-(1.0/4)*M_PI, (1.0/4)*M_PI, 0.1);
	begin = clock();
	offline_tube_generate(delta);
	end = clock();
	printf("generating lookup table time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	// unordered_map<string, vector<float>> map;
	// map = load_table();
    
	// float xp = 1;
    // float yp = 2;
    // float yaw = 0.3;
	// Point l2;
	// Point r2;
    // l2.x = -1;
	// l2.y = 1;
    // r2.x = 1;
	// r2.y = -1;
	// begin = clock();
	// online_check_for_single_grid(map, xp, yp, yaw, delta, l2, r2);
	// end = clock();
	// printf("single_point check time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

    // float x_min = -0.5;
    // float x_max = 0.5;
    // float y_min = -0.5;
    // float y_max = 0.5;
    // float yaw_min = -M_PI/2.0;
    // float yaw_max = M_PI/2.0;
    
	// begin = clock();
    // online_check_for_interval(map, x_min, x_max, y_min, y_max, yaw_min, y_max, delta, l2, r2);
	// end = clock();
	// printf("interval check time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	return 0;
}
	
