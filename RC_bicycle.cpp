#include "../flowstar/Continuous.h"
#include <tuple>
#include <unordered_map>


using namespace flowstar;
using namespace std;
#include <iostream>
#include <fstream>
#define M_PI 3.1415


struct Point
{
	float x;
	float y;
};


int overlap(Point l1, Point r1, Point l2, Point r2){
    if(l1.x >= r2.x or l2.x >= r1.x)
        return 0;
    if(l1.y <= r2.y or l2.y <= r1.y)
        return 0;
    return 1;
}


tuple<vector<double>, vector<double>, vector<double>> getPipes(double x_0_min, double x_0_max, double y_0_min, double y_0_max, double psi_0_min,  double psi_0_max,  double delta_0_min, double delta_0_max, double v_0_min, double v_0_max, double horizon, double stepsize, std::string L)
{
	unsigned int numVars = 6;

	std::vector<double> x_tube;
	std::vector<double> y_tube;
	std::vector<double> yaw_tube;
    
	Variables _stateVars;
	stateVars = _stateVars;
	int x_id = stateVars.declareVar("x");
	int y_id = stateVars.declareVar("y");
	int psi_id = stateVars.declareVar("psi");
	int delta_id = stateVars.declareVar("delta");
	int v_id = stateVars.declareVar("v");
	int t_id = stateVars.declareVar("t");





	// define the dynamics
	Expression_AST<Real> ode_expression_x("v * cos(psi)");
	Expression_AST<Real> ode_expression_y("v * sin(psi)");
	Expression_AST<Real> ode_expression_psi("v * sin(delta)/"+L+"*cos(delta)");
	//Expression_AST<Real> ode_expression_psi("v * delta/"+L);
	//Expression_AST<Real> ode_expression_psi("0");
	Expression_AST<Real> ode_expression_delta("0");
	Expression_AST<Real> ode_expression_v("0");
	Expression_AST<Real> ode_expression_t("1");



	vector<Expression_AST<Real> > ode_rhs(numVars);
	ode_rhs[x_id] = ode_expression_x;
	ode_rhs[y_id] = ode_expression_y;
	ode_rhs[psi_id] = ode_expression_psi;
	ode_rhs[delta_id] = ode_expression_delta;
	ode_rhs[v_id] = ode_expression_v;
	ode_rhs[t_id] = ode_expression_t;



	Deterministic_Continuous_Dynamics dynamics(ode_rhs);



	// set the reachability parameters
	Computational_Setting setting;

	// set the stepsize and the order
	setting.setFixedStepsize(stepsize, 3);

	// set the time horizon
	setting.setTime(horizon);

	// set the cutoff threshold
	setting.setCutoffThreshold(1e-8);

	// set the queue size for the symbolic remainder, it is 0 if symbolic remainder is not used
	setting.setQueueSize(0);

	// print out the computation steps
	setting.printOff();

	// set up the remainder estimation
	Interval I(-0.01, 0.01);
	vector<Interval> remainder_estimation(numVars, I);
	setting.setRemainderEstimation(remainder_estimation);

	// call this function when all of the parameters are defined
	setting.prepare();


	// define the initial set which is a box
	Interval init_x(x_0_min, x_0_max), init_y(y_0_min, y_0_max), init_psi(psi_0_min, psi_0_max),
			init_delta(delta_0_min, delta_0_max), init_v(v_0_min, v_0_max), init_t;

	vector<Interval> initial_box(numVars);
	initial_box[x_id] = init_x;
	initial_box[y_id] = init_y;
	initial_box[psi_id] = init_psi;
	initial_box[delta_id] = init_delta;
	initial_box[v_id] = init_v;
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
                

	//Plot_Setting plot_setting;
	//plot_setting.printOff();
	//plot_setting.setOutputDims(x_id, y_id);
        //int indicator = 1;
	//int indicator = plot_setting.plot_2D_interval_MATLAB("RC_bicycle", result);//Qin
        //plot_setting.plot_2D_interval_MATLAB("RC_bicycle", result);

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result.nonlinear_flowpipes.begin();
	std::vector<unsigned int> varIDs({0, 1, 2});
	
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

	}

	return make_tuple(x_tube, y_tube, yaw_tube);
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
        Output: store lookup table
	*/
	std::vector<float> yaw;
    yaw = rangeByStep(-1*M_PI/2, 1*M_PI/2, 0.1);
	cout<<"yaw-size"<<yaw.size()<<endl;
	float v = 10;
	float x_uncertainty = 0.05;
	float y_uncertainty = 0.05;
	float yaw_uncertainty = 0.05;
    float x0 = 0;
	float y0 = 0;
	std::vector<double> x_tube;
	std::vector<double> y_tube;
	std::vector<double> yaw_tube;
	ofstream lookup_tablefile;
    lookup_tablefile.open ("lookup.txt");
    cout<<"delta-size"<<delta.size()<<endl;
    for(unsigned i=0; i<yaw.size(); i++){
		for(unsigned j=0; j<delta.size(); j++){
            tie(x_tube, y_tube, yaw_tube) = getPipes(x0-x_uncertainty, x0+x_uncertainty, y0-y_uncertainty, y0+y_uncertainty, yaw.at(i)-yaw_uncertainty, yaw.at(i)+yaw_uncertainty, delta.at(j), delta.at(j), v, v, 3*0.025, 0.025, "2.728");
		    lookup_tablefile << std::to_string(yaw.at(i)-0.05)+" "+std::to_string(yaw.at(i)+0.05)+" "+std::to_string(delta.at(j))+" ";
			for(unsigned k=0; k<x_tube.size(); k++){
			    lookup_tablefile << std::to_string(x_tube.at(k))+" ";
			}
			for(unsigned k=0; k<y_tube.size(); k++){
				if(k<y_tube.size()-1)
			        lookup_tablefile << std::to_string(y_tube.at(k))+" ";
				else
				    lookup_tablefile << std::to_string(y_tube.at(k));	
			}
            lookup_tablefile << "\n";
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

	std::vector<float> x;
	std::vector<float> y;
	std::vector<float> data;
	int flag;
	string keyi;
	std::vector<float> xd;
	std::vector<float> yd;
    for(int i=0; i<delta.size(); i++){
		x.clear();
		y.clear();
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
	std::vector<float> x;
	std::vector<float> y;
	std::vector<float> data;
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
            if(map.find(key) != map.end()){
                data = map.at(key);
	            for(int di=0; di<data.size()/2; di++){
	            x.push_back(data.at(di));
		        y.push_back(data.at(data.size()/2+di));
	            }
			for(int l=0; l<x_r.size(); l++){
				for(int m=0; m<y_r.size(); m++){
					xp = x_r.at(l);
					yp = y_r.at(m);
				    flag = check_safety(x, y, xp, yp, l2, r2);
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
	std::vector<double> x;
	std::vector<double> y;
	clock_t begin, end;
	begin = clock();


	std::vector<float> delta;

	delta = rangeByStep(-(1.0/4)*M_PI, (1.0/4)*M_PI, 0.1);
	offline_tube_generate(delta);
	unordered_map<string, vector<float>> map;
	map = load_table();
    
	float xp = 1;
    float yp = 2;
    float yaw = 0.3;
	Point l2;
	Point r2;
    l2.x = -1;
	l2.y = 1;
    r2.x = 1;
	r2.y = -1;
	begin = clock();
	online_check_for_single_grid(map, xp, yp, yaw, delta, l2, r2);
	end = clock();
	printf("single_point check time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

    float x_min = -0.5;
    float x_max = 0.5;
    float y_min = -0.5;
    float y_max = 0.5;
    float yaw_min = -M_PI/2.0;
    float yaw_max = M_PI/2.0;
    
	begin = clock();
    online_check_for_interval(map, x_min, x_max, y_min, y_max, yaw_min, y_max, delta, l2, r2);
	end = clock();
	printf("interval check time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	return 0;
}
	