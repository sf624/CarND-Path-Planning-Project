#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);
	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];
	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// Distance in s-coordinate with sign.
// If sign is +, a is front of b.
double s_distance(double a, double b, double max_s){
  double ret = a - b;
  while(ret < -max_s/2.0) ret += max_s;
  while(ret >  max_s/2.0) ret -= max_s;
  return ret;
}

int main() {
  uWS::Hub h;

  OPT opt;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&max_s,&opt,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];

          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
            // Data format: [id,x,y,vx,vy,s,d]
          	auto sensor_fusion = j[1]["sensor_fusion"];


          	// Static variables used.
            static int target_lane = 1;               // Current target lane to follow. 0 ... left, 1 ... middle, 2 ... right
            static double ref_vel = 0;                // Current reference velocity in mph.
            static bool is_in_lane_change = false;    // Whether or not the car is in lane changing mode or not.

            // Extract first 50 points from previous unused trajectory.
            // Which means trajectory over 1 seconds are refined.
            int prev_size = previous_path_x.size() < 50 ? previous_path_x.size() : 50;
            if(prev_size > 0){
              vector<double> sd = getFrenet(previous_path_x[prev_size], previous_path_y[prev_size],
                                  deg2rad(car_yaw), map_waypoints_x, map_waypoints_y);
              end_path_s = sd[0];
              end_path_d = sd[1];
            }


            // Determine reference velocity to generate trajectory.
            double leading_car_speed = 50;
            double min_distance = 1e9;
            bool emergecy_brake = false;
            for(int i=0; i<sensor_fusion.size(); i++){
              float d = sensor_fusion[i][6];
              if(abs(end_path_d - d)<3){
                double vx = sensor_fusion[i][3];
                double vy = sensor_fusion[i][4];
                double check_speed = sqrt(vx*vx+vy*vy);
                double check_car_s = sensor_fusion[i][5];

                if(s_distance(check_car_s,car_s,max_s)>0 && s_distance(check_car_s,end_path_s,max_s)<0){
                  emergecy_brake = true;
                  break;
                }

                check_car_s += check_speed * 0.02 * prev_size;
                if(s_distance(check_car_s,end_path_s,max_s)>0){
                  if(s_distance(check_car_s,end_path_s,max_s) < min_distance){
                    min_distance = s_distance(check_car_s,end_path_s,max_s);
                    leading_car_speed = check_speed * 2.24;
                  }
                }
              }
            }


            if(emergecy_brake){
              ref_vel -= 0.224 * 1.5;
            }
            else if( (ref_vel > leading_car_speed
              && pow(ref_vel/2.24 - leading_car_speed/2.24,2)/10.0 > (min_distance - leading_car_speed / 2.24 * 2.0))
              || pow(ref_vel/2.24 - leading_car_speed/2.24,2)/10.0 < (leading_car_speed / 2.24 * 2.0 - min_distance) ){
              ref_vel -= 0.224;
              //cout << min_distance << endl;
            }
            else if(ref_vel < 49.5){
              ref_vel += 0.224;
            }



            // Extract last two points of unused trajectory and
            // use them to define reference frame where another new trajecty is generated.
            vector<double> ptsx;
            vector<double> ptsy;

            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);

            if(prev_size < 2){
              double prev_car_x = car_x - cos(car_yaw);
              double prev_car_y = car_y - sin(car_yaw);

              ptsx.push_back(prev_car_x);
              ptsx.push_back(car_x);

              ptsy.push_back(prev_car_y);
              ptsy.push_back(car_y);
            }
            else{
              ref_x = previous_path_x[prev_size-1];
              ref_y = previous_path_y[prev_size-1];

              double ref_x_prev = previous_path_x[prev_size-2];
              double ref_y_prev = previous_path_y[prev_size-2];
              ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

              ptsx.push_back(ref_x_prev);
              ptsx.push_back(ref_x);

              ptsy.push_back(ref_y_prev);
              ptsy.push_back(ref_y);
            }

            for(int i=0; i<ptsx.size(); i++){
              double shift_x = ptsx[i] - ref_x;
              double shift_y = ptsy[i] - ref_y;

              ptsx[i] = (shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw));
              ptsy[i] = (shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw));
            }



            // Generate three trajectory candidates which each corresponds to three differenct target lane.
            vector<vector<double>> x_candidates;
            vector<vector<double>> y_candidates;

            for(int lane=0; lane<3; lane++){
              vector<double> ptsx_lane;
              vector<double> ptsy_lane;

              // Extract 3 nearest waypoints and generate spline interpolated lanes' points.
              int closestWaypoint = ClosestWaypoint(ref_x, ref_y, map_waypoints_x, map_waypoints_y);
              for(int i=0; i<3; i++){
                int idx = (closestWaypoint + i - 1 + map_waypoints_s.size()) % map_waypoints_s.size();
                ptsx_lane.push_back(map_waypoints_x[idx] + (2+4*lane) * map_waypoints_dx[idx]);
                ptsy_lane.push_back(map_waypoints_y[idx] + (2+4*lane) * map_waypoints_dy[idx]);
              }

              for(int i=0; i<ptsx_lane.size(); i++){
                double shift_x = ptsx_lane[i] - ref_x;
                double shift_y = ptsy_lane[i] - ref_y;

                ptsx_lane[i] = (shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw));
                ptsy_lane[i] = (shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw));
              }

              tk::spline lane_spline;
              lane_spline.set_points(ptsx_lane, ptsy_lane);


              // Extract another 3 waypoints forward to ego car to generate
              // roughly spline interploted lanes' point.
              int temp = ptsx_lane.size();

              for(int i=3; i<6; i++){
                int idx = (closestWaypoint + i - 1 + map_waypoints_s.size()) % map_waypoints_s.size();
                ptsx_lane.push_back(map_waypoints_x[idx] + (2+4*lane) * map_waypoints_dx[idx]);
                ptsy_lane.push_back(map_waypoints_y[idx] + (2+4*lane) * map_waypoints_dy[idx]);
              }

              for(int i=temp; i<ptsx_lane.size(); i++){
                double shift_x = ptsx_lane[i] - ref_x;
                double shift_y = ptsy_lane[i] - ref_y;

                ptsx_lane[i] = (shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw));
                ptsy_lane[i] = (shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw));
              }

              tk::spline rough_lane_spline;
              rough_lane_spline.set_points(ptsx_lane, ptsy_lane);

              // Generate spline interpolated target trajectry, which ego vehicle should follow.
              // Tre trajectory intend to let ego car reach the target lane in 2 seconds when over 25m/s.
              vector<double> temp_ptsx;
              vector<double> temp_ptsy;
              copy(ptsx.begin(),ptsx.end(),back_inserter(temp_ptsx));
              copy(ptsy.begin(),ptsy.end(),back_inserter(temp_ptsy));

              for(int i=0; i<5; i++){
                double speed = car_speed < 25 ? 25 : car_speed;
                double x = 2*(i+1) * speed * 0.44704;
                if(lane == target_lane && !is_in_lane_change) x /= 2;
                temp_ptsx.push_back(x);
                temp_ptsy.push_back(lane_spline(x));
              }

              tk::spline trajectory_spline;
              trajectory_spline.set_points(temp_ptsx, temp_ptsy);

              // Generate candidate target trajectory in simulator's x-y frame.
              vector<double> next_x_vals;
              vector<double> next_y_vals;

              for(int i=0; i<prev_size; i++){
                next_x_vals.push_back(previous_path_x[i]);
                next_y_vals.push_back(previous_path_y[i]);
              }

              double target_x = 30.0;
              double target_y = trajectory_spline(target_x);
              double target_dist = sqrt(pow(target_x,2)+pow(target_y,2));

              double x_add_on = 0;

              for(int i=0; i<250-prev_size; i++){
                double N = (target_dist/(0.02*ref_vel/2.24));
                double x_point = x_add_on + target_x/N;
                double y_point = trajectory_spline(x_point);
                if(i >= 100 - prev_size){
                  y_point = rough_lane_spline(x_point);
                }

                x_add_on = x_point;

                double x_ref = x_point;
                double y_ref = y_point;

                x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
                y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

                x_point += ref_x;
                y_point += ref_y;

                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
              }

              x_candidates.push_back(next_x_vals);
              y_candidates.push_back(next_y_vals);
            }

            // Evaluate lanes.
            // Lane cost below stores general cost and special collision cost.
            vector<double> lane_cost{0,0,0};
            vector<double> lane_collision_cost{0,0,0};
            // Penalise lane which has distance between target lane.
            for(int lane=0; lane<lane_cost.size(); lane++){
              lane_cost[lane] += exp(-100) * abs(target_lane-lane);
            }
            // Penalise lane which has objective cars in front of ego car.
            for(int i=0; i<sensor_fusion.size(); i++){
              double d = sensor_fusion[i][6];
              int lane = (int)round((d-2.0)/4.0);
              lane = lane < 0 ? 0 : lane > 2 ? 2 : lane;

              double vx = sensor_fusion[i][3];
              double vy = sensor_fusion[i][4];
              double check_speed = sqrt(vx*vx+vy*vy);
              double check_car_s = sensor_fusion[i][5];

              if(s_distance(check_car_s,car_s,max_s)>0){
                lane_cost[lane] += exp(-abs(s_distance(check_car_s,car_s,max_s)));
              }
            }
            // Penalise lane which where collision is estimated.
            for(int lane=0; lane<lane_cost.size(); lane++){
              for(int t=0; t<x_candidates[lane].size(); t++){
                double x = x_candidates[lane][t];
                double y = y_candidates[lane][t];
                vector<double> frenet = getFrenet(x,y,ref_yaw,map_waypoints_x,map_waypoints_y);
                double car_s = frenet[0];
                double car_d = frenet[1];
                // Penalize collision
                for(int i=0; i<sensor_fusion.size(); i++){
                  float d = sensor_fusion[i][6];
                  if(abs(car_d - d) < 3){
                    double check_x = sensor_fusion[i][1];
                    double check_y = sensor_fusion[i][2];
                    double check_vx = sensor_fusion[i][3];
                    double check_vy = sensor_fusion[i][4];
                    double check_speed = sqrt(check_vx*check_vx+check_vy*check_vy);
                    double check_car_s = sensor_fusion[i][5];
                    check_car_s += check_speed * 0.02 * t;

                    if(abs(s_distance(car_s,check_car_s,max_s))<max(ref_vel/2.24,check_speed)*1){
                      lane_collision_cost[lane] += 1.0;
                    }
                  }
                }
              }
            }


            // Finite State Machine.
            // If the ego car is in lane changing mode, just look whether the ego car reached the target lane.
            // If not in lane chanding mode and collision-free lane with smallest lanes' cost
            // differs from current target lane, lane changing mode is exexuted.
            if(is_in_lane_change){
              // do nothing, but switch is_in_lane_change if accomplished.
              if(abs(car_d - (target_lane*4+2)) < 0.5){
                is_in_lane_change = false;
                cout << "DONE!" << endl;
              }
            }
            else{
              int argmin_lane = distance(lane_cost.begin(), min_element(lane_cost.begin(),lane_cost.end()) );
              // Note: Collision free lane collision cost must be smaller than 1.0
              if(argmin_lane != target_lane && lane_collision_cost[argmin_lane] < 1.0){
                is_in_lane_change = true;
                cout << "------------------------------" << endl;
                for(int lane=0;lane<3;lane++){
                  if(lane == target_lane) cout << "(now) ";
                  else if (lane == argmin_lane) cout << "(min) ";
                  else cout << "      ";
                  cout << "lane_cost[" << lane << "] = " << lane_cost[lane] << endl;
                }
                cout << endl;
                target_lane = argmin_lane;
                cout << "No collision estimated." << endl;
                cout << "Executing lane change ... " << flush;
              }
            }

            json msgJson;

            // Only first 50 target trajectory points are sent to simulator.
            // Which lasts for 1 seconds.
            //vector<double> next_x;
            //vector<double> next_y;
            //copy(x_candidates[target_lane].begin(),x_candidates[target_lane].begin()+100,back_inserter(next_x));
            //copy(y_candidates[target_lane].begin(),y_candidates[target_lane].begin()+100,back_inserter(next_y));

            msgJson["next_x"] = x_candidates[target_lane];
          	msgJson["next_y"] = y_candidates[target_lane];

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
