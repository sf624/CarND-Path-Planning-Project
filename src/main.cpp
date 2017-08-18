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
#include "opt.h"

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

  h.onMessage([&opt,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            static int lane = 1;
            static double ref_vel = 0;
            int prev_size = previous_path_x.size();
            //prev_size = prev_size < 2 ? prev_size : 2;

            // Sensor fusion START
            if(prev_size>0){
              car_s = end_path_s;
            }

            bool too_close = false;

            //find ref_v to use
            double leading_car_speed = 50;
            for(int i=0; i<sensor_fusion.size(); i++){
              float d = sensor_fusion[i][6];
              if(d < (2+4*lane+2) && d > (2+4*lane-2)){
                double vx = sensor_fusion[i][3];
                double vy = sensor_fusion[i][4];
                double check_speed = sqrt(vx*vx+vy*vy);
                double check_car_s = sensor_fusion[i][5];

                check_car_s += (double)prev_size * 0.02 * check_speed;

                if((check_car_s>car_s) && ((check_car_s-car_s)<30)){
                  leading_car_speed = check_speed * 2.24;
                  too_close = true;
                }
              }
            }

            if(too_close && ref_vel > leading_car_speed - 5){
              ref_vel -= 0.224;
            }
            else if(ref_vel < 49.5){
              ref_vel += 0.224;
            }

            // Sensor fusion END


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

            /*
             * Generate target-lane spline in reference frame, using 5 closest waypoint.
            */

            vector<vector<double>> x_candidates;
            vector<vector<double>> y_candidates;

            for(int lane=0; lane<3; lane++){
              tk::spline lane_spline;
              vector<double> ptsx_lane;
              vector<double> ptsy_lane;

              int closestWaypoint = ClosestWaypoint(ref_x, ref_y, map_waypoints_x, map_waypoints_y);
              for(int i=0; i<3; i++){
                int idx = (closestWaypoint + i - 1 + map_waypoints_s.size()) % map_waypoints_s.size();
                ptsx_lane.push_back(map_waypoints_x[idx] + (2+4*lane) * map_waypoints_dx[idx]);
                ptsy_lane.push_back(map_waypoints_y[idx] + (2+4*lane) * map_waypoints_dy[idx]);
                //cout << idx << "," << ptsx_lane[i] << endl;
              }

              for(int i=0; i<ptsx_lane.size(); i++){
                double shift_x = ptsx_lane[i] - ref_x;
                double shift_y = ptsy_lane[i] - ref_y;

                ptsx_lane[i] = (shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw));
                ptsy_lane[i] = (shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw));
                //cout << i << "," << ptsx_lane[i] << endl;
              }

              lane_spline.set_points(ptsx_lane, ptsy_lane);

              /*
               * Generate target-trajectory spline in reference frame, using target-lane spline.
               */

              // Intend to reach target-lane in 2 seconds.
              vector<double> temp_ptsx;
              vector<double> temp_ptsy;
              copy(ptsx.begin(),ptsx.end(),back_inserter(temp_ptsx));
              copy(ptsy.begin(),ptsy.end(),back_inserter(temp_ptsy));

              for(int i=0; i<3; i++){
                double speed = car_speed < 25 ? 25 : car_speed;
                double x = 2*(i+1) * speed * 0.44704;
                temp_ptsx.push_back(x);
                temp_ptsy.push_back(lane_spline(x));
              }

              tk::spline trajectory_spline;
              trajectory_spline.set_points(temp_ptsx, temp_ptsy);

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
              //double target_speed = i == lane ? ref_vel : (ref_vel > 49.5 ? ref_vel : ref_vel + 0.224);

              for(int i=0; i<250-prev_size; i++){
                double target_speed = i == lane ? ref_vel : (ref_vel > 49.5 ? ref_vel : ref_vel + 0.224);
                double N = (target_dist/(0.02*target_speed/2.24));
                double x_point = x_add_on + target_x/N;
                double y_point = trajectory_spline(x_point);

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

            /*
            * If too_close, attempt other trajectory.
            * In order to proceed, other trajectories are evaluated.
            */

            double min_cost = 1.0e19;
            int min_cost_idx = 0;
            for(int k=1; k<3; k++){
              int i = (lane+k)%3;
              double cost = 0;
              for(int t=0; t<x_candidates[i].size(); t++){
                double x = x_candidates[i][t];
                double y = y_candidates[i][t];
                vector<double> frenet = getFrenet(x,y,ref_yaw,map_waypoints_x,map_waypoints_y);
                double car_s = frenet[0];
                double car_d = frenet[1];

                // Penalize collision
                for(int i=0; i<sensor_fusion.size(); i++){
                  float d = sensor_fusion[i][6];
                  if(abs(car_d - d) < 2){
                    double check_x = sensor_fusion[i][1];
                    double check_y = sensor_fusion[i][2];
                    double check_vx = sensor_fusion[i][3];
                    double check_vy = sensor_fusion[i][4];
                    //double check_speed = sqrt(vx*vx+vy*vy);
                    //double check_car_s = sensor_fusion[i][5];

                    check_x += t * 0.02 * check_vx;
                    check_y += t * 0.02 * check_vy;

                    if(distance(x,y,check_x,check_y)<10){
                      cout << "Car is close." << endl;
                      cost += 1.0e5;
                    }
                  }
                }

                // Penalize lane change
                cost +=  pow((car_d - (2+4*i)),2);

                // Reward large velocity.
                //cost -= car_s;

                //if(i==lane && too_close) cost += 100;
              }

              if(cost < min_cost){
                min_cost = cost;
                min_cost_idx = i;
              }
            }

            if(too_close && min_cost < 1.0e4 && abs(car_d - lane*4 - 2) < 0.5){
              lane = min_cost_idx;
            }

            cout << "Target Lane: " << lane << " car_d = " << car_d << endl;
            cout << "====================================" << endl;

            json msgJson;

            x_candidates[lane].erase(x_candidates[lane].begin()+50,x_candidates[lane].end());
            y_candidates[lane].erase(y_candidates[lane].begin()+50,y_candidates[lane].end());

            msgJson["next_x"] = x_candidates[lane];
          	msgJson["next_y"] = y_candidates[lane];

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
